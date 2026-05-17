"""Bandwidth selection for kernel density estimation."""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from ._kernels import kd_K_fun, kd_bw_fun


_KERNEL_ALIASES = {
    "epa": "epa", "epanechnikov": "epa",
    "uni": "uni", "uniform": "uni",
}
_KERNEL_TYPE_NAMES = {"epa": "Epanechnikov", "uni": "Uniform"}


@dataclass
class KdbwselectResult:
    """Container for :func:`kdbwselect` results.

    Attributes
    ----------
    bws : pandas.DataFrame
        Bandwidths over the evaluation grid. When ``bwselect != "all"`` the
        frame has shape ``(neval, 3)`` with columns ``eval, h, b``. When
        ``bwselect == "all"`` it widens to shape ``(neval, 7)`` with columns
        ``eval, h.mse.dpi, b.mse.dpi, h.ce.dpi, b.ce.dpi, h.ce.rot,
        b.ce.rot``.
    bws_imse : pandas.DataFrame or None
        Only populated when ``bwselect == "all"``: a ``(2, 2)`` table indexed
        ``["h", "b"]`` with columns ``["IMSE-DPI", "IMSE-ROT"]`` carrying
        the IMSE bandwidths. ``None`` otherwise.
    opt : dict
        Options actually used. Keys: ``p, n, neval, kernel, bwselect``.
    call : dict
        Optional dict for storing the call metadata (default empty).
    """
    bws: pd.DataFrame
    bws_imse: Optional[pd.DataFrame]
    opt: Dict[str, Any]
    call: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"<KdbwselectResult n={self.opt['n']} kernel={self.opt['kernel']} "
                f"bwselect={self.opt['bwselect']} neval={self.opt['neval']}>")

    def tidy(self) -> pd.DataFrame: return self.bws.copy()

    def glance(self) -> pd.DataFrame:
        return pd.DataFrame([{
            "n": self.opt["n"], "neval": self.opt["neval"],
            "kernel": self.opt["kernel"],
            "bwselect": self.opt["bwselect"],
        }])


def _normalize_kernel(kernel):
    k = _KERNEL_ALIASES.get(kernel.lower())
    if k is None:
        raise ValueError(
            "kernel incorrectly specified. Supported kernels for kdbwselect: epa, uni."
        )
    return k


def kdbwselect(
    x,
    eval=None, neval=None,
    kernel: str = "epa",
    bwselect: str = "mse-dpi",
    bwcheck: Optional[int] = 21,
    imsegrid: int = 30,
    subset=None,
    data=None,
) -> KdbwselectResult:
    """Bandwidth selectors for kernel density estimation.

    Implements bandwidth selectors for kernel density point estimators and
    inference procedures developed in Calonico, Cattaneo and Farrell (2018),
    with related optimality results in Calonico, Cattaneo and Farrell (2022).
    It also implements other bandwidth selectors available in the literature.
    Companion command: :func:`kdrobust`.

    If ``data`` is a pandas DataFrame, ``x`` and ``subset`` may be passed as
    column-name strings.

    Parameters
    ----------
    x : array-like or column-name string
        Independent variable.
    eval : array-like or None
        Evaluation points. If ``None``, 30 quantile-spaced points (deciles
        0.1 through 0.9 in equal steps) over the support of ``x`` are used.
    neval : int or None
        Number of quantile-spaced evaluation points. Default 30.
    kernel : {"epa", "uni"}
        Kernel function. Default ``"epa"``.
    bwselect : str
        Bandwidth selector. Options are ``"mse-dpi"``, ``"imse-dpi"``,
        ``"imse-rot"``, ``"ce-dpi"``, ``"ce-rot"``, and ``"all"``.
        Default ``"mse-dpi"``.
    bwcheck : int or None
        If positive, selected bandwidths are enlarged so at least
        ``bwcheck`` observations are available at each evaluation point.
        Default 21.
    imsegrid : int
        Number of evaluation points used to compute IMSE bandwidth selectors.
        Default 30.
    subset : array-like or column-name string, optional
        Boolean rule selecting observations to use.
    data : pandas.DataFrame, optional
        Data source for column-name inputs.

    Returns
    -------
    result : KdbwselectResult
        An object with the following attributes:

        - ``bws`` (``pandas.DataFrame``): bandwidths over the evaluation
          grid. Shape ``(neval, 3)`` with columns ``eval, h, b`` when
          ``bwselect != "all"``; widens to shape ``(neval, 7)`` with
          columns ``eval, h.mse.dpi, b.mse.dpi, h.ce.dpi, b.ce.dpi,
          h.ce.rot, b.ce.rot`` when ``bwselect == "all"``.
        - ``bws_imse`` (``pandas.DataFrame`` or ``None``): only populated
          when ``bwselect == "all"`` - a ``(2, 2)`` table indexed
          ``["h", "b"]`` with columns ``["IMSE-DPI", "IMSE-ROT"]``
          carrying the IMSE bandwidths.
        - ``opt`` (``dict``): options actually used. Keys: ``p, n, neval,
          kernel, bwselect``.
    """
    if data is not None:
        from ._helpers import resolve_name
        x      = resolve_name(x,      data, "x")
        subset = resolve_name(subset, data, "subset")

    x = np.asarray(x, dtype=float).ravel()
    p = 2
    deriv = 0

    kernel = _normalize_kernel(kernel)
    bwselect = bwselect.lower()

    if subset is not None:
        idx = np.asarray(subset).ravel()
        if idx.dtype != bool:
            raise TypeError("subset must be a boolean array.")
        x = x[idx]
    x = x[np.isfinite(x)]

    N = len(x)
    if bwcheck is not None and bwcheck > N:
        warnings.warn(f"bwcheck ({bwcheck}) > N ({N}); reducing bwcheck to N.")
        bwcheck = N

    if eval is None:
        if neval is None:
            eval = np.quantile(x, np.linspace(0.1, 0.9, 30))
        else:
            eval = np.quantile(x, np.linspace(0.1, 0.9, int(neval)))
    eval = np.asarray(eval, dtype=float).ravel()
    neval = len(eval)

    # Kernel-specific constants
    if kernel == "epa":
        C_h, C_b = 2.34, 3.49
    else:
        C_h, C_b = 1.06, 1.0

    if bwselect == "all":
        bws = np.full((neval, 6), np.nan)
        colnames = ["h.mse.dpi", "b.mse.dpi", "h.ce.dpi", "b.ce.dpi",
                    "h.ce.rot", "b.ce.rot"]
        bws_imse = np.full((2, 2), np.nan)
    else:
        bws = np.full((neval, 2), np.nan)
        colnames = ["h", "b"]
        bws_imse = None

    # IMSE-ROT
    h_imse_rot = np.std(x, ddof=1) * C_h * N ** (-1.0 / (1 + 2 * p))
    b_imse_rot = np.std(x, ddof=1) * C_b * N ** (-1.0 / (1 + 2 * (p + 2) + 2 * p))

    if bwselect == "imse-rot":
        bws[:, 0] = h_imse_rot; bws[:, 1] = b_imse_rot

    h_imse_dpi = None
    if bwselect in ("imse-dpi", "all"):
        qseq_imse = np.linspace(0.1, 0.9, imsegrid)
        eval_imse = np.quantile(x, qseq_imse)
        B_h = np.empty(imsegrid); V_h = np.empty(imsegrid)
        for i in range(imsegrid):
            K_b_tup = kd_K_fun((x - eval_imse[i]) / b_imse_rot, v=p + 2, r=p, kernel=kernel)
            K_h_tup = kd_K_fun((x - eval_imse[i]) / h_imse_rot, v=p,     r=deriv, kernel=kernel)
            f_b = np.mean(K_b_tup[0]) / b_imse_rot ** (1 + p)
            f_h_rot = np.mean(K_h_tup[0]) / h_imse_rot
            B_h[i] = f_b * K_h_tup[1]
            V_h[i] = f_h_rot * K_h_tup[2]
        h_imse_dpi = kd_bw_fun(np.mean(V_h), np.mean(B_h), N, v=p, r=deriv)
        if bwselect == "imse-dpi":
            bws[:, 0] = h_imse_dpi
            bws[:, 1] = b_imse_rot

    if bwselect == "all":
        bws_imse[:, 0] = [h_imse_dpi, b_imse_rot]
        bws_imse[:, 1] = [h_imse_rot, b_imse_rot]

    if bwselect in ("all", "mse-dpi", "ce-dpi", "ce-rot"):
        for i in range(neval):
            if bwcheck is not None:
                bw_min = np.sort(np.abs(x - eval[i]))[bwcheck - 1]
            else:
                bw_min = 0.0

            K_b_tup = kd_K_fun((x - eval[i]) / b_imse_rot, v=p + 2, r=p, kernel=kernel)
            K_h_tup = kd_K_fun((x - eval[i]) / h_imse_rot, v=p,     r=deriv, kernel=kernel)
            f_b = np.mean(K_b_tup[0]) / b_imse_rot ** (1 + p)
            f_h_rot = np.mean(K_h_tup[0]) / h_imse_rot
            Bh = f_b * K_h_tup[1]
            Vh = f_h_rot * K_h_tup[2]
            h_mse_dpi = kd_bw_fun(Vh, Bh, N, v=p, r=deriv)
            b_mse_dpi = b_imse_rot

            if bwcheck is not None:
                h_mse_dpi = max(h_mse_dpi, bw_min)
                b_mse_dpi = max(b_mse_dpi, bw_min)

            h_ce_rot = b_ce_rot = h_ce_dpi = b_ce_dpi = np.nan

            if bwselect in ("ce-rot", "all"):
                h_ce_rot = h_mse_dpi * N ** (-(p - 2) / ((1 + 2 * p) * (1 + p + 2)))
                b_ce_rot = b_mse_dpi * N ** (-(p - 2) / ((1 + 2 * p) * (1 + p + 2)))
                if bwcheck is not None:
                    h_ce_rot = max(h_ce_rot, bw_min)
                    b_ce_rot = max(b_ce_rot, bw_min)

            if bwselect in ("ce-dpi", "all"):
                h_ce_dpi = _kd_cer(x, eval[i], h_mse_dpi, h_mse_dpi, p, kernel)
                b_ce_dpi = b_mse_dpi
                if bwcheck is not None:
                    h_ce_dpi = max(h_ce_dpi, bw_min)
                    b_ce_dpi = max(b_ce_dpi, bw_min)

            if bwselect == "all":
                bws[i, :] = [h_mse_dpi, b_mse_dpi, h_ce_dpi, b_ce_dpi,
                             h_ce_rot, b_ce_rot]
            elif bwselect == "mse-dpi":
                bws[i, :2] = [h_mse_dpi, b_mse_dpi]
            elif bwselect == "ce-rot":
                bws[i, :2] = [h_ce_rot, b_ce_rot]
            elif bwselect == "ce-dpi":
                bws[i, :2] = [h_ce_dpi, b_ce_dpi]

    bws_df = pd.DataFrame(np.column_stack([eval, bws]), columns=["eval", *colnames])
    bws_imse_df = None
    if bws_imse is not None:
        bws_imse_df = pd.DataFrame(bws_imse, columns=["IMSE-DPI", "IMSE-ROT"],
                                   index=["h", "b"])

    opt = dict(p=p, n=N, neval=neval,
               kernel=_KERNEL_TYPE_NAMES[kernel], bwselect=bwselect)
    return KdbwselectResult(bws=bws_df, bws_imse=bws_imse_df, opt=opt)


def _kd_cer(x, x0, h, b, v, kernel):
    """Coverage-error-optimal density bandwidth (port of kd.cer.fun)."""
    from scipy import integrate, optimize
    from scipy.stats import norm

    n = len(x)
    rho = h / b
    rng = float(np.max(x) - np.min(x))

    q_rot = np.std(x, ddof=1) * n ** (-1.0 / (1 + 2 * v + 2 * (v + 2)))
    K_q = kd_K_fun((x - x0) / q_rot, v=2, r=v + 2, kernel="gau")
    f_r_2 = np.mean(K_q[0]) / q_rot ** (1 + 2 * v)

    v_K = kd_K_fun(np.array([1.0]), v=v, r=0, kernel=kernel)[1]

    def M_fun(u):
        Ku = kd_K_fun(u, v=v, r=0, kernel=kernel)[0]
        Lu = kd_K_fun(rho * u, v=v + 2, r=v, kernel=kernel)[0]
        return Ku - rho ** (1 + v) * Lu * v_K

    def K_fun(u):
        return kd_K_fun(u, v=v, r=0, kernel=kernel)[0]

    def L_fun(u):
        return kd_K_fun(u, v=v + 2, r=v, kernel=kernel)[0]

    def v_int(power, f):
        return integrate.quad(lambda u: f(np.array([u]))[0] ** power if False
                              else f(u) ** power, -np.inf, np.inf)[0]

    def m_int(m, f):
        return integrate.quad(lambda u: u ** m * f(u), -np.inf, np.inf)[0]

    v_M_2 = v_int(2, M_fun)
    v_M_3 = v_int(3, M_fun)
    v_M_4 = v_int(4, M_fun)
    m_K_4 = m_int(4, K_fun)
    m_K_2 = m_int(2, K_fun)
    m_L_2 = m_int(2, L_fun)

    z = norm.ppf(0.975)
    q1 = v_M_4 * (z ** 2 - 3) / 6 - v_M_3 ** 2 * (z ** 4 - 4 * z ** 2 + 15) / 9
    q2 = f_r_2 ** 2 * (m_K_4 - rho ** (-2) * m_K_2 * m_L_2 / 12) ** 2 * v_M_2
    q3 = f_r_2     * (m_K_4 - rho ** (-2) * m_K_2 * m_L_2 / 12) * v_M_3 * (2 * z ** 2) / 3

    def obj(H):
        return (H ** (-1) * q1
                - H ** (1 + 2 * (v + 2)) * q2
                + H ** (v + 2) * q3) ** 2

    res = optimize.minimize_scalar(obj, bounds=(np.finfo(float).eps, rng), method="bounded", options={"xatol": 1e-10})
    return float(res.x) * n ** (-1.0 / (v + 3))
