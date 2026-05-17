"""Kernel density with robust bias-corrected inference."""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from ._kernels import kd_K_fun


_KERNEL_ALIASES = {
    "epa": "epa", "epanechnikov": "epa",
    "uni": "uni", "uniform": "uni",
}
_KERNEL_TYPE_NAMES = {"epa": "Epanechnikov", "uni": "Uniform"}


@dataclass
class KdrobustResult:
    """Container for :func:`kdrobust` results.

    Attributes
    ----------
    Estimate : pandas.DataFrame
        Density estimates and standard errors over the evaluation grid, of
        shape ``(neval, 8)``. Columns: ``eval, h, b, N, tau.us, tau.bc,
        se.us, se.rb``.
    opt : dict
        Options actually used. Keys: ``p, kernel, n, neval, bwselect``.
    call : dict
        Optional dict for storing the call metadata (default empty).
    """
    Estimate: pd.DataFrame
    opt: Dict[str, Any]
    call: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"<KdrobustResult n={self.opt['n']} kernel={self.opt['kernel']} "
                f"bwselect={self.opt['bwselect']} neval={self.opt['neval']}>")

    def summary(self, alpha: float = 0.05) -> pd.DataFrame:
        z = norm.ppf(1 - alpha / 2)
        df = self.Estimate.copy()
        df["CI.lower"] = df["tau.bc"] - z * df["se.rb"]
        df["CI.upper"] = df["tau.bc"] + z * df["se.rb"]
        return df[["eval", "h", "b", "N", "tau.us", "se.us",
                   "tau.bc", "se.rb", "CI.lower", "CI.upper"]]

    def tidy(self, conf_level: float = 0.95) -> pd.DataFrame:
        z = norm.ppf(1 - (1 - conf_level) / 2)
        e = self.Estimate
        return pd.DataFrame({
            "eval":      e["eval"].values,
            "estimate":  e["tau.us"].values,
            "std.error": e["se.us"].values,
            "tau.bc":    e["tau.bc"].values,
            "se.rb":     e["se.rb"].values,
            "h":         e["h"].values,
            "b":         e["b"].values,
            "n.eff":     e["N"].values,
            "conf.low":  e["tau.bc"].values - z * e["se.rb"].values,
            "conf.high": e["tau.bc"].values + z * e["se.rb"].values,
        })

    def glance(self) -> pd.DataFrame:
        return pd.DataFrame([{
            "n": self.opt["n"], "neval": self.opt["neval"],
            "kernel": self.opt["kernel"],
            "bwselect": self.opt["bwselect"],
        }])


def _normalize_kernel(kernel: str) -> str:
    k = _KERNEL_ALIASES.get(kernel.lower())
    if k is None:
        raise ValueError(
            "kernel incorrectly specified. Supported kernels for kdrobust: epa, uni."
        )
    return k


def kdrobust(
    x,
    eval=None, neval=None,
    h=None, b=None, rho: float = 1.0,
    kernel: str = "epa",
    bwselect: Optional[str] = None,
    bwcheck: Optional[int] = 21,
    imsegrid: int = 30,
    level: float = 95,
    subset=None,
    data=None,
) -> KdrobustResult:
    """Kernel density estimator with robust bias-corrected inference.

    Implements kernel density point estimators, with robust bias-corrected
    confidence intervals and inference procedures developed in Calonico,
    Cattaneo and Farrell (2018), with related optimality results in Calonico,
    Cattaneo and Farrell (2022). It also implements other estimation and
    inference procedures available in the literature. Companion command:
    :func:`kdbwselect`.

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
    h, b : float, array, or None
        Main bandwidth and bias-correction bandwidth. Each can be scalar or
        have the same length as ``eval``. If not supplied, bandwidths are
        selected by :func:`kdbwselect`.
    rho : float
        Sets ``b = h / rho`` when ``b`` is not supplied. Default 1.
    kernel : {"epa", "uni"}
        Kernel function. Default ``"epa"``.
    bwselect : str or None
        Bandwidth selector. Options are ``"mse-dpi"``, ``"imse-dpi"``,
        ``"imse-rot"``, ``"ce-dpi"``, and ``"ce-rot"``. Use
        :func:`kdbwselect` with ``bwselect="all"`` to report all available
        bandwidth selectors. Default is ``"mse-dpi"`` with one evaluation
        point and ``"imse-dpi"`` otherwise.
    bwcheck : int or None
        If positive, selected bandwidths are enlarged so at least
        ``bwcheck`` observations are available at each evaluation point.
        Default 21.
    imsegrid : int
        Number of evaluation points used to compute IMSE bandwidth selectors.
        Default 30.
    level : float
        Confidence level in percent. Default 95.
    subset : array-like or column-name string, optional
        Boolean rule selecting observations to use.
    data : pandas.DataFrame, optional
        Data source for column-name inputs.

    Returns
    -------
    result : KdrobustResult
        An object with the following attributes:

        - ``Estimate`` (``pandas.DataFrame``, shape ``(neval, 8)``):
          density estimates and standard errors over the evaluation grid.
          Columns: ``eval, h, b, N, tau.us, tau.bc, se.us, se.rb``.
        - ``opt`` (``dict``): options actually used. Keys: ``p, kernel, n,
          neval, bwselect``.
    """
    if data is not None:
        from ._helpers import resolve_name
        x      = resolve_name(x,      data, "x")
        subset = resolve_name(subset, data, "subset")

    x = np.asarray(x, dtype=float).ravel()
    p = 2
    deriv = 0

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

    if h is None and bwselect is None:
        bwselect = "mse-dpi" if neval == 1 else "imse-dpi"

    kernel = _normalize_kernel(kernel)
    if bwselect is not None: bwselect = bwselect.lower()

    if level <= 0 or level >= 100:
        raise ValueError("level must be in (0, 100).")
    if rho < 0:
        raise ValueError("rho must be non-negative.")

    kernel_type = _KERNEL_TYPE_NAMES[kernel]

    if h is not None:
        h = np.asarray(h, dtype=float).ravel()
        if b is None and rho > 0:
            b = h / rho
        if b is None and rho == 0:
            raise ValueError("When h is provided and rho=0, b must also be provided.")
        b = np.asarray(b, dtype=float).ravel()
        bwselect = "Manual"
    else:
        from .kdbwselect import kdbwselect
        r = kdbwselect(x, eval=eval, bwselect=bwselect, bwcheck=bwcheck,
                       imsegrid=imsegrid, kernel=kernel)
        h = r.bws.iloc[:, 1].to_numpy()
        b = r.bws.iloc[:, 2].to_numpy()
        if rho > 0: b = h / rho

    if len(h) == 1 and neval > 1: h = np.repeat(h, neval); b = np.repeat(b, neval)
    rho_vec = h / b

    ests = np.full((neval, 8), np.nan)
    for i in range(neval):
        h_i = h[i]; b_i = b[i]
        if bwcheck is not None:
            bw_min = np.sort(np.abs(x - eval[i]))[bwcheck - 1]
            h_i = max(h_i, bw_min)
            b_i = max(b_i, bw_min)
            rho_vec[i] = h_i / b_i

        u = (x - eval[i]) / h_i
        K_d = kd_K_fun(u, v=p, r=deriv, kernel=kernel)
        L_r = kd_K_fun(rho_vec[i] * u, v=p + 2, r=p, kernel=kernel)
        K = K_d[0]
        M = K - rho_vec[i] ** (1 + p) * L_r[0] * L_r[1]

        f_us = np.mean(K) / h_i
        f_bc = np.mean(M) / h_i
        se_us = np.sqrt((np.mean(K ** 2) - np.mean(K) ** 2) / (N * h_i ** 2))
        se_rb = np.sqrt((np.mean(M ** 2) - np.mean(M) ** 2) / (N * h_i ** 2))

        bw_eff = max(h_i, b_i)
        eN = int(np.sum(np.abs(x - eval[i]) <= bw_eff))

        ests[i, :] = [eval[i], h_i, b_i, eN, f_us, f_bc, se_us, se_rb]

    est_df = pd.DataFrame(ests, columns=["eval", "h", "b", "N",
                                          "tau.us", "tau.bc", "se.us", "se.rb"])
    opt = dict(p=p, kernel=kernel_type, n=N, neval=neval,
               bwselect=bwselect if bwselect else "Manual")
    return KdrobustResult(Estimate=est_df, opt=opt)
