"""Bandwidth selection for local polynomial regression."""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from ._bw import (
    lpbwselect_ce_dpi,
    lpbwselect_imse_dpi,
    lpbwselect_imse_rot,
    lpbwselect_mse_dpi,
    lpbwselect_mse_rot,
)

_KERNEL_ALIASES = {
    "epa": "epa", "epanechnikov": "epa",
    "uni": "uni", "uniform": "uni",
    "tri": "tri", "triangular": "tri",
    "gau": "gau", "gaussian": "gau",
}
_KERNEL_TYPE_NAMES = {"epa": "Epanechnikov", "uni": "Uniform",
                      "tri": "Triangular", "gau": "Gaussian"}


@dataclass
class LpbwselectResult:
    """Container for :func:`lpbwselect` results.

    Attributes
    ----------
    bws : pandas.DataFrame
        Bandwidths over the evaluation grid. When ``bwselect != "all"`` the
        frame has shape ``(neval, 3)`` with columns ``eval, h, b``. When
        ``bwselect == "all"`` it widens to shape ``(neval, 9)`` with columns
        ``eval, h.mse.dpi, b.mse.dpi, h.mse.rot, b.mse.rot, h.ce.dpi,
        b.ce.dpi, h.ce.rot, b.ce.rot``.
    bws_imse : pandas.DataFrame or None
        Only populated when ``bwselect == "all"``: a ``(2, 2)`` table indexed
        ``["h", "b"]`` with columns ``["IMSE-DPI", "IMSE-ROT"]`` carrying
        the IMSE bandwidths. ``None`` otherwise.
    opt : dict
        Options actually used. Keys: ``n, neval, p, q, deriv, kernel,
        bwselect``.
    call : dict
        Optional dict for storing the call metadata (default empty).
    """
    bws: pd.DataFrame
    bws_imse: Optional[pd.DataFrame]
    opt: Dict[str, Any]
    call: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"<LpbwselectResult n={self.opt['n']} p={self.opt['p']} "
                f"deriv={self.opt['deriv']} kernel={self.opt['kernel']} "
                f"bwselect={self.opt['bwselect']} neval={self.opt['neval']}>")

    def tidy(self) -> pd.DataFrame:
        return self.bws.copy()

    def glance(self) -> pd.DataFrame:
        return pd.DataFrame([{
            "n": self.opt["n"], "neval": self.opt["neval"],
            "p": self.opt["p"], "q": self.opt["q"],
            "deriv": self.opt["deriv"],
            "kernel": self.opt["kernel"],
            "bwselect": self.opt["bwselect"],
        }])


def _normalize_kernel(kernel):
    k = _KERNEL_ALIASES.get(kernel.lower())
    if k is None: raise ValueError(f"kernel {kernel!r} not recognised")
    return k


def lpbwselect(
    y, x,
    eval=None, neval=None, p=None, deriv=None,
    kernel: str = "epa",
    bwselect: str = "mse-dpi",
    bwcheck: Optional[int] = 21,
    bwregul: float = 1.0,
    imsegrid: int = 30,
    vce: str = "nn",
    cluster=None,
    nnmatch: int = 3,
    interior: bool = False,
    subset=None,
    weights=None,
    masspoints: str = "check",
    data=None,
) -> LpbwselectResult:
    """Bandwidth selection for local polynomial regression.

    Implements bandwidth selectors for local polynomial regression point
    estimators and inference procedures developed in Calonico, Cattaneo and
    Farrell (2018), with related optimality results in Calonico, Cattaneo and
    Farrell (2022). It also implements other bandwidth selectors available in
    the literature. Companion command: :func:`lprobust`.

    If ``data`` is a pandas DataFrame, ``y``, ``x``, ``cluster``, ``weights``,
    and ``subset`` may be passed as column-name strings.

    Parameters
    ----------
    y, x : array-like or column-name string
        Dependent and independent variables.
    eval : array-like or None
        Evaluation points. If ``None``, an equally spaced grid of ``neval``
        points over the support of ``x`` is used.
    neval : int or None
        Number of equally spaced evaluation points. Default 30.
    p, deriv : int or None
        Polynomial order and derivative order. Defaults are ``p = 1`` and
        ``deriv = 0``.
    kernel : {"epa", "tri", "uni", "gau"}
        Kernel function. Default ``"epa"``.
    bwselect : str
        Bandwidth selector. Options are ``"mse-dpi"``, ``"mse-rot"``,
        ``"imse-dpi"``, ``"imse-rot"``, ``"ce-dpi"``, ``"ce-rot"``, and
        ``"all"``. Default ``"mse-dpi"``.
    bwcheck : int or None
        If positive, selected bandwidths are enlarged so at least
        ``bwcheck`` observations are available at each evaluation point.
        Default 21.
    bwregul : float
        Scaling factor for the bandwidth regularization term. Setting this to
        0 removes the regularization term. Default 1.
    imsegrid : int
        Number of evaluation points used to compute IMSE bandwidth selectors.
        Default 30.
    vce : str
        Variance estimator. Without ``cluster``: ``"nn"`` (default, nearest
        neighbor with ``nnmatch`` matches), ``"hc0"``, ``"hc1"``, ``"hc2"``,
        or ``"hc3"``. With ``cluster``: ``"cr1"`` (default; Stata-style CR1
        multiplier ``((n-1)/(n-k)) * (G/(G-1))``), ``"cr2"`` (Bell-McCaffrey
        block-adjusted residuals), or ``"cr3"`` (block jackknife-style
        adjustment with ``(G-1)/G`` multiplier). With ``cluster``, ``hc0`` or
        ``hc1`` remap to ``cr1``, ``hc2`` remaps to ``cr2``, ``hc3`` remaps
        to ``cr3``, and ``nn`` silently defaults to ``cr1``.
    cluster : array-like or column-name string, optional
        Cluster identifier for cluster-robust variance estimation.
    nnmatch : int
        Minimum number of neighbors for ``vce="nn"``. Default 3.
    interior : bool
        If True, all evaluation points are assumed to be interior points.
        Default False.
    subset : array-like or column-name string, optional
        Boolean rule selecting observations to use.
    weights : array-like or column-name string, optional
        Non-negative observation weights, multiplicative with kernel weights
        in all bandwidth-selection steps.
    masspoints : {"check", "off"}
        How to handle evaluation points with few unique ``x`` values within a
        bandwidth. Default ``"check"``; ``"off"`` disables the check.
    data : pandas.DataFrame, optional
        Data source for column-name inputs.

    Returns
    -------
    result : LpbwselectResult
        An object with the following attributes:

        - ``bws`` (``pandas.DataFrame``): bandwidths over the evaluation
          grid. Shape ``(neval, 3)`` with columns ``eval, h, b`` when
          ``bwselect != "all"``; widens to shape ``(neval, 9)`` with
          columns ``eval, h.mse.dpi, b.mse.dpi, h.mse.rot, b.mse.rot,
          h.ce.dpi, b.ce.dpi, h.ce.rot, b.ce.rot`` when
          ``bwselect == "all"``.
        - ``bws_imse`` (``pandas.DataFrame`` or ``None``): only populated
          when ``bwselect == "all"`` - a ``(2, 2)`` table indexed
          ``["h", "b"]`` with columns ``["IMSE-DPI", "IMSE-ROT"]``
          carrying the IMSE bandwidths.
        - ``opt`` (``dict``): options actually used. Keys: ``n, neval, p,
          q, deriv, kernel, bwselect``.
    """
    if data is not None:
        from ._helpers import resolve_name
        y       = resolve_name(y,       data, "y")
        x       = resolve_name(x,       data, "x")
        cluster = resolve_name(cluster, data, "cluster")
        weights = resolve_name(weights, data, "weights")
        subset  = resolve_name(subset,  data, "subset")

    y = np.asarray(y, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    if cluster is not None: cluster = np.asarray(cluster).ravel()
    if weights is not None: weights = np.asarray(weights, dtype=float).ravel()

    if subset is not None:
        idx = np.asarray(subset).ravel()
        if idx.dtype != bool:
            raise TypeError("subset must be a boolean array.")
        x = x[idx]; y = y[idx]
        if cluster is not None: cluster = cluster[idx]
        if weights is not None: weights = weights[idx]

    na_mask = np.isfinite(x) & np.isfinite(y)
    if cluster is not None: na_mask &= ~pd.isna(cluster)
    if weights is not None: na_mask &= np.isfinite(weights) & (weights >= 0)
    x = x[na_mask]; y = y[na_mask]
    if cluster is not None: cluster = cluster[na_mask]
    if weights is None: weights = np.ones_like(x)
    else: weights = weights[na_mask]

    if masspoints not in ("check", "off"):
        raise ValueError("masspoints must be 'check' or 'off'.")

    N = len(x)
    if bwcheck is not None and bwcheck > N:
        warnings.warn(f"bwcheck ({bwcheck}) > N ({N}); reducing bwcheck to N.")
        bwcheck = N

    if deriv is not None and p is None: p = deriv + 1
    if p is None: p = 1
    if deriv is None: deriv = 0
    q = p + 1

    if eval is None:
        if neval is None:
            eval = np.linspace(x.min(), x.max(), 30)
        else:
            eval = np.linspace(x.min(), x.max(), int(neval))
    eval = np.asarray(eval, dtype=float).ravel()
    neval = len(eval)

    if bwselect in ("imse-dpi", "imse-rot"):
        neval = 1
        eval = np.array([1.0])  # unused; IMSE is grid-averaged internally

    even = (p - deriv) % 2 == 0

    kernel = _normalize_kernel(kernel)
    bwselect = bwselect.lower()
    vce = vce.lower()

    if vce not in ("nn", "hc0", "hc1", "hc2", "hc3", "cr1", "cr2", "cr3"):
        raise ValueError(f"vce {vce!r} not supported.")

    # Cluster vce normalization mirrors lprobust(): with cluster, only
    # cr1/cr2/cr3 are valid; without cluster, cr1/cr2/cr3 fall back to
    # hc1/hc2/hc3. Internal helpers only know about nn/hc0..hc3, so cr*
    # gets remapped to hc*.
    if cluster is not None:
        if vce == "nn":
            vce = "cr1"  # silent default
        elif vce in ("hc0", "hc1"):
            warnings.warn(f"vce={vce!r} is not a cluster option. Switching to vce='cr1'.")
            vce = "cr1"
        elif vce == "hc2":
            warnings.warn("vce='hc2' is not a cluster option. Switching to vce='cr2'.")
            vce = "cr2"
        elif vce == "hc3":
            warnings.warn("vce='hc3' is not a cluster option. Switching to vce='cr3'.")
            vce = "cr3"
    else:
        if vce == "cr1":
            warnings.warn("vce='cr1' requires a cluster variable. Falling back to vce='hc1'.")
            vce = "hc1"
        elif vce == "cr2":
            warnings.warn("vce='cr2' requires a cluster variable. Falling back to vce='hc2'.")
            vce = "hc2"
        elif vce == "cr3":
            warnings.warn("vce='cr3' requires a cluster variable. Falling back to vce='hc3'.")
            vce = "hc3"
    if vce == "cr1": vce = "hc1"
    elif vce == "cr2": vce = "hc2"
    elif vce == "cr3": vce = "hc3"

    if bwselect == "all":
        bws = np.full((neval, 8), np.nan)
        colnames = ["h.mse.dpi", "b.mse.dpi", "h.mse.rot", "b.mse.rot",
                    "h.ce.dpi", "b.ce.dpi", "h.ce.rot", "b.ce.rot"]
        bws_imse = np.full((2, 2), np.nan)
    else:
        bws = np.full((neval, 2), np.nan)
        colnames = ["h", "b"]
        bws_imse = None

    h_imse_dpi = b_imse_dpi = h_imse_rot = b_imse_rot = None
    if bwselect in ("imse-dpi", "all"):
        est = lpbwselect_imse_dpi(y, x, cluster, p, q, deriv, kernel,
                                bwcheck, bwregul, imsegrid, vce, nnmatch,
                                interior, weights)
        h_imse_dpi = est["h"]; b_imse_dpi = est["b"]
        if bwselect == "imse-dpi":
            bws[0, :2] = [h_imse_dpi, b_imse_dpi]

    if bwselect in ("imse-rot", "all"):
        h_imse_rot = lpbwselect_imse_rot(y, x, p, deriv, kernel, imsegrid)["h"]
        b_imse_rot = lpbwselect_imse_rot(y, x, q, p + 1, kernel, imsegrid)["h"]
        if bwselect == "imse-rot":
            bws[0, :2] = [h_imse_rot, b_imse_rot]

    if bwselect == "all":
        bws_imse[:, 0] = [h_imse_dpi, b_imse_dpi]
        bws_imse[:, 1] = [h_imse_rot, b_imse_rot]

    if bwselect in ("all", "mse-dpi", "mse-rot", "ce-dpi", "ce-rot"):
        # per-eval loop; results written directly into the full `bws` matrix.
        for i in range(neval):
            h_mse_dpi = b_mse_dpi = h_mse_rot = b_mse_rot = np.nan
            h_ce_dpi = b_ce_dpi = h_ce_rot = b_ce_rot = np.nan

            if bwselect in ("mse-dpi", "ce-dpi", "ce-rot", "all"):
                est = lpbwselect_mse_dpi(y, x, cluster, eval[i], p, q, deriv,
                                       kernel, bwcheck, bwregul, vce, nnmatch,
                                       interior, weights)
                h_mse_dpi, b_mse_dpi = est["h"], est["b"]
                if bwselect == "mse-dpi":
                    bws[i, :2] = [h_mse_dpi, b_mse_dpi]

            if bwselect in ("mse-rot", "all"):
                h_mse_rot = lpbwselect_mse_rot(y, x, eval[i], p, deriv, kernel)["h"]
                b_mse_rot = lpbwselect_mse_rot(y, x, eval[i], q, p + 1, kernel)["h"]
                if bwselect == "mse-rot":
                    bws[i, :2] = [h_mse_rot, b_mse_rot]

            if bwselect in ("ce-dpi", "all"):
                if even:
                    h_ce_dpi = h_mse_dpi * N ** (-((p + 2) / ((2 * p + 5) * (p + 3))))
                    b_ce_dpi = b_mse_dpi * N ** (-((q)     / ((2 * q + 3) * (q + 3))))
                else:
                    est = lpbwselect_ce_dpi(y, x, h_mse_dpi, b_mse_dpi, eval[i], p, q,
                                          deriv, rho=1, kernel=kernel, vce=vce,
                                          nnmatch=nnmatch, interior=interior,
                                          bwregul=bwregul, weights=weights)
                    h_ce_dpi = est["h"]
                    b_ce_dpi = b_mse_dpi * N ** (-((q + 2) / ((2 * q + 5) * (q + 3))))
                if bwselect == "ce-dpi":
                    bws[i, :2] = [h_ce_dpi, b_ce_dpi]

            if bwselect in ("ce-rot", "all"):
                if even:
                    h_ce_rot = h_mse_dpi * N ** (-((p + 2) / ((2 * p + 5) * (p + 3))))
                    b_ce_rot = b_mse_dpi * N ** (-((q)     / ((2 * q + 3) * (q + 3))))
                else:
                    h_ce_rot = h_mse_dpi * N ** (-((p)     / ((2 * p + 3) * (p + 3))))
                    b_ce_rot = b_mse_dpi * N ** (-((q + 2) / ((2 * q + 5) * (q + 3))))
                if bwselect == "ce-rot":
                    bws[i, :2] = [h_ce_rot, b_ce_rot]

            if bwselect == "all":
                bws[i, :] = [h_mse_dpi, b_mse_dpi, h_mse_rot, b_mse_rot,
                             h_ce_dpi, b_ce_dpi, h_ce_rot, b_ce_rot]

    # Build DataFrame
    bws_df = pd.DataFrame(np.column_stack([eval, bws]),
                          columns=["eval", *colnames])
    bws_imse_df = None
    if bws_imse is not None:
        bws_imse_df = pd.DataFrame(bws_imse, columns=["IMSE-DPI", "IMSE-ROT"],
                                   index=["h", "b"])

    opt = dict(n=N, neval=neval, p=p, q=q, deriv=deriv,
               kernel=_KERNEL_TYPE_NAMES[kernel], bwselect=bwselect)

    return LpbwselectResult(bws=bws_df, bws_imse=bws_imse_df, opt=opt)
