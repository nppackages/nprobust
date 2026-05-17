"""Local polynomial regression with robust bias-corrected inference."""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from ._helpers import (
    build_dups,
    lprobust_cluster_meat,
    lprobust_res,
    lprobust_vce,
    qrXXinv,
    resolve_name,
)
from ._kernels import W_fun


_KERNEL_ALIASES = {
    "epa": "epa", "epanechnikov": "epa",
    "uni": "uni", "uniform": "uni",
    "tri": "tri", "triangular": "tri",
    "gau": "gau", "gaussian": "gau",
}
_KERNEL_TYPE_NAMES = {"epa": "Epanechnikov", "uni": "Uniform",
                      "tri": "Triangular", "gau": "Gaussian"}


@dataclass
class LprobustResult:
    """Container for :func:`lprobust` results.

    Attributes
    ----------
    Estimate : pandas.DataFrame
        Point estimates and standard errors over the evaluation grid, of
        shape ``(neval, 8)``. Columns: ``eval, h, b, N, tau.us, tau.bc,
        se.us, se.rb``.
    opt : dict
        Options actually used. Keys: ``p, q, deriv, kernel, n, neval,
        bwselect, vce``.
    cov_us : numpy.ndarray or None
        Covariance matrix of the conventional (un-bias-corrected) estimator
        across the evaluation grid, of shape ``(neval, neval)``. ``None``
        unless ``covgrid=True`` was passed to :func:`lprobust`.
    cov_rb : numpy.ndarray or None
        Covariance matrix of the robust-bias-corrected estimator across the
        evaluation grid, of shape ``(neval, neval)``. ``None`` unless
        ``covgrid=True`` was passed to :func:`lprobust`.
    call : dict
        Optional dict for storing the call metadata (default empty).
    """
    Estimate: pd.DataFrame
    opt: Dict[str, Any]
    cov_us: Optional[np.ndarray] = None
    cov_rb: Optional[np.ndarray] = None
    call: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (f"<LprobustResult n={self.opt['n']} p={self.opt['p']} "
                f"deriv={self.opt['deriv']} kernel={self.opt['kernel']} "
                f"bwselect={self.opt['bwselect']} neval={self.opt['neval']}>")

    def summary(self, alpha: float = 0.05) -> pd.DataFrame:
        """Return a tidy summary DataFrame; mirrors summary.lprobust in R."""
        z = norm.ppf(1 - alpha / 2)
        df = self.Estimate.copy()
        df["CI.lower"] = df["tau.bc"] - z * df["se.rb"]
        df["CI.upper"] = df["tau.bc"] + z * df["se.rb"]
        return df[["eval", "h", "b", "N", "tau.us", "se.us",
                   "tau.bc", "se.rb", "CI.lower", "CI.upper"]]

    def tidy(self, conf_level: float = 0.95) -> pd.DataFrame:
        """broom-style tidy output."""
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
            "p": self.opt["p"], "q": self.opt["q"],
            "deriv": self.opt["deriv"],
            "kernel": self.opt["kernel"],
            "bwselect": self.opt["bwselect"],
        }])


def _normalize_kernel(kernel: str) -> str:
    k = _KERNEL_ALIASES.get(kernel.lower())
    if k is None:
        raise ValueError(f"kernel {kernel!r} not recognised")
    return k


def lprobust(
    y,
    x,
    eval=None,
    neval=None,
    p=None,
    deriv=None,
    h=None,
    b=None,
    rho: float = 1.0,
    kernel: str = "epa",
    bwselect: Optional[str] = None,
    bwcheck: Optional[int] = 21,
    bwregul: float = 1.0,
    imsegrid: int = 30,
    vce: str = "nn",
    covgrid: bool = False,
    cluster=None,
    nnmatch: int = 3,
    level: float = 95,
    interior: bool = False,
    subset=None,
    weights=None,
    masspoints: str = "check",
    data=None,
) -> LprobustResult:
    """Local polynomial regression with robust bias-corrected inference.

    Implements local polynomial regression point estimators, with robust
    bias-corrected confidence intervals and inference procedures developed in
    Calonico, Cattaneo and Farrell (2018), with related optimality results in
    Calonico, Cattaneo and Farrell (2022). It also implements other estimation
    and inference procedures available in the literature. Companion commands:
    :func:`lpbwselect` for bandwidth selection and :func:`plot_lprobust` for
    plotting.

    Parameters
    ----------
    y, x : array-like or column-name string (with ``data=``)
        Outcome and covariate (length-n each).
    eval : array-like or None
        Evaluation points. If ``None``, an equally spaced grid of
        ``neval`` points spanning the support of ``x`` is used.
    neval : int or None
        Number of equally spaced evaluation points. Default 30.
    p, deriv : int or None
        Polynomial order ``p`` and derivative ``deriv``. Defaults
        ``p = deriv + 1`` and ``deriv = 0``.
    h, b : float, array, or None
        Point-estimation bandwidth ``h`` and bias-correction bandwidth
        ``b``. If ``None`` they are chosen by ``bwselect``.
    rho : float
        Ratio ``h/b`` used when ``b`` is not supplied. Default 1.
    kernel : {"epa", "tri", "uni", "gau"}
        Kernel function. Default ``"epa"``.
    bwselect : str or None
        Bandwidth selector. One of ``"mse-dpi"``, ``"mse-rot"``,
        ``"imse-dpi"``, ``"imse-rot"``, ``"ce-dpi"``, or ``"ce-rot"``.
        Use :func:`lpbwselect` with ``bwselect="all"`` to report all
        available bandwidth selectors. Default ``"imse-dpi"``.
    bwcheck : int or None
        Minimum unique observations on each side of an evaluation point
        for the pilot bandwidth. Default 21.
    bwregul : float
        Scaling factor for the bandwidth regularization term.
    imsegrid : int
        Grid size for IMSE-based selectors. Default 30.
    vce : str
        Variance estimator. Without ``cluster``: ``"nn"`` (default, nearest
        neighbor with ``nnmatch`` matches), ``"hc0"``, ``"hc1"``, ``"hc2"``,
        or ``"hc3"``. With ``cluster``: ``"cr1"`` (default; Stata-style CR1
        multiplier ``((n-1)/(n-k)) * (G/(G-1))``), ``"cr2"`` (Bell-McCaffrey
        block-adjusted residuals), or ``"cr3"`` (block jackknife-style
        adjustment with ``(G-1)/G`` multiplier). With ``cluster``, ``hc0`` or
        ``hc1`` remap to ``cr1``, ``hc2`` remaps to ``cr2``, ``hc3`` remaps
        to ``cr3``, and ``nn`` silently defaults to ``cr1``.
    cluster : array-like or column-name string (optional)
        Cluster identifier; activates cluster-robust variance.
    nnmatch : int
        Number of NN matches for ``vce="nn"``. Default 3.
    covgrid : bool
        If True, also return ``(neval, neval)`` cross-grid covariance
        matrices ``cov_us`` and ``cov_rb``. Default False.
    level : float
        Confidence level in percent. Default 95.
    interior : bool
        If True, force interior bias-correction series. Default False.
    subset : array-like or column-name string, optional
        Boolean rule selecting observations to use.
    weights : array-like or column-name string, optional
        Non-negative observation weights of the same length as ``x``. User
        weights multiply the kernel weights in estimation, bandwidth
        selection, and variance computation (weighted least squares
        interpretation).
    masspoints : {"check", "off"}
        How to handle evaluation points whose bandwidth window contains few
        unique ``x`` values. Default ``"check"`` warns when there are fewer
        than ``p + 5`` unique values inside the main bandwidth; ``"off"``
        disables the check.
    data : pandas.DataFrame (optional)
        When supplied, ``y``, ``x``, ``cluster``, ``weights``, ``subset``
        may be column-name strings.

    Returns
    -------
    result : LprobustResult
        An object with the following attributes:

        - ``Estimate`` (``pandas.DataFrame``, shape ``(neval, 8)``): point
          estimates and standard errors over the evaluation grid. Columns:
          ``eval, h, b, N, tau.us, tau.bc, se.us, se.rb``.
        - ``cov_us`` (``numpy.ndarray`` of shape ``(neval, neval)`` or
          ``None``): conventional-estimator covariance matrix across the
          evaluation grid; ``None`` unless ``covgrid=True``.
        - ``cov_rb`` (``numpy.ndarray`` of shape ``(neval, neval)`` or
          ``None``): robust-bias-corrected covariance matrix across the
          evaluation grid; ``None`` unless ``covgrid=True``.
        - ``opt`` (``dict``): options actually used. Keys: ``p, q, deriv,
          kernel, n, neval, bwselect, vce``.
    """
    if data is not None:
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
        if idx.dtype == bool:
            x = x[idx]; y = y[idx]
            if cluster is not None: cluster = cluster[idx]
            if weights is not None: weights = weights[idx]
        else:
            raise TypeError("subset must be a boolean array.")

    na_mask = np.isfinite(x) & np.isfinite(y)
    if cluster is not None:
        na_mask &= ~pd.isna(cluster)
    if weights is not None:
        na_mask &= np.isfinite(weights) & (weights >= 0)

    x = x[na_mask]; y = y[na_mask]
    if cluster is not None: cluster = cluster[na_mask]
    if weights is None:
        weights = np.ones_like(x)
    else:
        weights = weights[na_mask]

    if deriv is not None and p is None:
        p = deriv + 1
    if p is None: p = 1
    if deriv is None: deriv = 0
    q = p + 1

    N = len(x)
    if bwcheck is not None and bwcheck > N:
        warnings.warn(f"bwcheck ({bwcheck}) is larger than sample size ({N}); reducing to N.")
        bwcheck = N

    if eval is None:
        if neval is None:
            eval = np.linspace(x.min(), x.max(), 30)
        else:
            eval = np.linspace(x.min(), x.max(), int(neval))
    eval = np.asarray(eval, dtype=float).ravel()
    neval = len(eval)

    if h is None and bwselect is None:
        bwselect = "mse-dpi" if neval == 1 else "imse-dpi"

    kernel = _normalize_kernel(kernel)
    if bwselect is not None:
        bwselect = bwselect.lower()
    vce = vce.lower()

    if masspoints not in ("check", "off"):
        raise ValueError("masspoints must be 'check' or 'off'.")
    if kernel not in ("epa", "uni", "tri", "gau"):
        raise ValueError(f"kernel {kernel!r} not supported for lprobust.")
    if vce not in ("nn", "hc0", "hc1", "hc2", "hc3", "cr1", "cr2", "cr3"):
        raise ValueError(f"vce {vce!r} not supported.")
    if p < 0 or deriv < 0 or nnmatch <= 0:
        raise ValueError("p, deriv must be >= 0 and nnmatch must be > 0.")
    if deriv > p:
        raise ValueError("deriv must be <= p.")
    if level <= 0 or level >= 100:
        raise ValueError("level must be in (0, 100).")
    if rho < 0:
        raise ValueError("rho must be non-negative.")

    # Cluster vce normalization: with cluster, only cr1/cr2/cr3 are valid;
    # without cluster, cr1/cr2/cr3 fall back to hc1/hc2/hc3.
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

    # User-facing display label preserved for the result object; internal
    # cr1/cr2/cr3 -> hc1/hc2/hc3 remap is deferred until AFTER lpbwselect()
    # is called, so lpbwselect sees the user-facing cr* label and does not
    # re-fire the cluster-vce warning the user has already seen.
    vce_user = vce

    if vce == "nn":
        order = np.argsort(x, kind="mergesort")
        x = x[order]; y = y[order]; weights = weights[order]
        if cluster is not None: cluster = cluster[order]

    # Bandwidth selection if not provided
    if h is not None:
        bwselect = "Manual"
        if rho > 0 and b is None:
            b = h / rho
        if rho == 0 and b is None:
            raise ValueError("When h is provided and rho=0, b must also be provided.")
        h = np.asarray(h, dtype=float).ravel()
        b = np.asarray(b, dtype=float).ravel()
        if len(h) == 1 and neval > 1:
            h = np.repeat(h, neval); b = np.repeat(b, neval)
    else:
        from .lpbwselect import lpbwselect
        lpbws = lpbwselect(
            y, x, eval=eval, deriv=deriv, p=p, vce=vce, cluster=cluster,
            bwselect=bwselect, interior=interior, kernel=kernel,
            bwcheck=bwcheck, bwregul=bwregul, imsegrid=imsegrid,
            subset=None, weights=weights, masspoints="off",
        )
        h = lpbws.bws.iloc[:, 1].to_numpy()
        b = lpbws.bws.iloc[:, 2].to_numpy()
        if rho > 0: b = h / rho
        if len(h) == 1 and neval > 1:
            h = np.repeat(h, neval); b = np.repeat(b, neval)

    # cr1/cr2/cr3 -> hc1/hc2/hc3 for downstream residual paths.
    if vce == "cr1": vce = "hc1"
    elif vce == "cr2": vce = "hc2"
    elif vce == "cr3": vce = "hc3"

    # Near-neighbour duplicates (for vce=nn)
    dups = np.zeros(N, dtype=int); dupsid = np.zeros(N, dtype=int)
    if vce == "nn":
        dups, dupsid = build_dups(x)

    ests = np.full((neval, 8), np.nan)
    cov_us = np.full((neval, neval), np.nan) if covgrid else None
    cov_rb = np.full((neval, neval), np.nan) if covgrid else None

    # covgrid builds its own (i,j) quantities below (mirrors the R version).

    def _compute_at(i: int):
        h_i = h[i]; b_i = b[i]
        if bwcheck is not None:
            bw_min = np.sort(np.abs(x - eval[i]))[bwcheck - 1]
            h_i = max(h_i, bw_min)
            b_i = max(b_i, bw_min)

        w_h = W_fun((x - eval[i]) / h_i, kernel) / h_i * weights
        w_b = W_fun((x - eval[i]) / b_i, kernel) / b_i * weights
        ind_h = w_h > 0; ind_b = w_b > 0
        ind = ind_b if h_i <= b_i else ind_h

        if masspoints == "check":
            n_unique = len(np.unique(x[ind_h]))
            if n_unique < p + 5:
                warnings.warn(
                    f"Only {n_unique} unique x values within bandwidth at eval="
                    f"{eval[i]:.4f} (p+5={p+5}); local polynomial may be unreliable. "
                    f"Set masspoints='off' to silence.")

        eN = int(ind.sum())
        eY = y[ind]; eX = x[ind]
        W_h = w_h[ind]; W_b = w_b[ind]
        eC = cluster[ind] if cluster is not None else None

        edups = np.zeros(1, dtype=int); edupsid = np.zeros(1, dtype=int)
        if vce == "nn":
            edups = dups[ind]; edupsid = dupsid[ind]

        u = (eX - eval[i]) / h_i
        R_q = (eX - eval[i])[:, None] ** np.arange(q + 1)
        R_p = R_q[:, : p + 1]

        L = (R_p * W_h[:, None]).T @ (u ** (p + 1))        # (p+1,)
        invG_q = qrXXinv(R_q * np.sqrt(W_b)[:, None])
        invG_p = qrXXinv(R_p * np.sqrt(W_h)[:, None])
        e_p1 = np.zeros(q + 1); e_p1[p + 1] = 1.0
        Q_q = ((R_p * W_h[:, None]).T
               - h_i ** (p + 1) * np.outer(L, e_p1) @ invG_q @ (R_q * W_b[:, None]).T).T

        beta_p = invG_p @ (R_p * W_h[:, None]).T @ eY
        beta_q = invG_q @ (R_q * W_b[:, None]).T @ eY
        beta_bc = invG_p @ Q_q.T @ eY

        tau_cl = math.factorial(deriv) * beta_p[deriv]
        tau_bc = math.factorial(deriv) * beta_bc[deriv]

        predicts_p = np.zeros(eN); predicts_q = np.zeros(eN)
        hii = np.zeros(1)
        if vce in ("hc0", "hc1", "hc2", "hc3"):
            predicts_p = R_p @ beta_p
            predicts_q = R_q @ beta_q
            if vce in ("hc2", "hc3") and eC is None:
                hii = np.einsum("ij,jk,ik->i", R_p, invG_p, R_p * W_h[:, None])

        if eC is None:
            res_h = lprobust_res(eX, eY, predicts_p, hii, vce, nnmatch, edups, edupsid, p + 1)
            if vce == "nn":
                res_b = res_h
            else:
                res_b = lprobust_res(eX, eY, predicts_q, hii, vce, nnmatch, edups, edupsid, q + 1)
            meat_cl = lprobust_vce(R_p * W_h[:, None], res_h, None)
            meat_bc = lprobust_vce(Q_q, res_b, None)
        else:
            cr_type = {"hc0": "CR0", "hc1": "CR1", "hc2": "CR2", "hc3": "CR3",
                       "nn": "CR1"}.get(vce, "CR1")
            if vce == "nn":
                res_h_raw = lprobust_res(eX, eY, predicts_p, hii, "nn",
                                          nnmatch, edups, edupsid, p + 1)
                res_b_raw = res_h_raw
            else:
                res_h_raw = (eY - predicts_p).reshape(-1, 1)
                res_b_raw = (eY - predicts_q).reshape(-1, 1)
            sqrtW = np.sqrt(W_h)
            X_std_h = R_p * sqrtW[:, None]
            r_std_h = res_h_raw.ravel() * sqrtW
            # k_override = q+1 aligns the CR1 df correction with the q-regression
            # that produced res_b_raw. Without it, k=ncol(Q_q)=p+1 was being used.
            meat_cl = lprobust_cluster_meat(X_std_h, r_std_h, eC, invG_p, cr_type)
            meat_bc = lprobust_cluster_meat(Q_q, res_b_raw.ravel(), eC, invG_p, cr_type, k_override=q + 1)

        V_cl = invG_p @ meat_cl @ invG_p
        V_bc = invG_p @ meat_bc @ invG_p
        se_cl = float(math.factorial(deriv) * np.sqrt(V_cl[deriv, deriv]))
        se_rb = float(math.factorial(deriv) * np.sqrt(V_bc[deriv, deriv]))

        return dict(
            eN=eN, tau_cl=tau_cl, tau_bc=tau_bc, se_cl=se_cl, se_rb=se_rb,
            h_i=h_i, b_i=b_i,
        )

    for i in range(neval):
        info = _compute_at(i)
        ests[i, :] = [eval[i], info["h_i"], info["b_i"], info["eN"],
                      info["tau_cl"], info["tau_bc"], info["se_cl"], info["se_rb"]]

    if covgrid:
        for i in range(neval):
            for j in range(i, neval):
                h_i, b_i = h[i], b[i]
                h_j, b_j = h[j], b[j]

                w_h_i = W_fun((x - eval[i]) / h_i, kernel) / h_i * weights
                w_b_i = W_fun((x - eval[i]) / b_i, kernel) / b_i * weights
                w_h_j = W_fun((x - eval[j]) / h_j, kernel) / h_j * weights
                w_b_j = W_fun((x - eval[j]) / b_j, kernel) / b_j * weights

                ind_i = (w_h_i > 0) | (w_b_i > 0)
                ind_j = (w_h_j > 0) | (w_b_j > 0)
                ind   = ind_i | ind_j
                eN_ij = int(ind.sum())

                eY = y[ind]; eX = x[ind]
                W_h_i = w_h_i[ind]; W_b_i = w_b_i[ind]
                W_h_j = w_h_j[ind]; W_b_j = w_b_j[ind]
                eC_ij = cluster[ind] if cluster is not None else None

                edups_ij = np.zeros(1, dtype=int); edupsid_ij = np.zeros(1, dtype=int)
                if vce == "nn":
                    edups_ij = dups[ind]; edupsid_ij = dupsid[ind]

                u_i = (eX - eval[i]) / h_i
                R_q_i = (eX - eval[i])[:, None] ** np.arange(q + 1)
                R_p_i = R_q_i[:, : p + 1]
                u_j = (eX - eval[j]) / h_j
                R_q_j = (eX - eval[j])[:, None] ** np.arange(q + 1)
                R_p_j = R_q_j[:, : p + 1]

                e_p1_ij = np.zeros(q + 1); e_p1_ij[p + 1] = 1.0

                L_i = (R_p_i * W_h_i[:, None]).T @ (u_i ** (p + 1))
                invG_q_i = qrXXinv(R_q_i * np.sqrt(W_b_i)[:, None])
                invG_p_i = qrXXinv(R_p_i * np.sqrt(W_h_i)[:, None])
                Q_q_i = ((R_p_i * W_h_i[:, None]).T
                         - h_i ** (p + 1) * np.outer(L_i, e_p1_ij) @ invG_q_i @ (R_q_i * W_b_i[:, None]).T).T
                beta_p_i = invG_p_i @ (R_p_i * W_h_i[:, None]).T @ eY
                beta_q_i = invG_q_i @ (R_q_i * W_b_i[:, None]).T @ eY

                L_j = (R_p_j * W_h_j[:, None]).T @ (u_j ** (p + 1))
                invG_q_j = qrXXinv(R_q_j * np.sqrt(W_b_j)[:, None])
                invG_p_j = qrXXinv(R_p_j * np.sqrt(W_h_j)[:, None])
                Q_q_j = ((R_p_j * W_h_j[:, None]).T
                         - h_j ** (p + 1) * np.outer(L_j, e_p1_ij) @ invG_q_j @ (R_q_j * W_b_j[:, None]).T).T
                beta_p_j = invG_p_j @ (R_p_j * W_h_j[:, None]).T @ eY
                beta_q_j = invG_q_j @ (R_q_j * W_b_j[:, None]).T @ eY

                hii_i = np.zeros(1); hii_j = np.zeros(1)
                predicts_p_i = predicts_q_i = np.zeros(eN_ij)
                predicts_p_j = predicts_q_j = np.zeros(eN_ij)
                if vce in ("hc0", "hc1", "hc2", "hc3"):
                    predicts_p_i = R_p_i @ beta_p_i
                    predicts_q_i = R_q_i @ beta_q_i
                    predicts_p_j = R_p_j @ beta_p_j
                    predicts_q_j = R_q_j @ beta_q_j
                    if vce in ("hc2", "hc3") and eC_ij is None:
                        hii_i = np.einsum("ij,jk,ik->i", R_p_i, invG_p_i, R_p_i * W_h_i[:, None])
                        hii_j = np.einsum("ij,jk,ik->i", R_p_j, invG_p_j, R_p_j * W_h_j[:, None])

                res_h_i = lprobust_res(eX, eY, predicts_p_i, hii_i, vce, nnmatch, edups_ij, edupsid_ij, p + 1)
                res_b_i = res_h_i if vce == "nn" else lprobust_res(eX, eY, predicts_q_i, hii_i, vce, nnmatch, edups_ij, edupsid_ij, q + 1)
                res_h_j = lprobust_res(eX, eY, predicts_p_j, hii_j, vce, nnmatch, edups_ij, edupsid_ij, p + 1)
                res_b_j = res_h_j if vce == "nn" else lprobust_res(eX, eY, predicts_q_j, hii_j, vce, nnmatch, edups_ij, edupsid_ij, q + 1)

                # NOTE: V_us_i carries factorial(deriv) (NOT factorial(deriv)^2)
                # so that cov(tau_i, tau_j) = (V_us_i @ V_us_j.T)[d, d] picks up
                # factorial(deriv)^2 from the cross-product -- matching the
                # scaling of se_us = factorial(deriv) * sqrt(V_cl[d,d]) in the
                # main loop. The previous implementation had factorial(deriv)^2
                # here, which gave factorial(deriv)^4 after the cross-product
                # (a factor factorial(deriv)^2 too large for deriv >= 2).
                fact_d = math.factorial(deriv)
                V_us_i = fact_d * invG_p_i @ (res_h_i.ravel()[:, None] * R_p_i * W_h_i[:, None]).T
                V_us_j = fact_d * invG_p_j @ (res_h_j.ravel()[:, None] * R_p_j * W_h_j[:, None]).T
                V_rb_i = fact_d * invG_p_i @ (res_b_i.ravel()[:, None] * Q_q_i).T
                V_rb_j = fact_d * invG_p_j @ (res_b_j.ravel()[:, None] * Q_q_j).T

                cov_us[i, j] = (V_us_i @ V_us_j.T)[deriv, deriv]
                cov_rb[i, j] = (V_rb_i @ V_rb_j.T)[deriv, deriv]
                cov_us[j, i] = cov_us[i, j]
                cov_rb[j, i] = cov_rb[i, j]

    est_df = pd.DataFrame(ests, columns=["eval", "h", "b", "N",
                                         "tau.us", "tau.bc", "se.us", "se.rb"])

    opt = dict(p=p, q=q, deriv=deriv,
               kernel=_KERNEL_TYPE_NAMES[kernel], n=N, neval=neval,
               bwselect=bwselect if bwselect else "Manual",
               vce=vce_user)

    return LprobustResult(Estimate=est_df, opt=opt, cov_us=cov_us, cov_rb=cov_rb)
