"""Residual / variance helpers for nprobust.

Port of :func:`lprobust.res`, :func:`lprobust.vce`, and
:func:`lprobust.cluster.meat` from the R package.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import linalg


def resolve_name(arg, data, argname):
    """Resolve a column-name string `arg` against a pandas DataFrame.

    If `arg` is a single string and `data` is a DataFrame, return
    ``data[arg]`` (a Series). Otherwise return `arg` unchanged.
    Used by lprobust/lpbwselect/kdrobust/kdbwselect so callers can pass
    ``y="yname"`` + ``data=df`` as an alternative to ``y=df["yname"]``.
    """
    if arg is None:
        return None
    if isinstance(arg, str):
        if data is None:
            raise ValueError(
                f"String `{argname}` requires `data = ` (a pandas DataFrame)."
            )
        if arg not in data.columns:
            raise ValueError(
                f"Column {arg!r} not found in `data` (for `{argname}`)."
            )
        return data[arg]
    return arg


# ---------------------------------------------------------------------------
def qrXXinv(X: np.ndarray) -> np.ndarray:
    """Compute (X'X)^{-1} via Cholesky, matching R's qrXXinv."""
    XtX = X.T @ X
    try:
        L = linalg.cholesky(XtX, lower=True)
    except linalg.LinAlgError:
        return linalg.pinv(XtX)
    invL = linalg.solve_triangular(L, np.eye(L.shape[0]), lower=True)
    return invL.T @ invL


# ---------------------------------------------------------------------------
def lprobust_res(
    X: np.ndarray,
    y: np.ndarray,
    m: np.ndarray,
    hii: np.ndarray,
    vce: str,
    matches: int,
    dups: np.ndarray,
    dupsid: np.ndarray,
    d: int,
) -> np.ndarray:
    """Compute residuals for the variance step.

    Port of :func:`lprobust.res` from R. Returns an (N, 1) column vector.
    """
    n = len(y)
    res = np.empty((n, 1), dtype=float)

    if vce == "nn":
        dups = np.asarray(dups, dtype=int)
        dupsid = np.asarray(dupsid, dtype=int)
        for pos in range(n):
            rpos = dups[pos] - dupsid[pos]
            lpos = dupsid[pos] - 1
            while (lpos + rpos) < min(matches, n - 1):
                if (pos - lpos - 1) < 0:
                    rpos += dups[pos + rpos + 1]
                elif (pos + rpos + 1) >= n:
                    lpos += dups[pos - lpos - 1]
                else:
                    left_gap  = X[pos] - X[pos - lpos - 1]
                    right_gap = X[pos + rpos + 1] - X[pos]
                    if left_gap > right_gap:
                        rpos += dups[pos + rpos + 1]
                    elif left_gap < right_gap:
                        lpos += dups[pos - lpos - 1]
                    else:
                        rpos += dups[pos + rpos + 1]
                        lpos += dups[pos - lpos - 1]
            start = pos - lpos
            stop  = min(n - 1, pos + rpos)
            idx = np.arange(start, stop + 1)
            yJ = np.sum(y[idx]) - y[pos]
            Ji = len(idx) - 1
            res[pos, 0] = np.sqrt(Ji / (Ji + 1)) * (y[pos] - yJ / Ji)
        return res

    # heteroskedasticity-consistent variants
    m = np.asarray(m).reshape(-1)
    if vce == "hc0":
        w = np.ones(n)
    elif vce == "hc1":
        w = np.full(n, np.sqrt(n / max(n - d, 1)))
    elif vce == "hc2":
        w = np.sqrt(1.0 / (1.0 - np.asarray(hii).reshape(-1)))
    elif vce == "hc3":
        w = 1.0 / (1.0 - np.asarray(hii).reshape(-1))
    else:
        raise ValueError(f"unknown vce {vce!r}")
    res[:, 0] = w * (y - m)
    return res


# ---------------------------------------------------------------------------
def lprobust_vce(RX: np.ndarray, res: np.ndarray, C: Optional[np.ndarray]) -> np.ndarray:
    """Return the meat matrix M for the HC / NN sandwich.

    For non-null ``C`` applies a CR1-style dof multiplier (legacy; the
    main entry points use :func:`lprobust_cluster_meat` with an explicit
    ``cr_type`` for CR0/CR1/CR2/CR3).
    """
    k = RX.shape[1]
    M = np.zeros((k, k), dtype=float)
    res = np.asarray(res).reshape(-1, 1)

    if C is None:
        M = (RX * res).T @ (RX * res)
        return M

    C = np.asarray(C)
    clusters = np.unique(C)
    g = len(clusters)
    n = len(C)
    w = ((n - 1) / max(n - k, 1)) * (g / max(g - 1, 1))
    for cl in clusters:
        idx = (C == cl)
        Xi = RX[idx, :]
        ri = res[idx, :]
        score = Xi.T @ ri
        M = M + score @ score.T
    return w * M


# ---------------------------------------------------------------------------
def lprobust_cluster_meat(
    X: np.ndarray,
    r: np.ndarray,
    C: np.ndarray,
    invG: np.ndarray,
    cr_type: str,
    k_override: Optional[int] = None,
) -> np.ndarray:
    """Cluster-robust meat for CR0/CR1/CR2/CR3.

    Mirror of :func:`lprobust.cluster.meat` in the R package.

    ``k_override``: when supplied, overrides the (N-1)/(N-k) df correction's
    k for CR1. Used when ``X`` is a "score-like" matrix (e.g. ``Q_q``) whose
    ncol is smaller than the effective number of regressors that produced
    ``r``. Without it, the bias-corrected CR1 path used k = p+1 from ncol(Q_q),
    inconsistent with the q-regression that produced the residuals (k = q+1).
    """
    C = np.asarray(C)
    clusters = np.unique(C)
    G = len(clusters)
    N = len(r)
    k    = X.shape[1]
    k_df = k if k_override is None else k_override
    r = np.asarray(r).reshape(-1)

    M = np.zeros((k, k), dtype=float)
    for cl in clusters:
        idx = (C == cl)
        Xg = X[idx, :]
        rg = r[idx]

        if cr_type in ("CR0", "CR1"):
            r_adj = rg
        else:
            Hgg = Xg @ invG @ Xg.T
            Hgg = 0.5 * (Hgg + Hgg.T)
            I_H = np.eye(Xg.shape[0]) - Hgg
            if cr_type == "CR2":
                vals, vecs = linalg.eigh(I_H)
                vals = np.maximum(vals, 1e-12)
                r_adj = vecs @ ((vecs.T @ rg) / np.sqrt(vals))
            elif cr_type == "CR3":
                try:
                    r_adj = linalg.solve(I_H, rg, assume_a="sym")
                except linalg.LinAlgError:
                    vals, vecs = linalg.eigh(I_H)
                    vals = np.maximum(vals, 1e-12)
                    r_adj = vecs @ ((vecs.T @ rg) / vals)
            else:
                raise ValueError(f"unknown cr_type {cr_type!r}")

        score = Xg.T @ r_adj
        M += np.outer(score, score)

    if cr_type == "CR0":
        mult = 1.0
    elif cr_type == "CR1":
        mult = ((N - 1) / max(N - k_df, 1)) * (G / max(G - 1, 1))
    elif cr_type == "CR2":
        mult = 1.0
    elif cr_type == "CR3":
        mult = (G - 1) / max(G, 1)
    else:
        raise ValueError(f"unknown cr_type {cr_type!r}")
    return mult * M


# ---------------------------------------------------------------------------
def build_dups(x_sorted: np.ndarray):
    """Return (dups, dupsid) arrays from sorted input in O(N).

    Equivalent to R's rle-based implementation: ``dups[j]`` is the size of
    the run containing ``x_sorted[j]``; ``dupsid[j]`` counts 1..dups[j]
    within each run.
    """
    x = np.asarray(x_sorted)
    n = len(x)
    if n == 0:
        return np.empty(0, dtype=int), np.empty(0, dtype=int)
    # Boundaries between runs on sorted input
    boundaries = np.flatnonzero(x[1:] != x[:-1]) + 1
    starts  = np.concatenate([[0], boundaries])
    lengths = np.diff(np.concatenate([starts, [n]]))
    dups   = np.repeat(lengths, lengths)
    dupsid = np.concatenate([np.arange(1, L + 1) for L in lengths])
    return dups, dupsid
