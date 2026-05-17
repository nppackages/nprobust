"""Bandwidth-selection internals.

Ports :func:`lp.bw.fun`, :func:`lprobust.bw`, and the MSE/IMSE/CE
helpers from ``npfunctions.R``.
"""
from __future__ import annotations

import math
from typing import Dict

import numpy as np
from scipy import integrate, linalg, optimize

from ._helpers import (
    build_dups,
    lprobust_res,
    lprobust_vce,
    qrXXinv,
)
from ._kernels import W_fun
from ._lpbwce import lpbwce


_EPS = np.finfo(float).eps


# ---------------------------------------------------------------------------
def lp_bw_fun(V: float, Bsq: float, p: int, v: int, N: int, kernel: str) -> Dict[str, float]:
    """Rule-of-thumb / plug-in constant."""
    def k_fun(u):
        return W_fun(u, kernel)

    def m1(i, j, k):
        return integrate.quad(lambda x: x ** i * x ** j * k(x), 0.0, np.inf)[0]

    def m2(i, j, k):
        return integrate.quad(lambda x: x ** i * x ** j * k(x) ** 2, 0.0, np.inf)[0]

    def GAMMA(p0, k):
        G = np.empty((p0 + 1, p0 + 1))
        for i in range(p0 + 1):
            for j in range(p0 + 1):
                G[i, j] = m1(i, j, k)
        return G

    def NU(p0, k):
        out = np.empty((p0 + 1, 1))
        for i in range(p0 + 1):
            out[i, 0] = m1(i, p0 + 1, k)
        return out

    def PSI(p0, k):
        P = np.empty((p0 + 1, p0 + 1))
        for i in range(p0 + 1):
            for j in range(p0 + 1):
                P[i, j] = m2(i, j, k)
        return P

    def C1_fun(p0, vv):
        Sinv = linalg.inv(GAMMA(p0, k_fun))
        return float((Sinv @ NU(p0, k_fun))[vv, 0])

    def C2_fun(p0, vv):
        Sinv = linalg.inv(GAMMA(p0, k_fun))
        return float((Sinv @ PSI(p0, k_fun) @ Sinv)[vv, vv])

    C1 = C1_fun(p, v)
    C2 = C2_fun(p, v)
    bw = (((2 * v + 1) * C2 * V) / (2 * (p + 1 - v) * C1 ** 2 * Bsq * N)) ** (1.0 / (2 * p + 3))
    return {"bw": bw, "C1": C1, "C2": C2}


# ---------------------------------------------------------------------------
def lprobust_bw(
    Y, X, cluster, c, o, nu, oB, hV, hB1, hB2, scale, vce, nnmatch, kernel,
    dups, dupsid, weights=None,
):
    """Port of :func:`lprobust.bw` (the pilot-variance builder)."""
    Y = np.asarray(Y, dtype=float).ravel()
    X = np.asarray(X, dtype=float).ravel()
    Nfull = len(X)
    if weights is None:
        weights = np.ones(Nfull)
    else:
        weights = np.asarray(weights, dtype=float).ravel()

    ## Variance block (at h.V)
    w_arr = W_fun((X - c) / hV, kernel) / hV * weights
    ind_V = w_arr > 0
    eY = Y[ind_V]; eX = X[ind_V]; eW = w_arr[ind_V]
    R_V = eX[:, None] - c
    R_V = R_V ** np.arange(o + 1)
    invG_V = qrXXinv(R_V * np.sqrt(eW)[:, None])
    beta_V = invG_V @ (R_V * eW[:, None]).T @ eY

    dups_V = dupsid_V = np.zeros(1, dtype=int)
    predicts_V = np.zeros(1)
    hii = np.zeros(1)

    eC = None
    if cluster is not None:
        eC = cluster[ind_V]

    if vce == "nn":
        dups_V   = dups[ind_V]
        dupsid_V = dupsid[ind_V]

    if vce in ("hc0", "hc1", "hc2", "hc3"):
        predicts_V = R_V @ beta_V
        if vce in ("hc2", "hc3"):
            hii = np.einsum("ij,jk,ik->i", R_V, invG_V, R_V * eW[:, None])

    res_V = lprobust_res(eX, eY, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o + 1)
    meat_V = lprobust_vce(R_V * eW[:, None], res_V, eC)
    V_V    = float((invG_V @ meat_V @ invG_V)[nu, nu])

    ## Bias block
    Hp = hV ** np.arange(o + 1)
    v1 = (R_V * eW[:, None]).T @ ((eX - c) / hV) ** (o + 1)
    v2 = (R_V * eW[:, None]).T @ ((eX - c) / hV) ** (o + 2)
    BConst1 = float((Hp * (invG_V @ v1))[nu])
    BConst2 = float((Hp * (invG_V @ v2))[nu])

    w_arr = W_fun((X - c) / hB1, kernel) * weights
    ind = w_arr > 0
    eY = Y[ind]; eX = X[ind]; eW = w_arr[ind]

    eC2 = None
    if cluster is not None:
        eC2 = cluster[ind]

    R_B1 = eX[:, None] - c
    R_B1 = R_B1 ** np.arange(oB + 1)
    invG_B1 = qrXXinv(R_B1 * np.sqrt(eW)[:, None])
    beta_B1 = invG_B1 @ (R_B1 * eW[:, None]).T @ eY

    BWreg = 0.0
    if scale > 0:
        dups_B = dupsid_B = np.zeros(1, dtype=int)
        hii_B = np.zeros(1)
        predicts_B = np.zeros(1)
        if vce == "nn":
            dups_B   = dups[ind]
            dupsid_B = dupsid[ind]
        if vce in ("hc0", "hc1", "hc2", "hc3"):
            predicts_B = R_B1 @ beta_B1
            if vce in ("hc2", "hc3"):
                hii_B = np.einsum("ij,jk,ik->i", R_B1, invG_B1, R_B1 * eW[:, None])
        res_B = lprobust_res(eX, eY, predicts_B, hii_B, vce, nnmatch, dups_B, dupsid_B, oB + 1)
        meat_B = lprobust_vce(R_B1 * eW[:, None], res_B, eC2)
        V_B = float((invG_B1 @ meat_B @ invG_B1)[o + 1, o + 1])
        BWreg = 3 * BConst1 ** 2 * V_B

    w_arr = W_fun((X - c) / hB2, kernel) * weights
    ind = w_arr > 0
    eY = Y[ind]; eX = X[ind]; eW = w_arr[ind]
    R_B2 = eX[:, None] - c
    R_B2 = R_B2 ** np.arange(oB + 2)
    invG_B2 = qrXXinv(R_B2 * np.sqrt(eW)[:, None])
    beta_B2 = invG_B2 @ (R_B2 * eW[:, None]).T @ eY

    N = Nfull
    B1 = BConst1 * beta_B1[o + 1]
    B2 = BConst2 * beta_B2[o + 2]
    V  = N * hV ** (2 * nu + 1) * V_V
    R  = BWreg
    r  = 1.0 / (2 * o + 3)
    rB = 2 * (o + 1 - nu)
    rV = 2 * nu + 1
    bw = ((rV * V) / (N * rB * (B1 ** 2 + scale * R))) ** r
    return {"V": V, "B1": float(B1), "B2": float(B2), "R": R, "r": r,
            "rB": rB, "rV": rV, "bw": float(bw)}


# ---------------------------------------------------------------------------
def lpbwselect_mse_rot(y, x, eval_pt, p, deriv, kernel):
    even = (p - deriv) % 2 == 0
    if kernel == "uni":
        C_c = 1.843
    elif kernel == "tri":
        C_c = 2.576
    else:
        C_c = 2.34
    x_iq = np.quantile(x, 0.75) - np.quantile(x, 0.25)
    x_max, x_min = np.max(x), np.min(x)
    rng = x_max - x_min
    N = len(x)

    c_bw = C_c * min(np.std(x, ddof=1), x_iq / 1.349) * N ** (-1.0 / 5)
    n_h1 = int(np.sum(np.abs(x - eval_pt) <= c_bw))
    f_hat = n_h1 / max(2 * N * c_bw, 1e-12)

    # Pilot sigma^2 from high-order polynomial fit (R uses k+1 = p+4 columns,
    # x^0..x^(p+3); see npfunctions.R:851-855). Scale x to avoid
    # ill-conditioning on the Vandermonde — R lm() uses pivoted QR, which
    # handles this; scipy.lstsq uses SVD and returns minimum-norm.
    k = p + 3
    s_x = float(np.max(np.abs(x))) if np.max(np.abs(x)) > 0 else 1.0
    u = x / s_x
    r_k = u[:, None] ** np.arange(k + 1)
    gamma, *_ = linalg.lstsq(r_k, y)
    fit = r_k @ gamma
    resid = y - fit
    s2_hat = float(np.sum(resid ** 2) / max(len(y) - (k + 1), 1))

    # m_{p+1}(eval_pt), m_{p+2}(eval_pt) via a local-poly call with the full range.
    # Delayed import to avoid a circular dep on lprobust.
    from .lprobust import lprobust  # noqa: WPS433

    mp1 = lprobust(y, x, eval=np.array([eval_pt]), h=rng,
                   p=p + 3, deriv=p + 1, kernel=kernel, vce="nn",
                   masspoints="off").Estimate["tau.us"].iloc[0]
    mp2 = lprobust(y, x, eval=np.array([eval_pt]), h=rng,
                   p=p + 3, deriv=p + 2, kernel=kernel, vce="nn",
                   masspoints="off").Estimate["tau.us"].iloc[0]

    bw = lp_bw_fun(s2_hat / f_hat, (mp1 / math.factorial(p + 1)) ** 2, p, deriv, N, kernel)
    B1 = bw["C1"] * mp1 / math.factorial(p + 1)
    B2 = bw["C1"] * mp2 / math.factorial(p + 2)
    V  = bw["C2"] * s2_hat / f_hat

    if not even:
        h_rot = bw["bw"]
        B = B1
    else:
        def obj(H):
            return abs(H ** (2 * p + 2 - 2 * deriv) * (B1 + H * B2) ** 2
                       + V / (N * H ** (1 + 2 * deriv)))
        # R optimize() uses .Machine$double.eps^0.25 by default. Match that
        # tolerance here because this CE-DPI objective can be very flat near
        # the minimum, so over-optimizing produces tiny cross-language drift.
        r_opt_tol = np.finfo(float).eps ** 0.25
        res = optimize.minimize_scalar(obj, bounds=(_EPS, rng), method="bounded", options={"xatol": r_opt_tol})
        h_rot = float(res.x)
        B = B1 + h_rot * B2
    return {"h": float(h_rot), "V": float(V), "B": float(B)}


# ---------------------------------------------------------------------------
def lpbwselect_imse_rot(y, x, p, deriv, kernel, imsegrid):
    eval_grid = np.quantile(x, np.arange(0.05, 0.9501, 0.025))
    N = len(x)
    x_max, x_min = np.max(x), np.min(x)
    rng = x_max - x_min
    even = (p - deriv) % 2 == 0
    V = np.empty(len(eval_grid))
    B = np.empty(len(eval_grid))
    for i, e in enumerate(eval_grid):
        r = lpbwselect_mse_rot(y, x, e, p, deriv, kernel)
        V[i] = r["V"]; B[i] = r["B"]
    if not even:
        h_rot = (np.mean((1 + 2 * deriv) * V)
                 / (N * np.mean(2 * (p + 1 - deriv) * B ** 2))) ** (1.0 / (2 * p + 3))
    else:
        def obj(H):
            return abs(H ** (2 * p + 2 - 2 * deriv) * np.mean(B ** 2)
                       + np.mean(V) / (N * H ** (1 + 2 * deriv)))
        r_opt_tol = np.finfo(float).eps ** 0.25
        res = optimize.minimize_scalar(obj, bounds=(_EPS, rng), method="bounded", options={"xatol": r_opt_tol})
        h_rot = float(res.x)
    return {"h": float(h_rot)}


# ---------------------------------------------------------------------------
def lpbwselect_mse_dpi(
    y, x, cluster, eval_pt, p, q, deriv, kernel, bwcheck, bwregul, vce, nnmatch,
    interior, weights=None,
):
    even = (p - deriv) % 2 == 0

    if kernel == "uni":
        C_c = 1.843
    elif kernel == "tri":
        C_c = 2.576
    elif kernel == "gau":
        C_c = 1.06
    else:
        C_c = 2.34

    x_iq  = np.quantile(x, 0.75) - np.quantile(x, 0.25)
    x_max, x_min = np.max(x), np.min(x)
    rng = x_max - x_min
    N = len(x)

    bw_max = max(abs(eval_pt - x_min), abs(eval_pt - x_max))
    c_bw = C_c * min(np.std(x, ddof=1), x_iq / 1.349) * N ** (-1.0 / 5)
    c_bw = min(c_bw, bw_max)

    dups = np.zeros(1, dtype=int)
    dupsid = np.zeros(1, dtype=int)
    if vce == "nn":
        order = np.argsort(x, kind="mergesort")
        x = x[order]; y = y[order]
        if cluster is not None:
            cluster = cluster[order]
        if weights is not None:
            weights = weights[order]
        dups, dupsid = build_dups(x)

    bw_min_val = None
    if bwcheck is not None:
        bw_min_val = np.sort(np.abs(x - eval_pt))[bwcheck - 1]
        c_bw = max(c_bw, bw_min_val)

    Cd1 = lprobust_bw(y, x, cluster, eval_pt, o=q + 1, nu=q + 1, oB=q + 2,
                     hV=c_bw, hB1=rng, hB2=rng, scale=0, vce=vce,
                     nnmatch=nnmatch, kernel=kernel, dups=dups, dupsid=dupsid,
                     weights=weights)
    if (not even) or interior:
        bw_mp2 = Cd1["bw"]
    else:
        def obj1(H):
            return abs(H ** (2 * (q + 1) + 2 - 2 * (q + 1)) * (Cd1["B1"] + H * Cd1["B2"]) ** 2
                       + Cd1["V"] / (N * H ** (1 + 2 * (q + 1))))
        bw_mp2 = float(optimize.minimize_scalar(obj1, bounds=(_EPS, rng), method="bounded").x)

    Cd2 = lprobust_bw(y, x, cluster, eval_pt, o=q + 2, nu=q + 2, oB=q + 3,
                     hV=c_bw, hB1=rng, hB2=rng, scale=0, vce=vce,
                     nnmatch=nnmatch, kernel=kernel, dups=dups, dupsid=dupsid,
                     weights=weights)
    if (not even) or interior:
        bw_mp3 = Cd2["bw"]
    else:
        def obj2(H):
            return abs(H ** (2 * (q + 2) + 2 - 2 * (q + 2)) * (Cd2["B1"] + H * Cd2["B2"]) ** 2
                       + Cd2["V"] / (N * H ** (1 + 2 * (q + 2))))
        bw_mp3 = float(optimize.minimize_scalar(obj2, bounds=(_EPS, rng), method="bounded").x)

    bw_mp2 = min(bw_mp2, bw_max)
    bw_mp3 = min(bw_mp3, bw_max)
    if bw_min_val is not None:
        bw_mp2 = max(bw_mp2, bw_min_val)
        bw_mp3 = max(bw_mp3, bw_min_val)

    Cb = lprobust_bw(y, x, cluster, eval_pt, o=q, nu=p + 1, oB=q + 1,
                    hV=c_bw, hB1=bw_mp2, hB2=bw_mp3, scale=bwregul, vce=vce,
                    nnmatch=nnmatch, kernel=kernel, dups=dups, dupsid=dupsid,
                    weights=weights)
    if (not even) or interior:
        b_mse = Cb["bw"]
    else:
        def objb(H):
            return abs(H ** (2 * q + 2 - 2 * (p + 1))
                       * (Cb["B1"] + H * Cb["B2"] + bwregul * Cb["R"]) ** 2
                       + Cb["V"] / (N * H ** (1 + 2 * (p + 1))))
        b_mse = float(optimize.minimize_scalar(objb, bounds=(_EPS, rng), method="bounded").x)

    b_mse = min(b_mse, bw_max)
    if bw_min_val is not None:
        b_mse = max(b_mse, bw_min_val)

    bw_mp1 = b_mse

    Ch = lprobust_bw(y, x, cluster, eval_pt, o=p, nu=deriv, oB=q,
                    hV=c_bw, hB1=bw_mp1, hB2=bw_mp2, scale=bwregul, vce=vce,
                    nnmatch=nnmatch, kernel=kernel, dups=dups, dupsid=dupsid,
                    weights=weights)
    if (not even) or interior:
        h_mse = Ch["bw"]
    else:
        def objh(H):
            return abs(H ** (2 * p + 2 - 2 * deriv)
                       * (Ch["B1"] + H * Ch["B2"] + bwregul * Ch["R"]) ** 2
                       + Ch["V"] / (N * H ** (1 + 2 * deriv)))
        h_mse = float(optimize.minimize_scalar(objh, bounds=(_EPS, rng), method="bounded").x)

    h_mse = min(h_mse, bw_max)
    if bw_min_val is not None:
        h_mse = max(h_mse, bw_min_val)

    if (not even) or interior:
        V_h = Ch["rV"] * Ch["V"]; B_h = Ch["rB"] * Ch["B1"] ** 2
        V_b = Cb["rV"] * Cb["V"]; B_b = Cb["rB"] * Cb["B1"] ** 2
    else:
        V_h = Ch["V"]; B_h = (Ch["B1"] + h_mse * Ch["B2"]) ** 2
        V_b = Cb["V"]; B_b = (Cb["B1"] + b_mse * Cb["B2"]) ** 2

    return {"h": h_mse, "b": b_mse, "V.h": V_h, "B.h": B_h, "V.b": V_b, "B.b": B_b}


# ---------------------------------------------------------------------------
def lpbwselect_imse_dpi(y, x, cluster, p, q, deriv, kernel, bwcheck, bwregul,
                       imsegrid, vce, nnmatch, interior, weights=None):
    N = len(x)
    x_max, x_min = np.max(x), np.min(x)
    rng = x_max - x_min
    eval_grid = np.linspace(x_min, x_max, imsegrid)
    even = (p - deriv) % 2 == 0
    V_h = np.empty(imsegrid); B_h = np.empty(imsegrid)
    V_b = np.empty(imsegrid); B_b = np.empty(imsegrid)
    for i, e in enumerate(eval_grid):
        r = lpbwselect_mse_dpi(y, x, cluster, e, p, q, deriv, kernel,
                             bwcheck, bwregul, vce, nnmatch, interior, weights)
        V_h[i] = r["V.h"]; B_h[i] = r["B.h"]
        V_b[i] = r["V.b"]; B_b[i] = r["B.b"]

    if (not even) or interior:
        b_imse = (np.mean(V_b) / (N * np.mean(B_b))) ** (1.0 / (2 * q + 3))
        h_imse = (np.mean(V_h) / (N * np.mean(B_h))) ** (1.0 / (2 * p + 3))
    else:
        def obj_b(H):
            return abs(H ** (2 * q + 2 - 2 * (p + 1)) * np.mean(B_b)
                       + np.mean(V_b) / (N * H ** (1 + 2 * (p + 1))))
        b_imse = float(optimize.minimize_scalar(obj_b, bounds=(_EPS, rng), method="bounded").x)
        def obj_h(H):
            return abs(H ** (2 * p + 2 - 2 * deriv) * np.mean(B_h)
                       + np.mean(V_h) / (N * H ** (1 + 2 * deriv)))
        h_imse = float(optimize.minimize_scalar(obj_h, bounds=(_EPS, rng), method="bounded").x)
    return {"h": h_imse, "b": b_imse}


# ---------------------------------------------------------------------------
def lpbwselect_ce_dpi(y, x, h, b, eval_pt, p, q, deriv, rho, kernel, vce,
                    nnmatch, interior, bwregul, weights=None):
    # R hard-codes both rho=1 and bwregul=0 inside lpbwselect.ce.dpi
    # (npfunctions.R:528-529). Mirror that to match R numerically.
    rho = 1.0
    bwregul = 0.0

    N = len(x)
    rng = float(np.max(x) - np.min(x))
    if weights is None:
        weights = np.ones(N)

    if vce == "nn":
        order = np.argsort(x, kind="mergesort")
        x = x[order]; y = y[order]; weights = weights[order]

    b = h / rho
    X_h = (x - eval_pt) / h
    X_b = (x - eval_pt) / b
    K_h = W_fun(X_h, kernel) * weights
    L_b = W_fun(X_b, kernel) * weights
    ind_h = K_h > 0
    ind_b = L_b > 0
    ind = ind_h | ind_b
    eN = int(ind.sum())
    eY = y[ind]; eX = x[ind]
    eX_h = X_h[ind]; eX_b = X_b[ind]
    eK_h = K_h[ind]; eL_b = L_b[ind]

    W_p = eK_h / h
    W_q = eL_b / b
    R_p_2 = eX_h[:, None] ** np.arange(p + 3)
    R_p   = R_p_2[:, : p + 1]
    R_q   = eX_b[:, None] ** np.arange(q + 1)

    L_p_1 = (R_p * W_p[:, None]).T @ (eX_h ** (p + 1)) / eN
    L_p_2 = (R_p * W_p[:, None]).T @ (eX_h ** (p + 2)) / eN
    L_p_3 = (R_p * W_p[:, None]).T @ (eX_h ** (p + 3)) / eN
    L_q_1 = (R_q * W_q[:, None]).T @ (eX_b ** (q + 1)) / eN
    L_q_2 = (R_q * W_q[:, None]).T @ (eX_b ** (q + 2)) / eN

    invG_p = eN * qrXXinv(R_p * np.sqrt(W_p)[:, None])
    invG_q = eN * qrXXinv(R_q * np.sqrt(W_q)[:, None])

    edups = np.zeros(1, dtype=int)
    edupsid = np.zeros(1, dtype=int)
    if vce == "nn":
        edups, edupsid = build_dups(eX)

    hii = np.zeros(1)
    predicts = np.zeros(eN)
    if vce in ("hc0", "hc1", "hc2", "hc3"):
        H_q = b ** (-np.arange(q + 1))
        beta_q = H_q * (invG_q @ (R_q * W_q[:, None]).T @ eY) / eN
        r_q = (eX - eval_pt)[:, None] ** np.arange(q + 1)
        predicts = r_q @ beta_q
        if vce in ("hc2", "hc3"):
            hii = np.einsum("ij,jk,ik->i",
                            R_p, invG_p, (R_p * W_p[:, None])) / eN

    res_q = lprobust_res(eX, eY, predicts.reshape(-1, 1), hii, vce, nnmatch,
                        edups, edupsid, q + 1)

    # Bias constants via high-order polynomial fit. R lm() uses pivoted QR
    # which is numerically stable on the Vandermonde matrix; scipy.lstsq
    # uses SVD and returns the minimum-norm solution, which differs when
    # x has large scale (e.g. cond ≈ 1e20 for x ≈ 300, k=4). Scale x first.
    k = p + 3
    s_x = float(np.max(np.abs(x))) if np.max(np.abs(x)) > 0 else 1.0
    u = x / s_x
    r_k = u[:, None] ** np.arange(k + 3)
    try:
        gamma, *_ = linalg.lstsq(r_k, y)
        beta = gamma / s_x ** np.arange(k + 3)
    except linalg.LinAlgError:
        beta = np.full(k + 3, np.nan)
    m_p_3 = (beta[p + 3] * math.factorial(p + 3)
             + beta[p + 4] * math.factorial(p + 4) * eval_pt
             + beta[p + 5] * math.factorial(p + 5) * eval_pt ** 2 / 2)
    r_k2 = u[:, None] ** np.arange(k + 2)
    try:
        gamma2, *_ = linalg.lstsq(r_k2, y)
        beta2 = gamma2 / s_x ** np.arange(k + 2)
    except linalg.LinAlgError:
        beta2 = np.full(k + 2, np.nan)
    m_p_2 = (beta2[p + 2] * math.factorial(p + 2)
             + beta2[p + 3] * math.factorial(p + 3) * eval_pt
             + beta2[p + 4] * math.factorial(p + 4) * eval_pt ** 2 / 2)

    # Use lprobust fallback if polynomial rank-deficient
    if not np.isfinite(m_p_3):
        from .lprobust import lprobust
        m_p_3 = lprobust(y, x, eval=np.array([eval_pt]), h=rng,
                         p=p + 4, deriv=p + 3, kernel=kernel, vce=vce,
                         masspoints="off").Estimate["tau.us"].iloc[0]
    if not np.isfinite(m_p_2):
        from .lprobust import lprobust
        m_p_2 = lprobust(y, x, eval=np.array([eval_pt]), h=rng,
                         p=p + 4, deriv=p + 2, kernel=kernel, vce=vce,
                         masspoints="off").Estimate["tau.us"].iloc[0]

    e_p_1 = np.zeros(q + 1); e_p_1[p + 1] = 1.0
    e_0   = np.zeros(p + 1); e_0[0] = 1.0

    # Pilot variance V.reg from a (p+2)-order WLS fit + sandwich. Mirrors
    # R lpbwselect.ce.dpi (npfunctions.R:620-626). Required to compute
    # the regularization term Reg used inside the optimize objective;
    # the previous Python code dropped V.reg entirely (Reg=0).
    V_reg = 0.0
    if bwregul > 0:
        o_reg = p + 2
        # Use (eX-eval_pt)^j directly to avoid catastrophic cancellation
        R_reg = (eX[:, None] - eval_pt) ** np.arange(o_reg + 1)
        invG_reg = qrXXinv(R_reg * np.sqrt(eK_h)[:, None])
        meat_reg = lprobust_vce(R_reg * eK_h[:, None], res_q, None)
        V_reg = float((invG_reg @ meat_reg @ invG_reg)[o_reg, o_reg])

    q_terms = lpbwce(
        eY, eX, eK_h, eL_b, res_q.ravel(), eval_pt, p, q, h, b, deriv,
        math.factorial(deriv),
    )
    q1_rbc = q_terms["q1rbc"]
    q2_rbc = q_terms["q2rbc"]
    q3_rbc = q_terms["q3rbc"]

    if interior:
        eta_bc1 = (e_0 @ invG_p) @ (
            (m_p_2 / math.factorial(p + 2)) * L_p_2 / h
            + (m_p_3 / math.factorial(p + 3)) * L_p_3
        )
        eta_bc2 = rho ** (-2) * b ** (q - p - 1) * (
            (e_0 @ invG_p) @ L_p_1 * (e_p_1 @ invG_q) @ (
                (m_p_2 / math.factorial(p + 2)) * L_q_1 / b
                + (m_p_3 / math.factorial(p + 3)) * L_q_2
            )
        )
        eta_bc = float(eta_bc1 - eta_bc2)
        # R: Reg = 3*(e_0' invG_p (L_p_2/h)^2) * V_reg  (npfunctions.R:638)
        Reg = float(3.0 * ((e_0 @ invG_p) @ (L_p_2 / h) ** 2) * V_reg)
        def obj(H):
            return abs(H ** (-1) * q1_rbc
                       + H ** (1 + 2 * (p + 3)) * (eta_bc ** 2 + bwregul * Reg) * q2_rbc
                       + H ** (p + 3) * (eta_bc + bwregul * Reg) * q3_rbc)
        r_opt_tol = np.finfo(float).eps ** 0.25
        res = optimize.minimize_scalar(obj, bounds=(_EPS, rng), method="bounded", options={"xatol": r_opt_tol})
        h_ce = float(res.x) * N ** (-1.0 / (p + 4))
    else:
        eta_bc_1 = (e_0 @ invG_p) @ (
            L_p_2 - rho ** (-1) * L_p_1 * ((e_p_1 @ invG_q) @ L_q_1)
        ) * (m_p_2 / math.factorial(p + 2))
        eta_bc_2 = (e_0 @ invG_p) @ (
            L_p_3 - rho ** (-2) * L_p_1 * ((e_p_1 @ invG_q) @ L_q_2)
        ) * (m_p_3 / math.factorial(p + 3))
        # R: Reg = 3 * ((e_0' invG_p) (L_p_2 - rho^-1 L_p_1 (e_p_1' invG_q) L_q_1))^2 * V_reg
        # (npfunctions.R:648).
        Reg = float(3.0 * ((e_0 @ invG_p) @ (
            L_p_2 - rho ** (-1) * L_p_1 * ((e_p_1 @ invG_q) @ L_q_1)
        )) ** 2 * V_reg)
        from scipy.stats import norm
        phi_z = norm.pdf(1.96)
        E1 = q1_rbc * phi_z
        E2 = float(eta_bc_1)
        E3 = float(eta_bc_2)
        E4 = q3_rbc * phi_z * E2
        E5 = q3_rbc * phi_z * E3
        def obj(H):
            return abs(
                E1 / (N * H)
                + N * H ** (2 * p + 5) * q2_rbc * phi_z * (E2 + H * E3 + bwregul * Reg) ** 2
                + H ** (p + 2) * (E4 + H * E5 + bwregul * Reg)
            )
        r_opt_tol = np.finfo(float).eps ** 0.25
        res = optimize.minimize_scalar(obj, bounds=(_EPS, rng), method="bounded", options={"xatol": r_opt_tol})
        h_ce = float(res.x)

    return {"h": h_ce}
