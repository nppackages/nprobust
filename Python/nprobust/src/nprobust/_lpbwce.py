"""Coverage-error DPI helper.

Pure NumPy port of the (already pure-R) :func:`lpbwce` in
``npfunctions.R``. See the R source for the derivation.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy import linalg


def lpbwce(
    y: np.ndarray,
    x: np.ndarray,
    K: np.ndarray,
    L: np.ndarray,
    res: np.ndarray,
    c: float,
    p: int,
    q: int,
    h: float,
    b: float,
    deriv: int,
    fact: int,
) -> Dict[str, float]:
    """Return the three Edgeworth-expansion terms q1rbc, q2rbc, q3rbc."""
    y = np.asarray(y, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    K = np.asarray(K, dtype=float).ravel()
    L = np.asarray(L, dtype=float).ravel()
    res = np.asarray(res, dtype=float).ravel()

    N = len(y)
    rho = h / b

    Xh = (x - c) / h
    Xb = (x - c) / b

    Rp = Xh[:, None] ** np.arange(p + 1)
    Rq = Xb[:, None] ** np.arange(q + 1)
    Wp = K / h
    Wq = L / b

    Lp1 = (Rp * Wp[:, None]).T @ (Xh ** (p + 1)) / N          # (p+1,)
    Gp  = (Rp * Wp[:, None]).T @ Rp / N
    Gq  = (Rq * Wq[:, None]).T @ Rq / N
    iGp = linalg.inv(Gp)
    iGq = linalg.inv(Gq)

    ep1 = np.zeros(q + 1); ep1[p + 1] = 1.0
    e0  = np.zeros(p + 1); e0[deriv] = float(fact)

    # lbc0 (length N)
    lus0 = (e0 @ iGp) @ (Rp * K[:, None]).T
    scale_bc = float((e0 @ iGp) @ Lp1)
    tail_bc  = (ep1 @ iGq) @ (Rq * L[:, None]).T
    lbc0 = lus0 - rho ** (p + 1) * scale_bc * tail_bc

    vx = res ** 2
    s2 = float(np.sum(lbc0 ** 2 * vx) / (N * h))

    Krrp   = (Rp * K[:, None]).T @ Rp
    Lrrq   = (Rq * L[:, None]).T @ Rq
    Krxip  = Rp.T @ (K * Xh ** (p + 1))
    sumKRp = Rp.T @ K
    Sxp1   = np.sum(Xh ** (p + 1))
    Krxp   = sumKRp * Sxp1 - Krxip

    EKrrp      = Krrp / N
    EKrxip_vec = Krxip / N
    EKrxp_vec  = Krxp / (N * (N - 1))
    ELrrq      = Lrrq / N

    a_row = float(fact) * iGp[deriv, :]
    u_a   = a_row.copy()
    v_a   = a_row @ EKrrp @ iGp

    RpiGp  = Rp @ iGp
    RqiGq  = Rq @ iGq
    quadRp = np.einsum("ij,ij->i", RpiGp, Rp)
    quadRq = np.einsum("ij,ij->i", RqiGq, Rq)

    Rp_va = Rp @ v_a
    Rp_ua = Rp @ u_a

    ## Single-index sums
    q1  = np.sum(lbc0 ** 3 * res ** 3)
    q3  = np.sum(lbc0 ** 4 * (res ** 4 - vx ** 2))
    q3a = q1
    q8  = np.sum((lbc0 * res) ** 4)
    q9  = np.sum((lbc0 ** 2 * vx - h * s2) * (lbc0 * res) ** 2)
    q12 = np.sum((lbc0 ** 2 * vx - h * s2) ** 2)
    q4  = np.sum(lbc0 ** 2 * L * quadRq * res ** 2)

    q5a = (lbc0 ** 3 * res ** 2) @ RqiGq        # (q+1,)
    q5b = (lbc0 * res ** 2 * L)   @ Rq          # (q+1,)
    q7a = (lbc0 * res ** 2 * L)   @ RqiGq       # (q+1,)
    q7b = (Rq * lbc0[:, None]).T @ (Rq * lbc0[:, None]) @ iGq
    q7c = q5b.copy()

    ## lbc1 at (i, i) for q2
    lus1_diag = K * Rp_va - K ** 2 * Rp_ua * quadRp

    C1_scalar = float(a_row @ EKrrp @ iGp @ Lp1)
    v_iGq_ep1 = iGq @ ep1
    C2_vec    = Rq @ v_iGq_ep1
    v_lp      = iGp @ Lp1
    D2_vec    = Rp @ v_lp

    T1_diag = L * C2_vec * (C1_scalar - K * Rp_ua * D2_vec)

    dot_a_krxip = float(u_a @ EKrxip_vec)
    T2_diag = L * C2_vec * (K * Xh ** (p + 1) * Rp_ua - dot_a_krxip)

    a_Lp1 = float(a_row @ Lp1)
    u_T3  = a_Lp1 * (ep1 @ iGq)
    u_T3_ELrrq_iGq = u_T3 @ ELrrq @ iGq
    Rq_uT3         = Rq @ u_T3
    Rq_uT3_ELrrq   = Rq @ u_T3_ELrrq_iGq

    T3_diag = L * Rq_uT3_ELrrq - L ** 2 * Rq_uT3 * quadRq

    lbc1_diag = lus1_diag - rho ** (p + 1) * (T1_diag + T2_diag + T3_diag)
    q2 = float(np.sum(lbc1_diag * lbc0 * res ** 2))

    ## Double-sum terms (i != j): N x N matrices.
    dot_a_krxp = float(u_a @ EKrxp_vec)

    M_iGp = Rp @ iGp @ Rp.T
    M_iGq = Rq @ iGq @ Rq.T

    A_mat     = np.broadcast_to(K * Rp_va, (N, N)).T            # [i,j] = K[i]*Rp_va[i]
    KjRp_ua_M = np.broadcast_to(K * Rp_ua, (N, N)) * M_iGp.T    # [i,j] = K[j]*Rp_ua[j]*M_iGp[j,i]
    Ki_mat    = np.broadcast_to(K, (N, N)).T
    B_mat     = Ki_mat * KjRp_ua_M

    Li_C2_mat = np.broadcast_to(L * C2_vec, (N, N)).T
    KjRp_uaD2 = np.broadcast_to(K * Rp_ua * D2_vec, (N, N))
    T1_mat    = Li_C2_mat * (C1_scalar - KjRp_uaD2)

    Xhpp1     = Xh ** (p + 1)
    first_T2  = np.outer(Xhpp1, K * Rp_ua)
    T2_mat    = Li_C2_mat * (first_T2 - dot_a_krxp)

    first_T3  = np.broadcast_to(L * Rq_uT3_ELrrq, (N, N)).T
    Lj_Rq_uT3 = np.broadcast_to(L * Rq_uT3, (N, N))
    Li_mat    = np.broadcast_to(L, (N, N)).T
    T3_mat    = first_T3 - Li_mat * Lj_Rq_uT3 * M_iGq.T

    lbc1_mat = A_mat - B_mat - rho ** (p + 1) * (T1_mat + T2_mat + T3_mat)

    ressq = res ** 2

    q6mat = np.outer(lbc0 ** 2, L ** 2 * ressq) * M_iGq ** 2
    np.fill_diagonal(q6mat, 0.0)
    q6 = float(q6mat.sum())

    q10mat = lbc1_mat * np.outer(lbc0 * ressq, lbc0 ** 2 * ressq)
    np.fill_diagonal(q10mat, 0.0)
    q10 = float(q10mat.sum())

    q11mat = lbc1_mat * np.outer(lbc0 * ressq, lbc0 ** 2 * ressq - h * s2)
    np.fill_diagonal(q11mat, 0.0)
    q11 = float(q11mat.sum())

    ## Expectations
    Eq1  = (q1 / (N * h)) ** 2
    Eq2  = q2 / (N * h)
    Eq3  = q3 / (N * h)
    Eq4  = q4 / (N * h)
    Eq5  = float(q5a @ q5b) / (N * h) ** 2
    Eq6  = q6 / (N * (N - 1) * h ** 2)
    Eq7  = float(q7a @ q7b @ q7c) / (N * h) ** 3
    Eq8  = q8 / (N * h)
    Eq9  = q9 / (N * h)
    Eq10 = q10 / (N * (N - 1) * h ** 2)
    Eq11 = q11 / (N * (N - 1) * h ** 2)
    Eq12 = q12 / (N * h)

    z  = 1.959964
    pz = 0.05844507

    q1bc = pz * (
        Eq1  * (z ** 3 / 3 + 7 * z / 4 + s2 * z * (z ** 2 - 3) / 4) / s2 ** 3
      + Eq2  * (-z * (z ** 2 - 3) / 2)                              / s2
      + Eq3  * (z * (z ** 2 - 3) / 8)                               / s2 ** 2
      - Eq4  * (z * (z ** 2 - 1) / 2)                               / s2
      - Eq5  * (z * (z ** 2 - 1))                                   / s2 ** 2
      + Eq6  * (z * (z ** 2 - 1) / 4)                               / s2
      + Eq7  * (z * (z ** 2 - 1) / 2)                               / s2 ** 2
      + Eq8  * (-z * (z ** 2 - 3) / 24)                             / s2 ** 2
      + Eq9  * (z * (z ** 2 - 1) / 4)                               / s2 ** 2
      + Eq10 * (z * (z ** 2 - 3))                                   / s2 ** 2
      + Eq11 * (-z)                                                 / s2 ** 2
      + Eq12 * (-z * (z ** 2 + 1) / 8)                              / s2 ** 2
    )
    q2bc = -pz * z / (2 * s2)
    Eq3a = q3a / (N * h)
    q3bc = pz * Eq3a / s2 ** 2 * (z ** 3 / 3)

    return {
        "q1rbc": 2 * q1bc / pz,
        "q2rbc": 2 * q2bc / pz,
        "q3rbc": 2 * q3bc / pz,
    }
