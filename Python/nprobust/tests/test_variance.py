"""HC0/HC1/HC2/HC3 and CR0/CR1/CR2/CR3 SE match the manual sandwich."""
import numpy as np
import pytest
from scipy import linalg

from nprobust import lprobust


def _sandwich_se(y, x, eval_pt, h, vce, p=1):
    u = (x - eval_pt) / h
    sub = np.abs(u) <= 1
    xc = x[sub] - eval_pt
    ys = y[sub]
    W = np.full(sub.sum(), 0.5 / h)
    X = np.column_stack([np.ones_like(xc), xc])
    XWX_inv = linalg.inv((X * np.sqrt(W)[:, None]).T @ (X * np.sqrt(W)[:, None]))
    beta = XWX_inv @ (X * W[:, None]).T @ ys
    resid = ys - X @ beta
    H = X @ XWX_inv @ (X * W[:, None]).T
    hii = np.diag(H)
    n, k = X.shape
    if vce == "hc0":
        adj = np.ones(n)
    elif vce == "hc1":
        adj = np.full(n, np.sqrt(n / max(n - (p + 1), 1)))
    elif vce == "hc2":
        adj = np.sqrt(1 / (1 - hii))
    elif vce == "hc3":
        adj = 1 / (1 - hii)
    e = adj * resid
    B = (X * W[:, None] * e[:, None]).T @ (X * W[:, None] * e[:, None])
    vc = XWX_inv @ B @ XWX_inv
    return float(np.sqrt(vc[0, 0]))


@pytest.mark.parametrize("vce", ["hc0", "hc1", "hc2", "hc3"])
def test_hc_se(fixture_data, vce):
    x, y = fixture_data["x"], fixture_data["y"]
    m = lprobust(y, x, eval=np.array([0.0]), h=0.3, p=1, kernel="uni", vce=vce)
    assert m.Estimate["se.us"].iloc[0] == pytest.approx(
        _sandwich_se(y, x, 0.0, 0.3, vce), rel=1e-10, abs=1e-12)


def _cluster_se(y, x, eval_pt, h, cluster, cr_type, p=1):
    u = (x - eval_pt) / h
    sub = np.abs(u) <= 1
    xc = x[sub]; ys = y[sub]; cs = cluster[sub]
    W = np.full(sub.sum(), 0.5 / h)
    sqrtW = np.sqrt(W)
    X = np.column_stack([np.ones_like(xc), xc])
    X_std = X * sqrtW[:, None]
    XWX_inv = linalg.inv(X_std.T @ X_std)
    beta = XWX_inv @ (X * W[:, None]).T @ ys
    resid = ys - X @ beta
    r_std = resid * sqrtW

    clusters = np.unique(cs)
    G = len(clusters); n = len(ys); k = X.shape[1]
    M = np.zeros((k, k))
    for cl in clusters:
        idx = cs == cl
        Xg = X_std[idx]
        rg = r_std[idx]
        if cr_type in ("CR0", "CR1"):
            r_adj = rg
        else:
            Hgg = Xg @ XWX_inv @ Xg.T
            Hgg = 0.5 * (Hgg + Hgg.T)
            I_H = np.eye(Xg.shape[0]) - Hgg
            if cr_type == "CR2":
                vals, vecs = linalg.eigh(I_H)
                vals = np.maximum(vals, 1e-12)
                r_adj = vecs @ ((vecs.T @ rg) / np.sqrt(vals))
            else:  # CR3
                r_adj = linalg.solve(I_H, rg, assume_a="sym")
        score = Xg.T @ r_adj
        M += np.outer(score, score)
    if cr_type == "CR0":   mult = 1.0
    elif cr_type == "CR1": mult = ((n - 1) / max(n - k, 1)) * (G / max(G - 1, 1))
    elif cr_type == "CR2": mult = 1.0
    elif cr_type == "CR3": mult = (G - 1) / max(G, 1)
    vc = XWX_inv @ (mult * M) @ XWX_inv
    return float(np.sqrt(vc[0, 0]))


@pytest.mark.parametrize("vce,cr", [("cr1", "CR1"),
                                     ("cr2", "CR2"), ("cr3", "CR3")])
def test_cluster_se(fixture_data, vce, cr):
    x, y = fixture_data["x"], fixture_data["y"]
    rng = np.random.default_rng(11)
    cl = rng.integers(0, 15, size=fixture_data["n"])
    m = lprobust(y, x, eval=np.array([0.0]), h=0.3, p=1, kernel="uni",
                 vce=vce, cluster=cl)
    assert m.Estimate["se.us"].iloc[0] == pytest.approx(
        _cluster_se(y, x, 0.0, 0.3, cl, cr), rel=1e-9, abs=1e-11)


def test_cr1_with_singleton_clusters_matches_manual(fixture_data):
    """With G == n singleton clusters, CR1's df multiplier
    ((n-1)/(n-k)) * (G/(G-1)) collapses to n/(n-k); the cluster meat
    just sums per-obs scores (no within-cluster correlation), so the
    SE matches the manual CR1 computation at G=n."""
    x, y = fixture_data["x"], fixture_data["y"]
    cl = np.arange(len(x))
    m_cr = lprobust(y, x, eval=np.array([0.0]), h=0.3, p=1, kernel="uni",
                    vce="cr1", cluster=cl)
    se_ref = _cluster_se(y, x, 0.0, 0.3, cl, "CR1")
    assert m_cr.Estimate["se.us"].iloc[0] == pytest.approx(se_ref, abs=1e-10)
