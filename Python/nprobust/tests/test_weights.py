import numpy as np
import pytest

from nprobust import lprobust
from .conftest import w_kernel


def test_weights_equal_wls(fixture_data):
    x, y = fixture_data["x"], fixture_data["y"]
    rng = np.random.default_rng(3)
    w = rng.uniform(0.2, 2, size=len(x))
    h = 0.3
    m = lprobust(y, x, eval=np.array([0.0]), h=h, p=1, kernel="uni",
                 vce="hc0", weights=w)

    u = (x - 0.0) / h
    kw = w_kernel(u, "uni")
    comb = kw * w
    sub = comb > 0
    xc = x[sub]; ys = y[sub]; ws = comb[sub]
    X = np.column_stack([np.ones_like(xc), xc])
    Xs = X * np.sqrt(ws)[:, None]
    beta = np.linalg.lstsq(Xs, ys * np.sqrt(ws), rcond=None)[0]
    assert m.Estimate["tau.us"].iloc[0] == pytest.approx(beta[0], abs=1e-10)


def test_weights_ones_matches_unweighted(fixture_data):
    x, y = fixture_data["x"], fixture_data["y"][:400]
    x = x[:400]
    m0 = lprobust(y, x, eval=np.array([0.0]), h=0.3, p=1, kernel="epa", vce="hc0")
    m1 = lprobust(y, x, eval=np.array([0.0]), h=0.3, p=1, kernel="epa", vce="hc0",
                  weights=np.ones_like(x))
    assert m0.Estimate["tau.us"].iloc[0] == pytest.approx(
        m1.Estimate["tau.us"].iloc[0], abs=1e-12)


def test_zero_weights_match_subset(fixture_data):
    x, y = fixture_data["x"][:400], fixture_data["y"][:400]
    w = np.ones_like(x); w[x > 0.5] = 0
    mw = lprobust(y, x, eval=np.array([0.3]), h=0.3, p=1, kernel="uni",
                  vce="hc0", weights=w)
    ms = lprobust(y, x, eval=np.array([0.3]), h=0.3, p=1, kernel="uni",
                  vce="hc0", subset=(x <= 0.5))
    assert mw.Estimate["tau.us"].iloc[0] == pytest.approx(
        ms.Estimate["tau.us"].iloc[0], abs=1e-12)
