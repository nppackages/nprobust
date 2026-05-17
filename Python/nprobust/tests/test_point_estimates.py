"""Cross-check lprobust point estimates against WLS."""
import math

import numpy as np
import pytest

from nprobust import lprobust
from .conftest import w_kernel


@pytest.mark.parametrize("kernel", ["uni", "tri", "epa"])
@pytest.mark.parametrize("p", [1, 2, 3])
@pytest.mark.parametrize("eval_pt", [-0.3, 0.0, 0.4])
def test_point_estimate_matches_wls(fixture_data, kernel, p, eval_pt):
    x, y = fixture_data["x"], fixture_data["y"]
    h = 0.25
    m = lprobust(y, x, eval=np.array([eval_pt]), h=h, p=p,
                 kernel=kernel, vce="hc0")
    u = (x - eval_pt) / h
    w = w_kernel(u, kernel)
    sub = w > 0
    xc = x[sub] - eval_pt
    ys = y[sub]
    ws = w[sub]
    X = np.column_stack([np.ones_like(xc)] + [xc ** k for k in range(1, p + 1)])
    Xs = X * np.sqrt(ws)[:, None]
    beta = np.linalg.lstsq(Xs, ys * np.sqrt(ws), rcond=None)[0]
    assert m.Estimate["tau.us"].iloc[0] == pytest.approx(beta[0], abs=1e-10)


@pytest.mark.parametrize("deriv", [0, 1, 2])
def test_derivative_matches_factorial_coef(fixture_data, deriv):
    x, y = fixture_data["x"], fixture_data["y"]
    p = deriv + 1
    h = 0.3
    m = lprobust(y, x, eval=np.array([0.0]), h=h, p=p, deriv=deriv,
                 kernel="uni", vce="hc0")
    sub = np.abs(x) <= h
    xc = x[sub]; ys = y[sub]
    X = np.column_stack([xc ** k for k in range(p + 1)])
    beta = np.linalg.lstsq(X, ys, rcond=None)[0]
    assert m.Estimate["tau.us"].iloc[0] == pytest.approx(
        math.factorial(deriv) * beta[deriv], abs=1e-10)
