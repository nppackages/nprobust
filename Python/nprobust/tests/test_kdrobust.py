"""kdrobust numerical checks."""
import numpy as np
import pytest
from scipy.stats import norm

from nprobust import kdrobust
from .conftest import w_kernel


@pytest.mark.parametrize("kernel", ["epa", "uni"])
@pytest.mark.parametrize("eval_pt", [-0.5, 0.0, 0.7])
def test_kdrobust_matches_direct_density(kernel, eval_pt):
    rng = np.random.default_rng(3)
    x = rng.normal(size=2000)
    h = 0.35
    kd = kdrobust(x, eval=np.array([eval_pt]), h=h, kernel=kernel)
    ref = np.mean(w_kernel((x - eval_pt) / h, kernel)) / h
    assert kd.Estimate["tau.us"].iloc[0] == pytest.approx(ref, abs=1e-12)


def test_kdrobust_rejects_unsupported_kernel():
    rng = np.random.default_rng(3)
    x = rng.normal(size=200)
    with pytest.raises(ValueError, match="kernel"):
        kdrobust(x, kernel="tri")
    with pytest.raises(ValueError, match="kernel"):
        kdrobust(x, kernel="gau")


def test_kdrobust_recovers_normal_density_at_zero():
    rng = np.random.default_rng(1)
    x = rng.normal(size=30000)
    kd = kdrobust(x, eval=np.array([0.0]), h=0.3, kernel="epa")
    assert kd.Estimate["tau.us"].iloc[0] == pytest.approx(norm.pdf(0), abs=0.01)
    assert kd.Estimate["tau.bc"].iloc[0] == pytest.approx(norm.pdf(0), abs=0.01)


def test_kdrobust_integrates_to_one():
    rng = np.random.default_rng(1)
    x = rng.normal(size=20000)
    step = 0.1
    grid = np.arange(-4, 4 + step / 2, step)
    kd = kdrobust(x, eval=grid, h=0.4, kernel="epa")
    integral = kd.Estimate["tau.us"].sum() * step
    assert integral == pytest.approx(1.0, abs=0.02)
