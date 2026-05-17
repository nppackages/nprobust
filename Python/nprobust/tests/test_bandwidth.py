import numpy as np
import pytest

from nprobust import lprobust, lpbwselect


@pytest.mark.parametrize("bw", ["mse-dpi", "mse-rot", "imse-dpi",
                                 "imse-rot", "ce-dpi", "ce-rot"])
def test_bw_selectors_produce_positive(fixture_data, bw):
    b = lpbwselect(fixture_data["y"], fixture_data["x"],
                   eval=np.array([0.0]), bwselect=bw)
    assert np.all(np.isfinite(b.bws["h"])) and np.all(b.bws["h"] > 0)
    assert np.all(np.isfinite(b.bws["b"])) and np.all(b.bws["b"] > 0)


def test_interior_flag_runs(fixture_data):
    b1 = lpbwselect(fixture_data["y"], fixture_data["x"],
                    eval=np.array([0.0]), bwselect="mse-dpi", interior=True)
    b2 = lpbwselect(fixture_data["y"], fixture_data["x"],
                    eval=np.array([0.0]), bwselect="mse-dpi", interior=False)
    assert np.isfinite(b1.bws["h"].iloc[0])
    assert np.isfinite(b2.bws["h"].iloc[0])


def test_bwcheck_larger_than_n_warns(fixture_data):
    x = fixture_data["x"][:10]; y = fixture_data["y"][:10]
    with pytest.warns(UserWarning, match="bwcheck"):
        lprobust(y, x, eval=np.array([0.5]), bwcheck=50)


def test_rho_zero_with_h_no_b_errors(fixture_data):
    with pytest.raises(ValueError, match="b must also be provided"):
        lprobust(fixture_data["y"], fixture_data["x"],
                 eval=np.array([0.5]), h=0.1, rho=0)
