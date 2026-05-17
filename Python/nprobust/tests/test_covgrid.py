"""covgrid=True should produce symmetric covariance matrices whose diagonal
equals se.us^2 / se.rb^2 at each eval point, and whose off-diagonals are
finite real numbers (not 0 placeholders)."""
import numpy as np

from nprobust import lprobust


def test_covgrid_diagonal_matches_squared_se(fixture_data):
    eval_pts = np.linspace(0.2, 0.8, 5)
    m = lprobust(fixture_data["y"], fixture_data["x"], eval=eval_pts, covgrid=True)

    assert m.cov_us.shape == (5, 5)
    assert m.cov_rb.shape == (5, 5)
    assert np.allclose(m.cov_us, m.cov_us.T)
    assert np.allclose(m.cov_rb, m.cov_rb.T)

    np.testing.assert_allclose(
        np.diag(m.cov_us), m.Estimate["se.us"].values ** 2, rtol=1e-10)
    np.testing.assert_allclose(
        np.diag(m.cov_rb), m.Estimate["se.rb"].values ** 2, rtol=1e-10)


def test_covgrid_off_diagonals_are_non_zero(fixture_data):
    eval_pts = np.linspace(0.2, 0.8, 5)
    m = lprobust(fixture_data["y"], fixture_data["x"], eval=eval_pts, covgrid=True)

    off = m.cov_us - np.diag(np.diag(m.cov_us))
    assert np.any(np.abs(off) > 1e-10), (
        "off-diagonal cov_us entries are zero; covgrid not implemented")


def test_covgrid_is_positive_semidefinite(fixture_data):
    eval_pts = np.linspace(0.2, 0.8, 5)
    m = lprobust(fixture_data["y"], fixture_data["x"], eval=eval_pts, covgrid=True)

    # Robust-CI covariance should be PSD up to numerical tolerance.
    eig = np.linalg.eigvalsh((m.cov_rb + m.cov_rb.T) / 2)
    assert eig.min() > -1e-10


def test_covgrid_scaling_correct_for_deriv_ge_2(fixture_data):
    """Pre-fix: V_us_i carried factorial(deriv)**2, so cross-product gave
    factorial(deriv)**4 instead of factorial(deriv)**2 -- diag(cov_us) was
    4x too large at deriv=2 (36x at deriv=3). Latent because all other
    tests in this file use deriv=0 where factorial(0)==1."""
    eval_pts = np.array([0.2, 0.5, 0.7])
    for d in (0, 1, 2, 3):
        m = lprobust(
            fixture_data["y"], fixture_data["x"],
            p=d + 1, deriv=d, eval=eval_pts, covgrid=True,
        )
        np.testing.assert_allclose(
            np.diag(m.cov_us), m.Estimate["se.us"].values ** 2,
            rtol=1e-8,
            err_msg=f"deriv={d}: diag(cov_us) != se.us**2",
        )
        np.testing.assert_allclose(
            np.diag(m.cov_rb), m.Estimate["se.rb"].values ** 2,
            rtol=1e-8,
            err_msg=f"deriv={d}: diag(cov_rb) != se.rb**2",
        )
