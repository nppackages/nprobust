"""Tidy/glance outputs + plot smoke tests."""
import numpy as np
import pytest

from nprobust import lprobust, kdrobust, lpbwselect


def test_lprobust_tidy_glance(fixture_data):
    m = lprobust(fixture_data["y"], fixture_data["x"],
                 eval=np.array([0.0, 0.5]))
    td = m.tidy()
    expected_cols = ["eval", "estimate", "std.error", "tau.bc", "se.rb",
                     "h", "b", "n.eff", "conf.low", "conf.high"]
    for c in expected_cols:
        assert c in td.columns
    assert len(td) == 2
    assert len(m.glance()) == 1


def test_plot_returns_figure(fixture_data):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    from nprobust import plot_lprobust, plot_kdrobust
    m = lprobust(fixture_data["y"], fixture_data["x"], neval=10)
    fig = plot_lprobust(m)
    assert fig is not None
    kd = kdrobust(fixture_data["x"], neval=10)
    fig2 = plot_kdrobust(kd)
    assert fig2 is not None
