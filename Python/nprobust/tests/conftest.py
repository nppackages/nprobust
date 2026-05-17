import numpy as np
import pytest


@pytest.fixture
def fixture_data():
    rng = np.random.default_rng(7)
    n = 1000
    x = np.sort(rng.uniform(-1, 1, n))
    y = 1 + 2 * x + 3 * x ** 2 - 0.5 * x ** 3 + rng.normal(0, 0.3, n)
    return {"x": x, "y": y, "n": n}


def w_kernel(u, kernel):
    u = np.asarray(u)
    if kernel == "uni": return 0.5 * (np.abs(u) <= 1)
    if kernel == "tri": return (1 - np.abs(u)) * (np.abs(u) <= 1)
    if kernel == "epa": return 0.75 * (1 - u ** 2) * (np.abs(u) <= 1)
    raise ValueError(kernel)
