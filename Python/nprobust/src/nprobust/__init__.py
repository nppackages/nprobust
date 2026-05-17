"""nprobust: nonparametric robust estimation and inference.

Python implementation of estimation, inference, bandwidth selection, and
graphical procedures for kernel density and local polynomial regression
methods, including robust bias-corrected confidence intervals (Calonico,
Cattaneo and Farrell, JASA 2018; JSS 2019; Bernoulli 2022).

Public API:
    lprobust       - local polynomial point estimation + RBC inference
    lpbwselect     - bandwidth selection for local polynomial regression
    kdrobust       - kernel density point estimation + RBC inference
    kdbwselect     - bandwidth selection for kernel density estimation
    plot_lprobust  - matplotlib plot of lprobust results
    plot_kdrobust  - matplotlib plot of kdrobust results
"""
from .lprobust import lprobust, LprobustResult
from .lpbwselect import lpbwselect, LpbwselectResult
from .kdrobust import kdrobust, KdrobustResult
from .kdbwselect import kdbwselect, KdbwselectResult
from .plot import plot_lprobust, plot_kdrobust

__version__ = "1.0.0"

__all__ = [
    "lprobust", "LprobustResult",
    "lpbwselect", "LpbwselectResult",
    "kdrobust", "KdrobustResult",
    "kdbwselect", "KdbwselectResult",
    "plot_lprobust", "plot_kdrobust",
    "__version__",
]
