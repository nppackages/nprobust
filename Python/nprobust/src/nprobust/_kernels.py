"""Kernel weight functions.

Python port of :func:`W.fun` and :func:`kd.K.fun` from the R package.
"""
from __future__ import annotations

import math
from typing import Callable, Tuple

import numpy as np
from scipy import integrate, stats


_NORM_PDF = stats.norm.pdf


def W_fun(u: np.ndarray, kernel: str) -> np.ndarray:
    """Kernel weight on the scaled argument ``u``.

    Matches :func:`W.fun` in the R package.

    Supported kernels: ``"epa"``, ``"uni"``, ``"tri"``, ``"gau"``.
    """
    u = np.asarray(u, dtype=float)
    if kernel == "epa":
        return 0.75 * (1.0 - u ** 2) * (np.abs(u) <= 1)
    if kernel == "uni":
        return 0.5 * (np.abs(u) <= 1)
    if kernel == "tri":
        return (1.0 - np.abs(u)) * (np.abs(u) <= 1)
    if kernel == "gau":
        return _NORM_PDF(u)
    raise ValueError(f"kernel {kernel!r} not recognised. Supported: epa, uni, tri, gau.")


def _epa_v6_r0(u):
    return (np.abs(u) <= 1) * (35.0 / 256.0) * (-99 * u ** 6 + 189 * u ** 4 - 105 * u ** 2 + 15)


def _epa_v6_r2(u):
    return (np.abs(u) <= 1) * (315.0 / 64.0) * (77 * u ** 6 - 135 * u ** 4 + 63 * u ** 2 - 5)


def _kd_k_impl(v: int, r: int, kernel: str) -> Callable[[np.ndarray], np.ndarray]:
    """Return the univariate kernel k(u) for (v, r, kernel)."""
    if v == 2:
        if kernel == "gau":
            if r == 0:
                return lambda u: _NORM_PDF(u)
            if r == 2:
                return lambda u: (u ** 2 - 1) * _NORM_PDF(u)
            if r == 4:
                return lambda u: (u ** 4 - 6 * u ** 2 + 3) * _NORM_PDF(u)
        if kernel == "uni" and r == 0:
            return lambda u: 0.5 * (np.abs(u) <= 1)
        if kernel == "epa" and r == 0:
            return lambda u: 0.75 * (1.0 - u ** 2) * (np.abs(u) <= 1)
    if v == 4:
        if kernel == "uni":
            if r == 0:
                return lambda u: (np.abs(u) <= 1) * 3 * (-5 * u ** 2 + 3) / 8
            if r == 2:
                return lambda u: (np.abs(u) <= 1) * 15 * (3 * u ** 2 - 1) / 4
        if kernel == "epa":
            if r == 0:
                return lambda u: (np.abs(u) <= 1) * (15.0 / 32.0) * (7 * u ** 4 - 10 * u ** 2 + 3)
            if r == 2:
                return lambda u: (np.abs(u) <= 1) * (105.0 / 16.0) * (6 * u ** 2 - 5 * u ** 4 - 1)
    if v == 6:
        if kernel == "uni":
            if r == 0:
                return lambda u: (np.abs(u) <= 1) * 15 * (63 * u ** 4 - 70 * u ** 2 + 15) / 128
            if r == 2:
                return lambda u: (np.abs(u) <= 1) * 105 * (-45 * u ** 4 + 42 * u ** 2 - 5) / 32
        if kernel == "epa":
            if r == 0:
                return _epa_v6_r0
            if r == 2:
                return _epa_v6_r2
    raise ValueError(
        f"kd kernel not implemented for v={v}, r={r}, kernel={kernel!r}."
    )


def kd_K_fun(x: np.ndarray, v: int, r: int, kernel: str) -> Tuple[np.ndarray, float, float]:
    """Evaluate a higher-order kernel and return ``(Kx, k_v, R_v)``.

    Parameters
    ----------
    x     : array-like of kernel arguments.
    v     : kernel order.
    r     : derivative order.
    kernel: "epa", "uni", or "gau" (triangular not supported for v >= 2 in this family).

    Returns
    -------
    Kx   : kernel evaluated at ``x``.
    k_v  : (-1)^v * v-th moment of the kernel, divided by v!.
    R_v  : integrated squared kernel ``int k(u)^2 du``.
    """
    k = _kd_k_impl(v, r, kernel)
    Kx = k(np.asarray(x, dtype=float))
    k_v = integrate.quad(lambda u: ((-1) ** v) * (u ** v) * k(u) / math.factorial(v),
                         -np.inf, np.inf)[0]
    R_v = integrate.quad(lambda u: k(u) ** 2, -np.inf, np.inf)[0]
    return Kx, k_v, R_v


def kd_bw_fun(V: float, B: float, N: int, v: int, r: int) -> float:
    """The closed-form DPI bandwidth constant."""
    return ((1 + 2 * r) * V / (2 * v * N * B ** 2)) ** (1 / (1 + 2 * v + 2 * r))
