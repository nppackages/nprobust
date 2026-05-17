"""Plotting helpers for nprobust results (matplotlib)."""
from __future__ import annotations

from typing import Optional, Sequence

from scipy.stats import norm


def _import_plt():  # pragma: no cover - trivial
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for nprobust plotting. "
            "Install with: pip install nprobust[plot]"
        ) from e
    return plt


def _subplots_with_fallback(plt):
    try:
        return plt.subplots()
    except Exception as first_error:  # pragma: no cover - backend dependent
        try:
            import matplotlib

            matplotlib.use("Agg", force=True)
            import matplotlib.pyplot as fallback_plt

            return fallback_plt.subplots()
        except Exception:
            raise first_error


def _plot_series(
    estimates,
    *,
    alpha: float = 0.05,
    labels: Optional[Sequence[str]] = None,
    ci_type: str = "region",
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    ax=None,
):
    plt = _import_plt()
    if not isinstance(estimates, (list, tuple)):
        estimates = [estimates]
    if labels is None:
        labels = [f"Series {i + 1}" for i in range(len(estimates))]

    if ax is None:
        fig, ax = _subplots_with_fallback(plt)
    else:
        fig = ax.figure

    z = norm.ppf(1 - alpha / 2)

    for est, lbl in zip(estimates, labels):
        df = est.Estimate
        ci_l = df["tau.bc"] - z * df["se.rb"]
        ci_u = df["tau.bc"] + z * df["se.rb"]
        line, = ax.plot(df["eval"], df["tau.us"], label=lbl)
        color = line.get_color()
        if ci_type in ("region", "all"):
            ax.fill_between(df["eval"], ci_l, ci_u, alpha=0.2, color=color)
        if ci_type in ("line", "all"):
            ax.plot(df["eval"], ci_l, linestyle="--", color=color, alpha=0.5)
            ax.plot(df["eval"], ci_u, linestyle="--", color=color, alpha=0.5)

    ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if len(estimates) > 1 or (labels and labels[0] != "Series 1"):
        ax.legend()
    return fig


def plot_lprobust(
    *estimates,
    alpha: float = 0.05,
    labels: Optional[Sequence[str]] = None,
    ci_type: str = "region",
    title: str = "",
    xlabel: str = "x",
    ylabel: str = "m(x)",
    ax=None,
):
    """Plot one or more :class:`LprobustResult` objects.

    Parameters
    ----------
    *estimates : LprobustResult
        One or more fitted local polynomial results.
    alpha : float
        Significance level for the confidence band (default 0.05).
    labels : sequence of str, optional
        Legend labels, one per series.
    ci_type : {"region", "line", "all", "none"}
        How to draw the robust bias-corrected confidence interval.
    title, xlabel, ylabel : str
        Plot annotations.
    ax : matplotlib Axes, optional
        Axes to draw on. A new figure/axes is created if not provided.

    Returns
    -------
    matplotlib.figure.Figure
    """
    return _plot_series(list(estimates), alpha=alpha, labels=labels,
                        ci_type=ci_type, title=title, xlabel=xlabel,
                        ylabel=ylabel, ax=ax)


def plot_kdrobust(
    *estimates,
    alpha: float = 0.05,
    labels: Optional[Sequence[str]] = None,
    ci_type: str = "region",
    title: str = "",
    xlabel: str = "x",
    ylabel: str = "f(x)",
    ax=None,
):
    """Plot one or more :class:`KdrobustResult` objects."""
    return _plot_series(list(estimates), alpha=alpha, labels=labels,
                        ci_type=ci_type, title=title, xlabel=xlabel,
                        ylabel=ylabel, ax=ax)
