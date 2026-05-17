################################################################################
# NPROBUST Package
# Numerical Illustration
################################################################################
from __future__ import annotations

import os
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

from nprobust import (
    kdrobust,
    kdbwselect,
    lpbwselect,
    lprobust,
    plot_kdrobust,
    plot_lprobust,
)


def section(title: str) -> None:
    print(f"\n=== {title} ===")


def close_plot(fig) -> None:
    try:
        import matplotlib.pyplot as plt

        plt.close(fig)
    except Exception:
        pass


section("Setup: cholesterol trial data")
here = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="nprobust-mpl-"))
chole = pd.read_csv(os.path.join(here, "nprobust_data.csv"))
print(chole.describe())

t = chole["t"].to_numpy()
chol1 = chole["chol1"].to_numpy()
chol2 = chole["chol2"].to_numpy()
cholf = chole["cholf"].to_numpy()
comp = chole["comp"].to_numpy()
control = t == 0
treated = t == 1

section("1. Kernel density reports for baseline cholesterol")
# Default IMSE-DPI bandwidth on seven evaluation points.
print(kdrobust(chol1[control], neval=7).summary())

# Pointwise MSE-DPI bandwidth on seven evaluation points.
print(kdrobust(chol1[control], neval=7, bwselect="mse-dpi").summary())

# IMSE-DPI bandwidth on the default denser 30-point grid.
print(kdrobust(chol1[control], neval=30, bwselect="imse-dpi").summary())

section("2. Kernel density estimates on common grids")
# Common grids for cholesterol variables and treatment compliance.
grid_chol = np.linspace(250, 350, 21)
grid_comp = np.linspace(0, 100, 21)

# Kernel density estimates for baseline, follow-up, and compliance variables.
f0_chol1 = kdrobust(chol1[control], eval=grid_chol)
f1_chol1 = kdrobust(chol1[treated], eval=grid_chol)
f0_chol2 = kdrobust(chol2[control], eval=grid_chol)
f1_chol2 = kdrobust(chol2[treated], eval=grid_chol)
f0_cholf = kdrobust(cholf[control], eval=grid_chol)
f1_cholf = kdrobust(cholf[treated], eval=grid_chol)
f0_comp = kdrobust(comp[control], eval=grid_comp)
f1_comp = kdrobust(comp[treated], eval=grid_comp)

# Report the first five common-grid density estimates for baseline cholesterol.
print(f0_chol1.Estimate.head())
print(f1_chol1.Estimate.head())

# Kernel density plots for the main trial variables.
try:
    kd_figures = [
        plot_kdrobust(f0_chol1, f1_chol1, labels=["Control Group", "Treatment Group"],
                      xlabel="Cholesterol at Baseline 1", ylabel="Density"),
        plot_kdrobust(f0_chol2, f1_chol2, xlabel="Cholesterol at Baseline 2", ylabel="Density"),
        plot_kdrobust(f0_cholf, f1_cholf, xlabel="Cholesterol after Treatment", ylabel="Density"),
        plot_kdrobust(f0_comp, f1_comp, xlabel="Treatment Compliance", ylabel="Density"),
    ]
    for fig in kd_figures:
        close_plot(fig)
except ImportError as exc:
    print(f"Plotting skipped: {exc}")

section("3. Difference in means")
# Difference in means for outcome and compliance.
print(stats.ttest_ind(cholf[control], cholf[treated], equal_var=False))
print(stats.ttest_ind(comp[control], comp[treated], equal_var=False))

section("4. Local polynomial regression on common grids")
# Local polynomial regression estimates for outcomes and compliance.
m0_cholf_1 = lprobust(cholf[control], chol1[control], eval=grid_chol)
m1_cholf_1 = lprobust(cholf[treated], chol1[treated], eval=grid_chol)
m0_comp_1 = lprobust(comp[control], chol1[control], eval=grid_chol)
m1_comp_1 = lprobust(comp[treated], chol1[treated], eval=grid_chol)
m0_cholf_2 = lprobust(cholf[control], chol2[control], eval=grid_chol)
m1_cholf_2 = lprobust(cholf[treated], chol2[treated], eval=grid_chol)
m0_comp_2 = lprobust(comp[control], chol2[control], eval=grid_chol)
m1_comp_2 = lprobust(comp[treated], chol2[treated], eval=grid_chol)

# Report local polynomial regression for the first seven grid points.
grid7 = grid_chol[:7]
print(lprobust(cholf[control], chol1[control], eval=grid7).summary())

# Local polynomial regression plots for the main trial variables.
try:
    lp_figures = [
        plot_lprobust(m0_cholf_1, m1_cholf_1, labels=["Control Group", "Treatment Group"],
                      xlabel="Cholesterol at Baseline 1", ylabel="Cholesterol after Treatment"),
        plot_lprobust(m0_cholf_2, m1_cholf_2, xlabel="Cholesterol at Baseline 2",
                      ylabel="Cholesterol after Treatment"),
        plot_lprobust(m0_comp_1, m1_comp_1, xlabel="Cholesterol at Baseline 1",
                      ylabel="Treatment Compliance"),
        plot_lprobust(m0_comp_2, m1_comp_2, xlabel="Cholesterol at Baseline 2",
                      ylabel="Treatment Compliance"),
    ]
    for fig in lp_figures:
        close_plot(fig)
except ImportError as exc:
    print(f"Plotting skipped: {exc}")

section("5. Bandwidth selection")
# Local polynomial bandwidth selection with MSE-DPI, CE-DPI, and all selectors.
print(lpbwselect(cholf[control], chol1[control], eval=grid7).bws)
print(lpbwselect(cholf[control], chol1[control], eval=grid7, bwselect="ce-dpi").bws)
print(lpbwselect(cholf[control], chol1[control], eval=grid7, bwselect="all").bws)

# Kernel density bandwidth selection with MSE-DPI and IMSE-DPI selectors.
print(kdbwselect(chol1[control]).bws)
print(kdbwselect(chol1[control], bwselect="imse-dpi").bws)

section("6. Efron-Feldman compliance plot")
# Replication of the Efron and Feldman (1991) compliance plot.
y_ef = 0.25 * chol1 + 0.75 * chol2 - cholf
x_ef = comp
m0_ef = lprobust(y_ef[control], x_ef[control], neval=100)
m1_ef = lprobust(y_ef[treated], x_ef[treated], neval=100)

try:
    close_plot(
        plot_lprobust(
            m1_ef,
            m0_ef,
            labels=["Treatment Group", "Control Group"],
            xlabel="Compliance",
            ylabel="Cholesterol Difference",
        )
    )
except ImportError as exc:
    print(f"Plotting skipped: {exc}")

print("\nDone.")
