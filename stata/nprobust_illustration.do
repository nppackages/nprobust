********************************************************************************
** NPROBUST Package
** Numerical Illustration
********************************************************************************
clear all
set more off
set scheme sj

di _newline "=== Setup: cholesterol trial data ==="
findfile nprobust_data.dta
use "`r(fn)'", clear
summarize

********************************************************************************
** Common grids for cholesterol variables and treatment compliance.
********************************************************************************
gen grid_chol = 245 + 5 * _n if _n <= 21
gen grid_comp = -5 + 5 * _n if _n <= 21
gen grid7 = grid_chol if _n <= 7

di _newline "=== 1. Kernel density reports for baseline cholesterol ==="
** Default IMSE-DPI bandwidth on seven evaluation points.
kdrobust chol1 if t==0, neval(7)

** Pointwise MSE-DPI bandwidth on seven evaluation points.
kdrobust chol1 if t==0, neval(7) bwselect(mse-dpi)

** IMSE-DPI bandwidth on the default denser 30-point grid.
kdrobust chol1 if t==0, neval(30) bwselect(imse-dpi)

di _newline "=== 2. Kernel density estimates on common grids ==="
** Kernel density estimates and plots for baseline cholesterol.
kdrobust chol1 if t==0, eval(grid_chol) plot graph_options(name(kd_chol1_t0, replace)) bwselect(imse-dpi)
kdrobust chol1 if t==1, eval(grid_chol) plot graph_options(name(kd_chol1_t1, replace)) bwselect(imse-dpi)
graph combine kd_chol1_t0 kd_chol1_t1, name(kd_chol1, replace)

** Kernel density estimates and plots for second baseline cholesterol.
kdrobust chol2 if t==0, eval(grid_chol) plot graph_options(name(kd_chol2_t0, replace)) bwselect(imse-dpi)
kdrobust chol2 if t==1, eval(grid_chol) plot graph_options(name(kd_chol2_t1, replace)) bwselect(imse-dpi)
graph combine kd_chol2_t0 kd_chol2_t1, name(kd_chol2, replace)

** Kernel density estimates and plots for cholesterol after treatment.
kdrobust cholf if t==0, eval(grid_chol) plot graph_options(name(kd_cholf_t0, replace)) bwselect(imse-dpi)
kdrobust cholf if t==1, eval(grid_chol) plot graph_options(name(kd_cholf_t1, replace)) bwselect(imse-dpi)
graph combine kd_cholf_t0 kd_cholf_t1, name(kd_cholf, replace)

** Kernel density estimates and plots for treatment compliance.
kdrobust comp if t==0, eval(grid_comp) plot graph_options(name(kd_comp_t0, replace)) bwselect(imse-dpi)
kdrobust comp if t==1, eval(grid_comp) plot graph_options(name(kd_comp_t1, replace)) bwselect(imse-dpi)
graph combine kd_comp_t0 kd_comp_t1, name(kd_comp, replace)

di _newline "=== 3. Difference in means ==="
** Difference in means for outcome and compliance.
ttest cholf, by(t) unequal
ttest comp, by(t) unequal

di _newline "=== 4. Local polynomial regression on common grids ==="
** Local polynomial regression estimates and plots for cholf on chol1.
lprobust cholf chol1 if t==0, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_cholf1_t0, replace))
lprobust cholf chol1 if t==1, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_cholf1_t1, replace))
graph combine lp_cholf1_t0 lp_cholf1_t1, name(lp_cholf1, replace)

** Local polynomial regression estimates and plots for comp on chol1.
lprobust comp chol1 if t==0, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_comp1_t0, replace))
lprobust comp chol1 if t==1, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_comp1_t1, replace))
graph combine lp_comp1_t0 lp_comp1_t1, name(lp_comp1, replace)

** Local polynomial regression estimates and plots for cholf on chol2.
lprobust cholf chol2 if t==0, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_cholf2_t0, replace))
lprobust cholf chol2 if t==1, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_cholf2_t1, replace))
graph combine lp_cholf2_t0 lp_cholf2_t1, name(lp_cholf2, replace)

** Local polynomial regression estimates and plots for comp on chol2.
lprobust comp chol2 if t==0, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_comp2_t0, replace))
lprobust comp chol2 if t==1, eval(grid_chol) bwcheck(21) plot graph_options(name(lp_comp2_t1, replace))
graph combine lp_comp2_t0 lp_comp2_t1, name(lp_comp2, replace)

** Report local polynomial regression for the first seven grid points.
lprobust cholf chol1 if t==0, eval(grid7)

di _newline "=== 5. Bandwidth selection ==="
** Local polynomial bandwidth selection with MSE-DPI, CE-DPI, and all selectors.
lpbwselect cholf chol1 if t==0, eval(grid7)
lpbwselect cholf chol1 if t==0, eval(grid7) bwselect(ce-dpi)
lpbwselect cholf chol1 if t==0, eval(grid7) bwselect(all)

** Kernel density bandwidth selection with MSE-DPI and IMSE-DPI selectors.
kdbwselect chol1 if t==0
kdbwselect chol1 if t==0, bwselect(imse-dpi)

di _newline "=== 6. Efron-Feldman compliance plot ==="
** Replication of the Efron and Feldman (1991) compliance plot.
gen y_ef = 0.25 * chol1 + 0.75 * chol2 - cholf
lprobust y_ef comp if t==0, neval(100) plot graph_options(name(ef_t0, replace))
lprobust y_ef comp if t==1, neval(100) plot graph_options(name(ef_t1, replace))
graph combine ef_t1 ef_t0, name(ef_compliance, replace)
