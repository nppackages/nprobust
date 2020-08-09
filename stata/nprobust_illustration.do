******************************************************************************** 
** NPROBUST Stata Package
** Do-file for Empirical Illustration 
** Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell 
********************************************************************************
** hlp2winpdf, cdn(lprobust) replace
** hlp2winpdf, cdn(kdrobust) replace
** hlp2winpdf, cdn(lpbwselect) replace
** hlp2winpdf, cdn(kdbwselect) replace
********************************************************************************
** net install nprobust, from(https://sites.google.com/site/nppackages/nprobust/stata) replace
********************************************************************************
clear all
set more off
set scheme sj

use nprobust_data.dta, clear

********************************************************************************
** generate grid points from 250 to 350 by 5
********************************************************************************
gen grid = 245 + 5 * _n if _n <= 21

********************************************************************************
** Kernel Density Estimation with Robust Bias-Corrected Confidence Intervals using MSE-DPI optimal bandwidth.
********************************************************************************
kdrobust chol1 if t==0, eval(grid) plot graph_options(name(chol1_t0)) bwselect(imse-dpi)
kdrobust chol1 if t==1, eval(grid) plot graph_options(name(chol1_t1)) bwselect(imse-dpi)
graph combine chol1_t0 chol1_t1

kdrobust chol2 if t==0, eval(grid) plot graph_options(name(chol2_t0)) bwselect(imse-dpi)
kdrobust chol2 if t==1, eval(grid) plot graph_options(name(chol2_t1)) bwselect(imse-dpi)
graph combine chol2_t0 chol2_t1

kdrobust cholf if t==0, eval(grid) plot graph_options(name(cholf_t0)) bwselect(imse-dpi)
kdrobust cholf if t==1, eval(grid) plot graph_options(name(cholf_t1)) bwselect(imse-dpi)
graph combine cholf_t0 cholf_t1

kdrobust comp if t==0, plot graph_options(name(comp_t0)) bwselect(imse-dpi)
kdrobust comp if t==1, plot graph_options(name(comp_t1)) bwselect(imse-dpi)
graph combine comp_t0 comp_t1

********************************************************************************
** Local Polynomial Regression with Robust Bias-Corrected Confidence Intervals using MSE-DPI optimal bandwidth.
********************************************************************************
lprobust cholf chol1 if t==0, eval(grid) bwcheck(21) plot graph_options(name(cholf_t0))
lprobust cholf chol1 if t==1, eval(grid) bwcheck(21) plot graph_options(name(cholf_t1))
graph combine cholf_t0 cholf_t1

lprobust comp  chol1 if t==0, eval(grid) bwcheck(21) plot graph_options(name(comp_t0))
lprobust comp  chol1 if t==1, eval(grid) bwcheck(21) plot graph_options(name(comp_t1))
graph combine comp_t0 comp_t1

lprobust cholf chol2 if t==0, eval(grid) bwcheck(21) plot graph_options(name(cholf_t0))
lprobust cholf chol2 if t==1, eval(grid) bwcheck(21) plot graph_options(name(cholf_t1))
graph combine cholf_t0 cholf_t1

lprobust comp  chol2 if t==0, eval(grid) bwcheck(21) plot graph_options(name(comp_t0))
lprobust comp  chol2 if t==1, eval(grid) bwcheck(21) plot graph_options(name(comp_t1))
graph combine comp_t0 comp_t1

********************************************************************************
** Optimal Bandwidth selectors for Local Polynomial Regression 
********************************************************************************
lpbwselect cholf  chol1 if t==0, eval(grid)
lpbwselect cholf  chol1 if t==0, bwselect(ce-dpi) eval(grid)
lpbwselect cholf  chol1 if t==0, bwselect(all) eval(grid)


