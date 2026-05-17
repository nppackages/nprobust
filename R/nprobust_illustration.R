########################################################################
## NPROBUST Package
## Numerical Illustration
########################################################################
rm(list = ls(all = TRUE))
library(ggplot2)
library(nprobust)

cat("\n=== Setup: cholesterol trial data ===\n")
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)][1])
script_dir <- if (!is.na(script_file)) {
  dirname(normalizePath(script_file, mustWork = FALSE))
} else {
  getwd()
}
chole <- read.csv(file.path(script_dir, "nprobust_data.csv"))
summary(chole)

t <- chole$t
chol1 <- chole$chol1
chol2 <- chole$chol2
cholf <- chole$cholf
comp <- chole$comp
control <- t == 0
treated <- t == 1

cat("\n=== 1. Kernel density reports for baseline cholesterol ===\n")
# Default IMSE-DPI bandwidth on seven evaluation points.
summary(kdrobust(chol1, subset = control, neval = 7))

# Pointwise MSE-DPI bandwidth on seven evaluation points.
summary(kdrobust(chol1, subset = control, neval = 7, bwselect = "MSE-DPI"))

# IMSE-DPI bandwidth on the default denser 30-point grid.
summary(kdrobust(chol1, subset = control, neval = 30, bwselect = "IMSE-DPI"))

cat("\n=== 2. Kernel density estimates on common grids ===\n")
# Common grids for cholesterol variables and treatment compliance.
grid_chol <- seq(250, 350, length.out = 21)
grid_comp <- seq(0, 100, length.out = 21)

# Kernel density estimates for baseline, follow-up, and compliance variables.
f0_chol1 <- kdrobust(chol1, subset = control, eval = grid_chol)
f1_chol1 <- kdrobust(chol1, subset = treated, eval = grid_chol)
f0_chol2 <- kdrobust(chol2, subset = control, eval = grid_chol)
f1_chol2 <- kdrobust(chol2, subset = treated, eval = grid_chol)
f0_cholf <- kdrobust(cholf, subset = control, eval = grid_chol)
f1_cholf <- kdrobust(cholf, subset = treated, eval = grid_chol)
f0_comp <- kdrobust(comp, subset = control, eval = grid_comp)
f1_comp <- kdrobust(comp, subset = treated, eval = grid_comp)

# Report the first five common-grid density estimates for baseline cholesterol.
head(f0_chol1$Estimate, 5)
head(f1_chol1$Estimate, 5)

# Kernel density plots for the main trial variables.
p_kd_chol1 <- nprobust.plot(f0_chol1, f1_chol1, legendGroups = c("Control Group", "Treatment Group"), xlabel = "Cholesterol at Baseline 1", ylabel = "Density") + theme(legend.position = c(.4, .2))
p_kd_chol2 <- nprobust.plot(f0_chol2, f1_chol2, xlabel = "Cholesterol at Baseline 2", ylabel = "Density") + theme(legend.position = "none")
p_kd_cholf <- nprobust.plot(f0_cholf, f1_cholf, xlabel = "Cholesterol after Treatment", ylabel = "Density") + theme(legend.position = "none")
p_kd_comp <- nprobust.plot(f0_comp, f1_comp, xlabel = "Treatment Compliance", ylabel = "Density") + theme(legend.position = "none")

cat("\n=== 3. Difference in means ===\n")
# Difference in means for outcome and compliance.
t.test(cholf[control], cholf[treated])
t.test(comp[control], comp[treated])

cat("\n=== 4. Local polynomial regression on common grids ===\n")
# Local polynomial regression estimates for outcomes and compliance.
m0_cholf_1 <- lprobust(cholf, chol1, subset = control, eval = grid_chol)
m1_cholf_1 <- lprobust(cholf, chol1, subset = treated, eval = grid_chol)
m0_comp_1 <- lprobust(comp, chol1, subset = control, eval = grid_chol)
m1_comp_1 <- lprobust(comp, chol1, subset = treated, eval = grid_chol)
m0_cholf_2 <- lprobust(cholf, chol2, subset = control, eval = grid_chol)
m1_cholf_2 <- lprobust(cholf, chol2, subset = treated, eval = grid_chol)
m0_comp_2 <- lprobust(comp, chol2, subset = control, eval = grid_chol)
m1_comp_2 <- lprobust(comp, chol2, subset = treated, eval = grid_chol)

# Report local polynomial regression for the first seven grid points.
grid7 <- grid_chol[1:7]
summary(lprobust(cholf, chol1, subset = control, eval = grid7))

# Local polynomial regression plots for the main trial variables.
p_lp_cholf_1 <- nprobust.plot(m0_cholf_1, m1_cholf_1, legendGroups = c("Control Group", "Treatment Group"), xlabel = "Cholesterol at Baseline 1", ylabel = "Cholesterol after Treatment") + theme(legend.position = c(.3, .8))
p_lp_cholf_2 <- nprobust.plot(m0_cholf_2, m1_cholf_2, xlabel = "Cholesterol at Baseline 2", ylabel = "Cholesterol after Treatment") + theme(legend.position = "none")
p_lp_comp_1 <- nprobust.plot(m0_comp_1, m1_comp_1, xlabel = "Cholesterol at Baseline 1", ylabel = "Treatment Compliance") + theme(legend.position = "none")
p_lp_comp_2 <- nprobust.plot(m0_comp_2, m1_comp_2, xlabel = "Cholesterol at Baseline 2", ylabel = "Treatment Compliance") + theme(legend.position = "none")

cat("\n=== 5. Bandwidth selection ===\n")
# Local polynomial bandwidth selection with MSE-DPI, CE-DPI, and all selectors.
summary(lpbwselect(cholf, chol1, subset = control, eval = grid7))
summary(lpbwselect(cholf, chol1, subset = control, eval = grid7, bwselect = "CE-DPI"))
summary(lpbwselect(cholf, chol1, subset = control, eval = grid7, bwselect = "ALL"))

# Kernel density bandwidth selection with MSE-DPI and IMSE-DPI selectors.
summary(kdbwselect(chol1, subset = control))
summary(kdbwselect(chol1, subset = control, bwselect = "IMSE-DPI"))

cat("\n=== 6. Efron-Feldman compliance plot ===\n")
# Replication of the Efron and Feldman (1991) compliance plot.
y_ef <- 0.25 * chol1 + 0.75 * chol2 - cholf
x_ef <- comp
m0_ef <- lprobust(y_ef, x_ef, subset = control, neval = 100)
m1_ef <- lprobust(y_ef, x_ef, subset = treated, neval = 100)
p_ef <- nprobust.plot(m1_ef, m0_ef, legendGroups = c("Treatment Group", "Control Group"), xlabel = "Compliance", ylabel = "Cholesterol Difference")

# Print plots only in interactive sessions to avoid creating Rplots.pdf in batch runs.
if (interactive()) {
  print(p_kd_chol1); print(p_kd_chol2); print(p_kd_cholf); print(p_kd_comp)
  print(p_lp_cholf_1); print(p_lp_cholf_2); print(p_lp_comp_1); print(p_lp_comp_2)
  print(p_ef)
}
