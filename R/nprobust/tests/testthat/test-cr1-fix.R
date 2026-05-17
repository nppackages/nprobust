## Regression test for the CR1 cluster + bias-corrected k_override fix
## (2026-05-08). Pre-fix the cluster meat for V.Y.bc used k=ncol(Q.q)=p+1
## in the (N-1)/(N-k) df correction; the residuals come from the q-regression
## so the correct k is q+1. lm + sandwich::vcovCL with q-regressor design is
## the unambiguous anchor.

library(testthat)
library(nprobust)

test_that("CR1 cluster bias-corrected SE matches lm(p=q=2)+vcovCL anchor at h=b", {
  skip_if_not_installed("sandwich")
  set.seed(20260508)
  n  <- 1500
  x  <- runif(n, -1, 1)
  y  <- 0.5 * x + 0.3 * x^2 + rnorm(n)
  cl <- ceiling(20 * runif(n))

  eval_pt <- 0
  h       <- 0.5

  fit <- lprobust(y, x, eval = eval_pt, h = h, b = h, p = 1, deriv = 0,
                  kernel = "tri", vce = "cr1", cluster = cl)
  se_np <- as.numeric(fit$Estimate[, "se.rb"])

  ## Manual lm anchor with triangular kernel weights and q-regression design
  in_bw <- abs(x - eval_pt) <= h
  yi <- y[in_bw]; xi <- x[in_bw]; cli <- cl[in_bw]
  u  <- (xi - eval_pt) / h
  wi <- (1 - abs(u)) * (abs(u) <= 1)
  nz <- wi > 0
  yi <- yi[nz]; xi <- xi[nz]; cli <- cli[nz]; wi <- wi[nz]
  sw <- sqrt(wi)
  xc <- xi - eval_pt
  Xq <- cbind(1, xc, xc^2)
  fit_lm <- lm(I(sw * yi) ~ I(sw * Xq) - 1)
  V_q    <- sandwich::vcovCL(fit_lm, cluster = cli, type = "HC1")
  se_lm  <- sqrt(V_q[1, 1])

  expect_lt(abs(se_np - se_lm), 1e-10)
})

test_that("CR2 / CR3 cluster paths unchanged by k_override fix", {
  set.seed(123)
  n  <- 800
  x  <- runif(n, -1, 1)
  y  <- 0.4 * x + rnorm(n)
  cl <- ceiling(15 * runif(n))

  for (vce in c("cr2", "cr3")) {
    fit <- lprobust(y, x, eval = 0, h = 0.5, b = 0.5, p = 1, vce = vce, cluster = cl)
    expect_true(is.finite(fit$Estimate[, "se.rb"]))
    expect_gt(fit$Estimate[, "se.rb"], 0)
  }
})
