## kdrobust point estimate must equal the direct kernel density estimator
## mean(K((x-eval)/h))/h for kernels that have closed-form support.

test_that("kdrobust point estimate matches direct kernel density", {
  set.seed(3)
  xx <- rnorm(2000)
  for (kernel in c("epa", "uni")) {
    for (eval_pt in c(-0.5, 0, 0.7)) {
      h <- 0.35
      kd <- kdrobust(xx, eval = eval_pt, h = h, kernel = kernel)
      tau_ref <- mean(w_fun((xx - eval_pt) / h, kernel)) / h
      expect_equal(as.numeric(kd$Estimate[, "tau.us"]),
                   as.numeric(tau_ref),
                   tolerance = 1e-12,
                   label = sprintf("kernel=%s eval=%+.1f", kernel, eval_pt))
    }
  }
})

test_that("kdrobust rejects unsupported kernels with a clear error", {
  set.seed(3)
  xx <- rnorm(200)
  expect_error(kdrobust(xx, kernel = "tri"), "kernel incorrectly specified")
  expect_error(kdrobust(xx, kernel = "gau"), "kernel incorrectly specified")
})
