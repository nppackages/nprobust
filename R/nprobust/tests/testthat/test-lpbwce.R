## Regression test pinning lpbwce numerical output.
## Reference values were generated from the original O(N^2) matrix
## implementation; the current closed-form contraction must agree to
## machine precision.

test_that("lpbwce closed-form matches O(N^2) reference", {
  set.seed(42)
  n   <- 400
  x   <- sort(runif(n, -1, 1))
  y   <- 1 + 2 * x + 3 * x^2 + rnorm(n, sd = 0.3)
  res <- y - (1 + 2 * x + 3 * x^2)
  h   <- 0.6
  b   <- 0.9
  K   <- nprobust:::W.fun((x - 0) / h, "epa")
  L   <- nprobust:::W.fun((x - 0) / b, "epa")

  out <- nprobust:::lpbwce(y = y, x = x, K = K, L = L, res = res,
                           c = 0, p = 1, q = 2, h = h, b = b,
                           deriv = 0, fact = 1)

  expect_equal(out$q1rbc, -8.71677124715237, tolerance = 1e-12)
  expect_equal(out$q2rbc, -11.0080523759492, tolerance = 1e-12)
  expect_equal(out$q3rbc,   2.85555479540377, tolerance = 1e-12)
})
