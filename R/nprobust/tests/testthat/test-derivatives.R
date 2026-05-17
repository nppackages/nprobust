## For deriv=k, lprobust should return k! * beta_k from the local polynomial.

test_that("lprobust derivative estimates match deriv! * local poly coefficient", {
  fx <- make_fixture()
  kernel <- "uni"
  h <- 0.3
  eval <- 0

  for (deriv in 0:2) {
    for (p in (deriv + 1):(deriv + 2)) {
      est <- lprobust(fx$y, fx$x, eval = eval, h = h, p = p, deriv = deriv,
                      kernel = kernel, vce = "hc0")

      sub <- abs((fx$x - eval) / h) <= 1
      xc <- fx$x[sub] - eval
      ys <- fx$y[sub]
      X <- cbind(1, poly(xc, degree = p, raw = TRUE))
      beta <- solve(crossprod(X), crossprod(X, ys))
      tau_ref <- factorial(deriv) * beta[deriv + 1]

      expect_equal(as.numeric(est$Estimate[, "tau.us"]),
                   as.numeric(tau_ref),
                   tolerance = 1e-10,
                   label = sprintf("deriv=%d p=%d", deriv, p))
    }
  }
})
