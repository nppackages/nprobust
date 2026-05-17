## Point estimates from lprobust must match the intercept of the equivalent
## weighted least squares fit (uniform kernel = OLS on bandwidth subset).

test_that("lprobust point estimates match lm / WLS across kernels, p, eval", {
  fx <- make_fixture()

  for (kernel in c("uni", "tri", "epa")) {
    for (p in 1:3) {
      for (eval in c(-0.3, 0, 0.4)) {
        h <- 0.25
        est <- lprobust(fx$y, fx$x, eval = eval, h = h, p = p,
                        kernel = kernel, vce = "hc0")

        u  <- (fx$x - eval) / h
        w  <- w_fun(u, kernel)
        sub <- w > 0
        xc <- fx$x[sub] - eval
        ys <- fx$y[sub]
        ws <- w[sub]
        X  <- cbind(1, poly(xc, degree = p, raw = TRUE))
        beta <- solve(crossprod(X * sqrt(ws)), crossprod(X * ws, ys))

        expect_equal(as.numeric(est$Estimate[, "tau.us"]),
                     as.numeric(beta[1]),
                     tolerance = 1e-10,
                     label = sprintf("kernel=%s p=%d eval=%+.1f", kernel, p, eval))
      }
    }
  }
})
