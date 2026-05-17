## User weights should multiply kernel weights throughout: point estimates
## must match a weighted least squares fit with w = K(u) * user_weights.

test_that("weights arg matches WLS with kernel * user weights", {
  fx <- make_fixture(n = 600)
  set.seed(3)
  w <- runif(fx$n, 0.2, 2)

  h <- 0.3; eval <- 0; p <- 1
  est <- lprobust(fx$y, fx$x, eval = eval, h = h, p = p,
                  kernel = "uni", vce = "hc0", weights = w)

  u <- (fx$x - eval) / h
  kw <- w_fun(u, "uni")
  comb <- kw * w
  sub  <- comb > 0
  xc <- fx$x[sub] - eval
  ys <- fx$y[sub]
  ws <- comb[sub]
  X <- cbind(1, xc)
  beta <- solve(crossprod(X * sqrt(ws)), crossprod(X * ws, ys))

  expect_equal(as.numeric(est$Estimate[, "tau.us"]),
               as.numeric(beta[1]),
               tolerance = 1e-10)
})

test_that("weights=1 matches no weights", {
  fx <- make_fixture(n = 400)
  est0 <- lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "epa", vce = "hc0")
  est1 <- lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "epa", vce = "hc0",
                   weights = rep(1, fx$n))
  expect_equal(as.numeric(est0$Estimate), as.numeric(est1$Estimate),
               tolerance = 1e-12)
})

test_that("zero weights exclude observations", {
  fx <- make_fixture(n = 400)
  w <- rep(1, fx$n); w[fx$x > 0.5] <- 0
  est_w <- lprobust(fx$y, fx$x, eval = 0.3, h = 0.3, p = 1, kernel = "uni",
                    vce = "hc0", weights = w)
  est_s <- lprobust(fx$y, fx$x, eval = 0.3, h = 0.3, p = 1, kernel = "uni",
                    vce = "hc0", subset = fx$x <= 0.5)
  expect_equal(as.numeric(est_w$Estimate[, "tau.us"]),
               as.numeric(est_s$Estimate[, "tau.us"]),
               tolerance = 1e-12)
})
