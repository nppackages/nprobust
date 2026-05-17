## kdrobust should approximate the true density for a large sample and
## integrate close to 1 over the support.

test_that("kdrobust recovers standard normal density at 0", {
  set.seed(1)
  x <- rnorm(30000)
  kd <- kdrobust(x, eval = 0, h = 0.3, kernel = "epa")
  expect_equal(as.numeric(kd$Estimate[, "tau.us"]), dnorm(0), tolerance = 0.01)
  expect_equal(as.numeric(kd$Estimate[, "tau.bc"]), dnorm(0), tolerance = 0.01)
})

test_that("kdrobust integrates to ~1 over a wide grid", {
  set.seed(1)
  x <- rnorm(20000)
  grid_step <- 0.1
  eval <- seq(-4, 4, by = grid_step)
  kd <- kdrobust(x, eval = eval, h = 0.4, kernel = "epa")
  integral <- sum(kd$Estimate[, "tau.us"]) * grid_step
  expect_equal(integral, 1, tolerance = 0.02)
})

test_that("kdrobust point estimate is non-negative for epa/uni", {
  set.seed(1)
  x <- rnorm(500)
  kd <- kdrobust(x, eval = seq(-2, 2, length.out = 20), h = 0.3, kernel = "epa")
  expect_true(all(kd$Estimate[, "tau.us"] >= 0))
})

test_that("kdrobust standard errors are positive and finite", {
  set.seed(1)
  x <- rnorm(500)
  kd <- kdrobust(x, eval = seq(-1, 1, length.out = 10))
  expect_true(all(is.finite(kd$Estimate[, "se.us"])))
  expect_true(all(kd$Estimate[, "se.us"] > 0))
  expect_true(all(is.finite(kd$Estimate[, "se.rb"])))
  expect_true(all(kd$Estimate[, "se.rb"] > 0))
})
