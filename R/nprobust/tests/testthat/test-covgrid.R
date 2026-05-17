## covgrid=TRUE should return symmetric covariance matrices whose diagonal
## equals se.us^2 / se.rb^2 at each eval point.

test_that("covgrid diagonal matches squared SEs", {
  fx <- make_fixture(n = 400)
  est <- lprobust(fx$y, fx$x, eval = seq(0.2, 0.8, length.out = 5),
                  covgrid = TRUE)

  expect_equal(dim(est$cov.us), c(5, 5))
  expect_equal(dim(est$cov.rb), c(5, 5))
  expect_true(isSymmetric(est$cov.us))
  expect_true(isSymmetric(est$cov.rb))

  expect_equal(diag(est$cov.us), as.numeric(est$Estimate[, "se.us"]^2),
               tolerance = 1e-10)
  expect_equal(diag(est$cov.rb), as.numeric(est$Estimate[, "se.rb"]^2),
               tolerance = 1e-10)
})

test_that("covgrid scaling is correct for deriv >= 2", {
  ## Previously V.us.i had factorial(deriv)^2 baked in, so
  ## (V.us.i %*% t(V.us.j))[d+1,d+1] picked up factorial(deriv)^4 instead
  ## of factorial(deriv)^2 -- inflating cov entries by factorial(deriv)^2
  ## for deriv >= 2 (factor 4 at deriv=2; 36 at deriv=3).
  fx <- make_fixture(n = 600)
  for (d in c(0, 1, 2, 3)) {
    est <- lprobust(fx$y, fx$x, p = d + 1, deriv = d,
                    eval = c(0.2, 0.5, 0.7), covgrid = TRUE)
    expect_equal(diag(est$cov.us), as.numeric(est$Estimate[, "se.us"]^2),
                 tolerance = 1e-10,
                 label = sprintf("deriv=%d cov.us diag", d))
    expect_equal(diag(est$cov.rb), as.numeric(est$Estimate[, "se.rb"]^2),
                 tolerance = 1e-10,
                 label = sprintf("deriv=%d cov.rb diag", d))
  }
})
