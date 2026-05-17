## broom-style tidy / glance outputs have the expected columns, types,
## and CI interpretation (tau.bc +/- z * se.rb).

skip_if_no_broom <- function() {
  testthat::skip_if_not_installed("broom")
}

test_that("tidy.lprobust returns expected columns", {
  skip_if_no_broom()
  fx <- make_fixture(n = 300)
  m <- lprobust(fx$y, fx$x, eval = c(0, 0.5))
  td <- broom::tidy(m)
  expect_s3_class(td, "data.frame")
  expect_named(td, c("eval","estimate","std.error","tau.bc","se.rb",
                     "h","b","n.eff","conf.low","conf.high"))
  z <- qnorm(0.975)
  expect_equal(td$conf.low,  td$tau.bc - z * td$se.rb, tolerance = 1e-12)
  expect_equal(td$conf.high, td$tau.bc + z * td$se.rb, tolerance = 1e-12)
})

test_that("glance.lprobust returns a one-row summary", {
  skip_if_no_broom()
  fx <- make_fixture(n = 300)
  m <- lprobust(fx$y, fx$x, eval = 0.5)
  gl <- broom::glance(m)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("n","neval","p","deriv","kernel","bwselect") %in% names(gl)))
})

test_that("tidy/glance on kdrobust and lpbwselect work", {
  skip_if_no_broom()
  fx <- make_fixture(n = 300)
  kd <- kdrobust(fx$x, eval = 0)
  expect_s3_class(broom::tidy(kd), "data.frame")
  expect_equal(nrow(broom::glance(kd)), 1L)

  bw <- lpbwselect(fx$y, fx$x, eval = 0)
  expect_s3_class(broom::tidy(bw), "data.frame")
  expect_equal(nrow(broom::glance(bw)), 1L)
})
