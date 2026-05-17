## print/summary contracts:
##   - print(x) returns x invisibly
##   - summary(x) returns an S3 object (class "summary.<class>")
##   - print(summary(x)) returns its argument invisibly

test_that("print methods return their input invisibly", {
  fx <- make_fixture(n = 300)
  m  <- lprobust(fx$y, fx$x, eval = c(0, 0.5))
  bw <- lpbwselect(fx$y, fx$x, eval = c(0, 0.5))
  kd <- kdrobust(fx$x, eval = c(0, 0.5))
  kb <- kdbwselect(fx$x, eval = c(0, 0.5))

  silence <- file(tempfile(), open = "wt")
  sink(silence)
  on.exit({ sink(); close(silence) }, add = TRUE)

  for (obj in list(m, bw, kd, kb)) {
    res <- withVisible(print(obj))
    expect_false(res$visible)
    expect_identical(res$value, obj)
  }
})

test_that("summary methods return an S3 summary object", {
  fx <- make_fixture(n = 300)
  m  <- lprobust(fx$y, fx$x, eval = c(0, 0.5))
  bw <- lpbwselect(fx$y, fx$x, eval = c(0, 0.5))
  kd <- kdrobust(fx$x, eval = c(0, 0.5))
  kb <- kdbwselect(fx$x, eval = c(0, 0.5))

  capture.output({
    s.m  <- summary(m)
    s.bw <- summary(bw)
    s.kd <- summary(kd)
    s.kb <- summary(kb)
  })

  expect_s3_class(s.m,  "summary.lprobust")
  expect_s3_class(s.bw, "summary.lpbwselect")
  expect_s3_class(s.kd, "summary.kdrobust")
  expect_s3_class(s.kb, "summary.kdbwselect")

  ## CI bounds inside the summary match tau.bc +/- z*se.rb
  z <- qnorm(1 - 0.05 / 2)
  expect_equal(s.m$CI_l,  s.m$Estimate[, "tau.bc"]  - z * s.m$Estimate[, "se.rb"])
  expect_equal(s.m$CI_r,  s.m$Estimate[, "tau.bc"]  + z * s.m$Estimate[, "se.rb"])
  expect_equal(s.kd$CI_l, s.kd$Estimate[, "tau.bc"] - z * s.kd$Estimate[, "se.rb"])
})

test_that("summary respects alpha argument", {
  fx <- make_fixture(n = 300)
  m  <- lprobust(fx$y, fx$x, eval = c(0, 0.5))

  capture.output(s90 <- summary(m, alpha = 0.10))
  z90 <- qnorm(1 - 0.10 / 2)
  expect_equal(s90$CI_l, s90$Estimate[, "tau.bc"] - z90 * s90$Estimate[, "se.rb"])
  expect_equal(s90$alpha, 0.10)
})

test_that("cov.us / cov.rb are NULL when covgrid = FALSE", {
  fx <- make_fixture(n = 300)
  m  <- lprobust(fx$y, fx$x, eval = c(0, 0.25, 0.5))
  expect_null(m$cov.us)
  expect_null(m$cov.rb)
})
