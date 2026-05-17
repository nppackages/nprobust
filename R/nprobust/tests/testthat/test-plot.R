## S3 plot methods should dispatch to nprobust.plot and return ggplot objects.

test_that("plot.lprobust returns a ggplot", {
  fx <- make_fixture(n = 300)
  m  <- lprobust(fx$y, fx$x, neval = 10)
  g  <- plot(m)
  expect_s3_class(g, "ggplot")
})

test_that("plot.kdrobust returns a ggplot", {
  fx <- make_fixture(n = 300)
  kd <- kdrobust(fx$x, neval = 10)
  g  <- plot(kd)
  expect_s3_class(g, "ggplot")
})

test_that("nprobust.plot accepts multiple series", {
  fx <- make_fixture(n = 400)
  m0 <- lprobust(fx$y[fx$x <= 0], fx$x[fx$x <= 0], eval = seq(-0.8, 0, length.out = 10))
  m1 <- lprobust(fx$y[fx$x >  0], fx$x[fx$x >  0], eval = seq(0, 0.8, length.out = 10))
  g  <- nprobust.plot(m0, m1, legendGroups = c("Left", "Right"))
  expect_s3_class(g, "ggplot")
})
