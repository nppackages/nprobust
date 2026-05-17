test_that("fixed local-polynomial numerical baseline is stable", {
  fx <- make_fixture()
  est <- lprobust(fx$y, fx$x, eval = c(-0.3, 0, 0.4), h = 0.25,
                  p = 2, kernel = "tri", vce = "hc0")

  expect_equal(
    unname(as.matrix(est$Estimate[, c("tau.us", "tau.bc", "se.rb")])),
    matrix(c(
      0.638434993601916, 0.643508307149877, 0.0284533673080986,
      1.00801554482691, 1.00276733583938, 0.0281615189408019,
      2.21898254529333, 2.22208193444536, 0.0316295999817565
    ), ncol = 3, byrow = TRUE),
    tolerance = 1e-10
  )
})

test_that("fixed density and bandwidth numerical baselines are stable", {
  fx <- make_fixture()
  kd <- kdrobust(fx$x, eval = c(-0.5, 0, 0.7), h = 0.35, kernel = "epa")
  bw <- lpbwselect(fx$y, fx$x, eval = 0, bwselect = "mse-dpi")

  expect_equal(
    unname(as.matrix(kd$Estimate[, c("tau.us", "tau.bc", "se.rb")])),
    matrix(c(
      0.478701016953545, 0.464121681225866, 0.0306029948861648,
      0.524201444786818, 0.534315088652292, 0.0326839226564865,
      0.51741170181274, 0.535936195748363, 0.0323836078172713
    ), ncol = 3, byrow = TRUE),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(bw$bws[, c("h", "b")]),
    c(0.150319401795546, 0.99775844020769),
    tolerance = 1e-10
  )
})
