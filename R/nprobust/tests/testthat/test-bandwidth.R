## All bandwidth selectors must produce finite, strictly positive bandwidths.

test_that("all lpbwselect methods produce finite positive bandwidths", {
  fx <- make_fixture()
  for (bw in c("mse-dpi", "mse-rot", "imse-dpi", "imse-rot", "ce-dpi", "ce-rot")) {
    b <- lpbwselect(fx$y, fx$x, eval = 0, bwselect = bw)
    hh <- b$bws[, "h"]; bb <- b$bws[, "b"]
    expect_true(all(is.finite(hh) & hh > 0), label = paste("h for", bw))
    expect_true(all(is.finite(bb) & bb > 0), label = paste("b for", bw))
  }
})

test_that("interior=TRUE and interior=FALSE both succeed", {
  fx <- make_fixture()
  b1 <- lpbwselect(fx$y, fx$x, eval = 0, bwselect = "mse-dpi", interior = TRUE)
  b2 <- lpbwselect(fx$y, fx$x, eval = 0, bwselect = "mse-dpi", interior = FALSE)
  expect_true(is.finite(b1$bws[, "h"]))
  expect_true(is.finite(b2$bws[, "h"]))
})

test_that("bwcheck > N triggers a warning (bug #8 regression)", {
  fx <- make_fixture(n = 10)
  expect_warning(
    lprobust(fx$y, fx$x, eval = 0.5, bwcheck = 50),
    "bwcheck .* is larger than the sample size"
  )
})

test_that("rho=0 with manual h and no b is an error (bug #5 regression)", {
  fx <- make_fixture()
  expect_error(
    lprobust(fx$y, fx$x, eval = 0.5, h = 0.1, rho = 0),
    "b must also be provided"
  )
})
