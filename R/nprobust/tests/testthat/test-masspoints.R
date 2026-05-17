## masspoints="check" warns when too few unique x values fall inside the
## bandwidth window. masspoints="off" silences the warning.

test_that("masspoints='check' warns on discrete running variable", {
  set.seed(5)
  # discrete grid with few unique values
  x <- sample(c(-1, -0.5, 0, 0.5, 1), 300, replace = TRUE)
  y <- x + rnorm(300, sd = 0.1)
  expect_warning(
    lprobust(y, x, eval = 0, h = 0.75, p = 1, kernel = "uni",
             vce = "hc0", masspoints = "check"),
    "unique x values within bandwidth"
  )
})

test_that("masspoints='off' suppresses the warning", {
  set.seed(5)
  x <- sample(c(-1, -0.5, 0, 0.5, 1), 300, replace = TRUE)
  y <- x + rnorm(300, sd = 0.1)
  expect_silent(
    lprobust(y, x, eval = 0, h = 0.75, p = 1, kernel = "uni",
             vce = "hc0", masspoints = "off")
  )
})

test_that("masspoints rejects invalid values", {
  fx <- make_fixture(n = 100)
  expect_error(
    lprobust(fx$y, fx$x, eval = 0, h = 0.3, masspoints = "bogus"),
    "masspoints must be one of"
  )
})

test_that("masspoints='check' stays silent with abundant unique x", {
  fx <- make_fixture(n = 800)
  expect_silent(
    lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "uni", vce = "hc0")
  )
})
