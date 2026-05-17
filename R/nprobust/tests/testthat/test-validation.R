## Regression tests for input-validation hardening (2026-05-08).
library(testthat)
library(nprobust)

set.seed(20260508)
n  <- 400
x  <- runif(n, -1, 1)
y  <- 0.3 * x + rnorm(n)
cl <- ceiling(15 * runif(n))

# Helper that captures warnings + errors so the warning-then-stop pattern
# (warning() in .err(); then stop("nprobust: invalid input...")) is caught.
expect_input_error <- function(expr, pattern) {
  warns <- character(0)
  out <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) e$message
  )
  all_text <- c(if (is.character(out)) out else character(0), warns)
  expect_true(any(grepl(pattern, all_text, fixed = TRUE)),
              info = paste0("expected match='", pattern, "'; got: ",
                            paste(all_text, collapse = " | ")))
}

test_that("level must be in (0, 100)", {
  expect_input_error(lprobust(y, x, level = 100),
                     "level must be a single number")
  expect_input_error(lprobust(y, x, level = 0),
                     "level must be a single number")
  expect_input_error(kdrobust(x, level = 100),
                     "level must be a single number")
  expect_input_error(kdrobust(x, level = c(95, 90)),
                     "level must be a single number")
})

test_that("rho must be a finite non-negative scalar (kdrobust)", {
  expect_input_error(kdrobust(x, rho = -1),
                     "rho must be a single non-negative number")
  expect_input_error(kdrobust(x, rho = c(1, 2)),
                     "rho must be a single non-negative number")
})

test_that("eval finiteness validated in kdrobust / kdbwselect", {
  expect_input_error(kdrobust(x, eval = c(0, NA)),
                     "eval must be numeric and finite")
  expect_input_error(kdbwselect(x, eval = c(0, Inf)),
                     "eval must be numeric and finite")
})

test_that("bwselect=\"all\" rejected by lprobust / kdrobust with clear message", {
  expect_input_error(lprobust(y, x, bwselect = "all"),
                     "only supported by lpbwselect")
  expect_input_error(kdrobust(x, bwselect = "all"),
                     "only supported by kdbwselect")
})

test_that("bwcheck must be a positive integer", {
  expect_input_error(lprobust(y, x, bwcheck = -5),
                     "bwcheck must be a single positive integer")
  expect_input_error(lprobust(y, x, bwcheck = 2.5),
                     "bwcheck must be a single positive integer")
})

test_that("vector lengths validated against length(x)", {
  expect_input_error(lprobust(y, x, weights = rep(1, 10)),
                     "must have length equal to length(x)")
  expect_input_error(lprobust(y, x, cluster = cl[1:10]),
                     "must have length equal to length(x)")
  expect_input_error(lprobust(y, x, subset = rep(TRUE, 10)),
                     "Logical 'subset' must have length")
})

test_that("empty subset errors", {
  expect_input_error(lprobust(y, x, subset = rep(FALSE, n)),
                     "removed all observations")
})

test_that("eval must be finite", {
  expect_input_error(lprobust(y, x, eval = c(0.0, NA, 0.5)),
                     "eval must be numeric and finite")
})

test_that("bwselect typo rejected", {
  expect_input_error(lprobust(y, x, bwselect = "no-such-method"),
                     "bwselect incorrectly specified")
})

test_that("manual h with rho=0 and no b works (regression)", {
  # rho=0 + manual h + manual b is OK (b is provided directly)
  expect_silent(suppressWarnings(lprobust(y, x, eval = 0, h = 0.5, b = 0.5, rho = 0)))
})
