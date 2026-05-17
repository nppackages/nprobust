## Tests for cluster-vce surface: cr1/cr2/cr3 acceptance, nn+cluster
## redirects to cr1, hc-cluster aliases match cr* equivalents, and cr*
## without cluster falls back to hc* with a warning.

make_cluster_fixture <- function(n = 800, n_clusters = 30, seed = 11) {
  set.seed(seed)
  x  <- sort(runif(n, -1, 1))
  y  <- 1 + 2 * x + 3 * x^2 + rnorm(n, sd = 0.4)
  cl <- ceiling(n_clusters * runif(n))
  list(x = x, y = y, cl = cl)
}

test_that("cr1/cr2/cr3 require cluster (without cluster they warn + fall back)", {
  d <- make_cluster_fixture()
  expect_warning(
    m1 <- lprobust(d$y, d$x, eval = 0, vce = "cr1"),
    "requires a cluster variable"
  )
  expect_warning(
    m2 <- lprobust(d$y, d$x, eval = 0, vce = "cr2"),
    "requires a cluster variable"
  )
  expect_warning(
    m3 <- lprobust(d$y, d$x, eval = 0, vce = "cr3"),
    "requires a cluster variable"
  )
  ## After fallback: cr1->hc1, cr2->hc2, cr3->hc3 produces the same
  ## numbers as if the user had passed hc1/hc2/hc3 directly.
  m_hc1 <- lprobust(d$y, d$x, eval = 0, vce = "hc1")
  m_hc2 <- lprobust(d$y, d$x, eval = 0, vce = "hc2")
  m_hc3 <- lprobust(d$y, d$x, eval = 0, vce = "hc3")
  expect_equal(m1$Estimate[, "se.us"], m_hc1$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m2$Estimate[, "se.us"], m_hc2$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m3$Estimate[, "se.us"], m_hc3$Estimate[, "se.us"], tolerance = 1e-12)
})

test_that("vce='nn' + cluster silently defaults to cr1", {
  d <- make_cluster_fixture()
  ## No warning expected for the silent default.
  expect_silent(
    m_default <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "nn")
  )
  ## Same as explicit cr1.
  m_cr1 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr1")
  expect_equal(m_default$Estimate[, "se.us"], m_cr1$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m_default$Estimate[, "se.rb"], m_cr1$Estimate[, "se.rb"], tolerance = 1e-12)
})

test_that("hc{0,1,2,3} + cluster warns and remaps to cr1/cr1/cr2/cr3", {
  d <- make_cluster_fixture()
  expect_warning(
    m_hc0 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "hc0"),
    "Switching to vce='cr1'"
  )
  expect_warning(
    m_hc1 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "hc1"),
    "Switching to vce='cr1'"
  )
  expect_warning(
    m_hc2 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "hc2"),
    "Switching to vce='cr2'"
  )
  expect_warning(
    m_hc3 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "hc3"),
    "Switching to vce='cr3'"
  )
  m_cr1 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr1")
  m_cr2 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr2")
  m_cr3 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr3")
  expect_equal(m_hc0$Estimate[, "se.us"], m_cr1$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m_hc1$Estimate[, "se.us"], m_cr1$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m_hc2$Estimate[, "se.us"], m_cr2$Estimate[, "se.us"], tolerance = 1e-12)
  expect_equal(m_hc3$Estimate[, "se.us"], m_cr3$Estimate[, "se.us"], tolerance = 1e-12)
})

test_that("CR1, CR2, CR3 produce different SEs (sanity)", {
  d <- make_cluster_fixture()
  m1 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr1")
  m2 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr2")
  m3 <- lprobust(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr3")
  expect_false(isTRUE(all.equal(m1$Estimate[, "se.us"], m2$Estimate[, "se.us"])))
  expect_false(isTRUE(all.equal(m1$Estimate[, "se.us"], m3$Estimate[, "se.us"])))
  expect_false(isTRUE(all.equal(m2$Estimate[, "se.us"], m3$Estimate[, "se.us"])))
  ## CR3 >= CR1 elementwise (block jackknife is conservative)
  expect_true(m3$Estimate[, "se.us"] >= m1$Estimate[, "se.us"])
})

test_that("invalid vce names error", {
  d <- make_cluster_fixture()
  expect_error(
    suppressWarnings(lprobust(d$y, d$x, eval = 0, vce = "garbage")),
    "vce incorrectly specified|invalid input"
  )
})

test_that("lpbwselect mirrors lprobust vce surface", {
  d <- make_cluster_fixture()
  ## lpbwselect with cr1 (cluster supplied) returns same h as hc1 (which
  ## also gets remapped to cr1 with a warning).
  expect_silent(
    bw_cr1 <- lpbwselect(d$y, d$x, eval = 0, cluster = d$cl, vce = "cr1")
  )
  expect_warning(
    bw_hc1 <- lpbwselect(d$y, d$x, eval = 0, cluster = d$cl, vce = "hc1"),
    "Switching to vce='cr1'"
  )
  expect_equal(bw_cr1$bws[, "h"], bw_hc1$bws[, "h"], tolerance = 1e-12)
  expect_equal(bw_cr1$bws[, "b"], bw_hc1$bws[, "b"], tolerance = 1e-12)
})
