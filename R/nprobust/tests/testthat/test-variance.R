## HC0/HC1/HC2/HC3 standard errors must match the standard sandwich formulas
## evaluated at the corresponding WLS fit.

sandwich_se <- function(y, x, eval, h, vce, p = 1) {
  u <- (x - eval) / h
  sub <- abs(u) <= 1
  xc <- x[sub] - eval
  ys <- y[sub]
  W <- rep(0.5 / h, sum(sub))
  X <- cbind(1, xc)
  XWX_inv <- solve(crossprod(X * sqrt(W)))
  beta <- XWX_inv %*% crossprod(X * W, ys)
  resid <- ys - X %*% beta

  H_mat <- X %*% XWX_inv %*% t(X * W)
  hii <- diag(H_mat)

  nk <- nrow(X); kk <- ncol(X)
  adj <- switch(vce,
                "hc0" = rep(1, nk),
                "hc1" = rep(sqrt(nk / (nk - (p + 1))), nk),
                "hc2" = sqrt(1 / (1 - hii)),
                "hc3" = 1 / (1 - hii))
  e <- as.numeric(adj * resid)
  B_mat <- crossprod(X * W * e, X * W * e)
  vc <- XWX_inv %*% B_mat %*% XWX_inv
  sqrt(vc[1, 1])
}

test_that("HC0/HC1/HC2/HC3 SE match the manual sandwich", {
  fx <- make_fixture()
  for (vce in c("hc0", "hc1", "hc2", "hc3")) {
    est <- lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "uni", vce = vce)
    expect_equal(as.numeric(est$Estimate[, "se.us"]),
                 sandwich_se(fx$y, fx$x, eval = 0, h = 0.3, vce = vce),
                 tolerance = 1e-12,
                 label = paste("vce=", vce))
  }
})

## Cluster-robust variance tests. vce + cluster combinations map to:
##   vce="hc0" + cluster  -> CR0  (raw residuals, no DOF multiplier)
##   vce="hc1" + cluster  -> CR1  (raw residuals, ((n-1)/(n-k))*G/(G-1))
##   vce="hc2" + cluster  -> CR2  (Bell-McCaffrey block-adjusted residuals)
##   vce="hc3" + cluster  -> CR3  (block jackknife (I-H_gg)^{-1}, (G-1)/G multiplier)

cluster_se_manual <- function(y, x, eval, h, cluster, cr_type, p = 1) {
  u   <- (x - eval) / h
  sub <- abs(u) <= 1
  xc  <- x[sub]; ys <- y[sub]; cs <- cluster[sub]
  W   <- rep(0.5 / h, sum(sub))
  sqrtW <- sqrt(W)
  X   <- cbind(1, xc)
  X_std <- X * sqrtW

  XWX_inv <- solve(crossprod(X_std))
  beta <- XWX_inv %*% crossprod(X * W, ys)
  r_raw <- as.numeric(ys - X %*% beta)
  r_std <- r_raw * sqrtW

  clusters <- unique(cs); G <- length(clusters)
  nk <- nrow(X); kk <- ncol(X)

  M <- matrix(0, kk, kk)
  for (g in clusters) {
    ind <- cs == g
    Xg <- X_std[ind, , drop = FALSE]
    rg <- r_std[ind]
    if (cr_type %in% c("CR0", "CR1")) {
      r_adj <- rg
    } else {
      H_gg <- Xg %*% XWX_inv %*% t(Xg)
      H_gg <- (H_gg + t(H_gg)) / 2
      I_H <- diag(nrow(Xg)) - H_gg
      if (cr_type == "CR2") {
        ev <- eigen(I_H, symmetric = TRUE)
        vals <- pmax(ev$values, 1e-12)
        r_adj <- as.numeric(ev$vectors %*% (t(ev$vectors) %*% rg / sqrt(vals)))
      } else {
        r_adj <- drop(solve(I_H, rg))
      }
    }
    score <- crossprod(Xg, r_adj)
    M <- M + tcrossprod(score)
  }
  mult <- switch(cr_type,
                 "CR0" = 1,
                 "CR1" = ((nk - 1) / (nk - kk)) * (G / (G - 1)),
                 "CR2" = 1,
                 "CR3" = (G - 1) / G)
  vc <- XWX_inv %*% (mult * M) %*% XWX_inv
  sqrt(vc[1, 1])
}

test_that("cluster variance: CR1, CR2, CR3 match manual", {
  fx <- make_fixture()
  set.seed(11)
  cl <- sample(1:15, fx$n, replace = TRUE)

  spec <- list(
    list(vce = "cr1", cr = "CR1"),
    list(vce = "cr2", cr = "CR2"),
    list(vce = "cr3", cr = "CR3")
  )
  for (s in spec) {
    est <- lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "uni",
                    vce = s$vce, cluster = cl)
    se_ref <- cluster_se_manual(fx$y, fx$x, eval = 0, h = 0.3, cluster = cl,
                                cr_type = s$cr)
    expect_equal(as.numeric(est$Estimate[, "se.us"]),
                 se_ref, tolerance = 1e-10,
                 label = paste("vce=", s$vce, "->", s$cr))
  }
})

test_that("cluster=1:n with cr1 reduces to HC1-style SE (singleton-cluster sanity)", {
  fx <- make_fixture(n = 200)
  cl <- seq_len(fx$n)   # every obs in its own cluster
  est_cr <- lprobust(fx$y, fx$x, eval = 0, h = 0.3, p = 1, kernel = "uni",
                     vce = "cr1", cluster = cl)
  ## With G == n singleton clusters and CR1's ((n-1)/(n-k)) * (G/(G-1)) =
  ## ((n-1)/(n-k)) * (n/(n-1)) = n/(n-k), the cluster sandwich collapses to
  ## the HC1-style SE on the un-degenerate cluster expansion. Compare to
  ## the manual CR1 reference at G=n.
  se_ref <- cluster_se_manual(fx$y, fx$x, eval = 0, h = 0.3, cluster = cl,
                              cr_type = "CR1")
  expect_equal(as.numeric(est_cr$Estimate[, "se.us"]),
               se_ref, tolerance = 1e-10)
})
