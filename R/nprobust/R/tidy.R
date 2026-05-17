## broom-compatible tidy()/glance() methods.
##
## Registered in NAMESPACE via `S3method(broom::tidy, <class>)` so they
## dispatch when broom is loaded, without requiring broom at install time.
##
## Convention:
##   - `estimate`   = tau.us  (classical/conventional point estimate)
##   - `std.error`  = se.us   (its SE)
##   - `conf.low`/`conf.high` = robust bias-corrected CI:
##        tau.bc +/- z * se.rb
##   - additional columns tau.bc, se.rb, h, b, n.eff are kept for
##     downstream convenience.

tidy.lprobust <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  est <- as.data.frame(x$Estimate)
  out <- data.frame(
    eval       = est$eval,
    estimate   = est$tau.us,
    std.error  = est$se.us,
    tau.bc     = est$tau.bc,
    se.rb      = est$se.rb,
    h          = est$h,
    b          = est$b,
    n.eff      = est$N,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low  <- est$tau.bc - z * est$se.rb
    out$conf.high <- est$tau.bc + z * est$se.rb
  }
  out
}

glance.lprobust <- function(x, ...) {
  data.frame(
    n        = x$opt$n,
    neval    = x$opt$neval,
    p        = x$opt$p,
    q        = x$opt$q,
    deriv    = x$opt$deriv,
    kernel   = x$opt$kernel,
    bwselect = x$opt$bwselect,
    stringsAsFactors = FALSE
  )
}

tidy.kdrobust <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  est <- as.data.frame(x$Estimate)
  out <- data.frame(
    eval       = est$eval,
    estimate   = est$tau.us,
    std.error  = est$se.us,
    tau.bc     = est$tau.bc,
    se.rb      = est$se.rb,
    h          = est$h,
    b          = est$b,
    n.eff      = est$N,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low  <- est$tau.bc - z * est$se.rb
    out$conf.high <- est$tau.bc + z * est$se.rb
  }
  out
}

glance.kdrobust <- function(x, ...) {
  data.frame(
    n        = x$opt$n,
    neval    = x$opt$neval,
    kernel   = x$opt$kernel,
    bwselect = x$opt$bwselect,
    stringsAsFactors = FALSE
  )
}

tidy.lpbwselect <- function(x, ...) {
  bws <- as.data.frame(x$bws)
  rownames(bws) <- NULL
  bws
}

glance.lpbwselect <- function(x, ...) {
  data.frame(
    n        = x$opt$n,
    neval    = x$opt$neval,
    p        = x$opt$p,
    q        = x$opt$q,
    deriv    = x$opt$deriv,
    kernel   = x$opt$kernel,
    bwselect = x$opt$bwselect,
    stringsAsFactors = FALSE
  )
}

tidy.kdbwselect <- function(x, ...) {
  bws <- as.data.frame(x$bws)
  rownames(bws) <- NULL
  bws
}

glance.kdbwselect <- function(x, ...) {
  data.frame(
    n        = x$opt$n,
    neval    = x$opt$neval,
    kernel   = x$opt$kernel,
    bwselect = x$opt$bwselect,
    stringsAsFactors = FALSE
  )
}
