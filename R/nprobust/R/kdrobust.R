kdrobust <- function(x, eval=NULL, neval=NULL, h=NULL, b=NULL, rho=1, kernel="epa",
                    bwselect=NULL, bwcheck=21, imsegrid=30, level=95, subset = NULL,
                    data = NULL) {

  if (!is.null(data)) {
    mc <- match.call()
    caller_env <- parent.frame()
    .lookup <- function(arg) {
      expr <- mc[[arg]]
      if (is.null(expr)) return(NULL)
      eval(expr, envir = data, enclos = caller_env)
    }
    x <- .lookup("x")
    if ("subset" %in% names(mc)) subset <- .lookup("subset")
  }

  p <- 2
  deriv <- 0
  if (!is.null(subset)) x <- x[subset]

  ## UX prechecks: catch invalid h/b/bwcheck/imsegrid with clear messages
  .nperrs <- character()
  .bwbad <- function(v, name) {
    if (is.null(v))         return(NULL)
    if (!is.numeric(v))     return(paste0(name, " must be numeric."))
    if (any(!is.finite(v))) return(paste0(name, " contains non-finite values (NA/Inf)."))
    if (any(v <= 0))        return(paste0(name, " must be strictly positive."))
    NULL
  }
  .nperrs <- c(.nperrs, .bwbad(h, "h"), .bwbad(b, "b"))
  if (!is.null(bwcheck) && (!is.numeric(bwcheck) || length(bwcheck) != 1 || !is.finite(bwcheck) || bwcheck <= 0))
    .nperrs <- c(.nperrs, "bwcheck must be a single positive finite integer.")
  if (!is.numeric(imsegrid) || length(imsegrid) != 1 || !is.finite(imsegrid) || imsegrid <= 0)
    .nperrs <- c(.nperrs, "imsegrid must be a single positive finite integer.")
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) || level <= 0 || level >= 100)
    .nperrs <- c(.nperrs, "level must be a single number in (0, 100).")
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho) || rho < 0)
    .nperrs <- c(.nperrs, "rho must be a single non-negative number.")
  if (!is.null(eval)) {
    if (!is.numeric(eval) || any(!is.finite(eval)))
      .nperrs <- c(.nperrs, "eval must be numeric and finite.")
    if (length(eval) == 0L)
      .nperrs <- c(.nperrs, "eval must have at least one element.")
  }
  if (length(.nperrs) > 0) {
    for (.m in .nperrs) warning(.m, call. = FALSE)
    stop("nprobust: invalid input (see warnings above).", call. = FALSE)
  }

  na.ok <- complete.cases(x)
  x <- x[na.ok]
  
  x.min <- min(x);  x.max <- max(x)
  N <- length(x)

  if (!is.null(bwcheck)) {
    if (bwcheck > N) {
      warning("bwcheck (", bwcheck, ") is larger than the sample size (", N,
              "); reducing bwcheck to N.")
      bwcheck <- N
    }
  }

  if (is.null(eval)) {
    if (is.null(neval)) {
      qseq <- seq(0.1,0.9,length.out=30)
      eval <- quantile(x, qseq)
      #eval <- unique(x)
      #qseq <- seq(0,1,1/(20+1))
      #eval <- quantile(x, qseq[2:(length(qseq)-1)])
      #eval <- seq(x.min, x.max, length.out=30)
    }
    else {
      qseq <- seq(0.1,0.9,length.out=neval)
      eval <- quantile(x, qseq)
      #eval <- seq(x.min,x.max,length.out=neval)
      #qseq <- seq(0,1,1/(neval+1))
      #eval <- quantile(x, qseq[2:(length(qseq)-1)])
      #eval <- seq(x.min, x.max, length.out=neval)
    }
  }
  neval <- length(eval)

  ## Precheck (continued): h/b length must be 1 or neval
  .nperrs2 <- character()
  if (!is.null(h) && length(h) != 1L && length(h) != neval)
    .nperrs2 <- c(.nperrs2, paste0("h must have length 1 or neval (=", neval, ")."))
  if (!is.null(b) && length(b) != 1L && length(b) != neval)
    .nperrs2 <- c(.nperrs2, paste0("b must have length 1 or neval (=", neval, ")."))
  if (length(.nperrs2) > 0) {
    for (.m in .nperrs2) warning(.m, call. = FALSE)
    stop("nprobust: invalid input (see warnings above).", call. = FALSE)
  }

  if (is.null(h) & is.null(bwselect) & neval==1) bwselect="mse-dpi"
  if (is.null(h) & is.null(bwselect) & neval>1)  bwselect="imse-dpi"
  
  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)

  #####################################################   CHECK ERRORS
  if (!(kernel %in% c("epa","epanechnikov","uni","uniform"))) {
    stop("kernel incorrectly specified. Supported kernels for kdrobust: epa, uni.")
  }
  if (kernel %in% c("epanechnikov")) kernel <- "epa"
  if (kernel %in% c("uniform"))      kernel <- "uni"

  if (identical(bwselect, "all"))
    stop("bwselect=\"all\" is only supported by kdbwselect; kdrobust requires a single method (e.g. \"mse-dpi\", \"imse-dpi\", \"ce-dpi\").",
         call. = FALSE)

  if (!is.null(h)) bwselect <- "Manual"

  kernel.type <- if (kernel == "epa") "Epanechnikov" else "Uniform"

  if (!is.null(h) & rho>0 & is.null(b)) {
    #rho <- rep(1,neval)
    b <- h/rho
  }
  #if (!is.null(h) & !is.null(rho) ) b <- h/rho
  if (is.null(h)) {
      kdbws <- kdbwselect(x=x, eval=eval, bwselect=bwselect, bwcheck=bwcheck, imsegrid=imsegrid, kernel=kernel)
      h <- kdbws$bws[,2]
      b <- kdbws$bws[,3]
      if (rho>0) b <- h/rho
      rho <- h/b  
  }
  
  if (length(h)==1 & neval>1) {
    h <- rep(h,neval)
    b <- rep(b,neval)
    rho   <- h/b
  }

  
  Estimate<-matrix(NA,neval,8)
  colnames(Estimate)<-c("eval","h","b","N","tau.us","tau.bc","se.us","se.rb")
  
  for (i in 1:neval) {
    
    if (!is.null(bwcheck)) {
      bw.min   <- sort(abs(x-eval[i]))[bwcheck]
      h[i]     <- max(h[i], bw.min)
      b[i]     <- max(b[i], bw.min)
      rho[i]   <- h[i]/b[i]
    }
    
    u   <- (x-eval[i])/h[i]
    K.d <- kd.K.fun(u        , v=p,   r=deriv, kernel=kernel)
    L.r <- kd.K.fun(rho[i]*u , v=p+2, r=p,     kernel=kernel)
  
    K     <- K.d$Kx
    M     <- K - rho[i]^(1+p)*L.r$Kx*L.r$k.v
    f.us  <- mean(K)/h[i]
    f.bc  <- mean(M)/h[i]
    se.us <- sqrt((mean((K^2)) - mean(K)^2)/(N*h[i]^2))
    se.rb <- sqrt((mean((M^2)) - mean(M)^2)/(N*h[i]^2))

    if (kernel == "gau") {
      eN <- N
    } else {
      eN <- sum(abs(x - eval[i]) <= max(h[i], b[i]))
    }
    
    Estimate[i,] <- c(eval[i], h[i], b[i], eN, f.us, f.bc, se.us, se.rb) 
  }
  out<-list(Estimate=Estimate, opt=list(p=p, kernel=kernel.type, n=N, neval=neval, bwselect=bwselect))
  out$call <- match.call()
  class(out) <- "kdrobust"
  return(out)
}

print.kdrobust <- function(x,...){
  cat("Call: kdrobust\n\n")

  cat(paste("Sample size (n)                            =     ", x$opt$n,        "\n", sep=""))
  cat(paste("Kernel order for point estimation (p)      =     ", x$opt$p,        "\n", sep=""))
  cat(paste("Kernel function                            =     ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                           =     ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  invisible(x)
}

summary.kdrobust <- function(object, alpha = 0.05, sep = 5, ...) {
  z    <- qnorm(1 - alpha / 2)
  CI_l <- object$Estimate[, "tau.bc"] - object$Estimate[, "se.rb"] * z
  CI_r <- object$Estimate[, "tau.bc"] + object$Estimate[, "se.rb"] * z

  out <- list(opt      = object$opt,
              Estimate = object$Estimate,
              alpha    = alpha,
              sep      = sep,
              CI_l     = CI_l,
              CI_r     = CI_r)
  class(out) <- "summary.kdrobust"
  out
}

print.summary.kdrobust <- function(x, ...) {
  cat("Call: kdrobust\n\n")

  cat(paste("Sample size (n)                            =     ", x$opt$n,        "\n", sep=""))
  cat(paste("Kernel order for point estimation (p)      =     ", x$opt$p,        "\n", sep=""))
  cat(paste("Kernel function                            =     ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth selection method                 =     ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  alpha <- x$alpha
  sep   <- x$sep

  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  cat(format(" ", width= 14 ))
  cat(format(" ", width= 10 ))
  cat(format(" ", width= 8  ))
  cat(format("Point", width= 10, justify="right"))
  cat(format("Std." , width= 10, justify="right"))
  cat(format("Robust B.C.", width=25, justify="centre"))
  cat("\n")

  cat(format("eval"            , width=14, justify="right"))
  cat(format("bw"              , width=10, justify="right"))
  cat(format("Eff.n"           , width=8 , justify="right"))
  cat(format("Est."            , width=10, justify="right"))
  cat(format("Error"           , width=10, justify="right"))
  cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep="")
             , width=25, justify="centre"))
  cat("\n")

  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  for (j in 1:nrow(x$Estimate)) {
    cat(format(toString(j), width=4))
    cat(format(sprintf("%3.3f", x$Estimate[j, "eval"]), width=10, justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[j, "h"])  , width=10, justify="right"))
    cat(format(sprintf("%3.0f", x$Estimate[j, "N"])  , width=8 , justify="right"))
    cat(format(sprintf("%3.3f", x$Estimate[j, "tau.us"]) , width=10, justify="right"))
    cat(format(paste(sprintf("%3.3f", x$Estimate[j, "se.us"]), sep=""), width=10, justify="right"))
    cat(format(paste("[", sprintf("%3.3f", x$CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", x$CI_r[j]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
      cat(paste(rep("-", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
    }
  }

  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")

  invisible(x)
}