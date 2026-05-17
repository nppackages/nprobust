lprobust = function(y, x, eval=NULL, neval=NULL, p=NULL, deriv=NULL, h=NULL, b=NULL, rho=1,
                    kernel="epa", bwselect=NULL, bwcheck=21, bwregul=1, imsegrid=30, vce="nn", covgrid = FALSE,
                    cluster=NULL, nnmatch=3, level=95, interior = FALSE, subset = NULL,
                    weights = NULL, masspoints = "check", data = NULL) {

  if (!is.null(data)) {
    mc <- match.call()
    caller_env <- parent.frame()
    .lookup <- function(arg) {
      expr <- mc[[arg]]
      if (is.null(expr)) return(NULL)
      eval(expr, envir = data, enclos = caller_env)
    }
    y <- .lookup("y")
    x <- .lookup("x")
    if ("cluster" %in% names(mc)) cluster <- .lookup("cluster")
    if ("weights" %in% names(mc)) weights <- .lookup("weights")
    if ("subset"  %in% names(mc)) subset  <- .lookup("subset")
  }

  ## Validate aux-vector lengths against length(x) BEFORE subset filtering,
  ## so wrong-length inputs error explicitly instead of silently recycling.
  .n_orig <- length(x)
  if (length(y) != .n_orig)
    stop(sprintf("'y' and 'x' must have equal length (got y=%d, x=%d).", length(y), .n_orig), call. = FALSE)
  if (!is.null(cluster) && length(cluster) != .n_orig)
    stop(sprintf("'cluster' must have length equal to length(x) (got %d, expected %d).", length(cluster), .n_orig), call. = FALSE)
  if (!is.null(weights) && length(weights) != .n_orig)
    stop(sprintf("'weights' must have length equal to length(x) (got %d, expected %d).", length(weights), .n_orig), call. = FALSE)
  if (!is.null(subset)) {
    if (is.logical(subset)) {
      if (length(subset) != .n_orig)
        stop(sprintf("Logical 'subset' must have length equal to length(x) (got %d, expected %d).", length(subset), .n_orig), call. = FALSE)
    } else if (is.numeric(subset)) {
      if (any(!is.finite(subset)) || any(subset < 1) || any(subset > .n_orig) || any(subset != round(subset)))
        stop(sprintf("Numeric 'subset' must contain integer indices in 1..%d.", .n_orig), call. = FALSE)
    } else {
      stop("'subset' must be logical or integer.", call. = FALSE)
    }
  }

  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
    if (!is.null(cluster)) cluster <- cluster[subset]
    if (!is.null(weights)) weights <- weights[subset]
    if (length(x) == 0L) stop("'subset' removed all observations.", call. = FALSE)
  }

  ## UX prechecks: catch invalid h/b/bwcheck/imsegrid/weights/cluster/eval/level/rho/etc.
  .nperrs <- character()
  .bwbad <- function(v, name) {
    if (is.null(v))         return(NULL)
    if (!is.numeric(v))     return(paste0(name, " must be numeric."))
    if (any(!is.finite(v))) return(paste0(name, " contains non-finite values (NA/Inf)."))
    if (any(v <= 0))        return(paste0(name, " must be strictly positive."))
    NULL
  }
  .nperrs <- c(.nperrs, .bwbad(h, "h"), .bwbad(b, "b"))
  if (!is.null(bwcheck) && (!is.numeric(bwcheck) || length(bwcheck) != 1 || !is.finite(bwcheck) || bwcheck <= 0 || bwcheck != round(bwcheck)))
    .nperrs <- c(.nperrs, "bwcheck must be a single positive integer.")
  if (!is.numeric(imsegrid) || length(imsegrid) != 1 || !is.finite(imsegrid) || imsegrid <= 0 || imsegrid != round(imsegrid))
    .nperrs <- c(.nperrs, "imsegrid must be a single positive integer.")
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
  if (!is.null(weights)) {
    if (!is.numeric(weights))               .nperrs <- c(.nperrs, "weights must be numeric.")
    if (length(weights) != length(x))       .nperrs <- c(.nperrs, "weights length must equal length(x) (after subset).")
    if (sum(weights, na.rm = TRUE) <= 0)    .nperrs <- c(.nperrs, "weights must have a strictly positive sum.")
  }
  if (length(.nperrs) > 0) {
    for (.m in .nperrs) warning(.m, call. = FALSE)
    stop("nprobust: invalid input (see warnings above).", call. = FALSE)
  }

  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(cluster)) na.ok <- na.ok & complete.cases(cluster)
  if (!is.null(weights)) na.ok <- na.ok & complete.cases(weights) & weights >= 0

  x <- x[na.ok]
  y <- y[na.ok]
  if (!is.null(cluster)) cluster <- cluster[na.ok]
  if (!is.null(weights)) weights <- weights[na.ok]
  if (is.null(weights))  weights <- rep(1, length(x))
  
  if (!is.null(deriv) & is.null(p)) p <- deriv+1
  if (is.null(p))         p <- 1
  if (is.null(deriv)) deriv <- 0
  q <- p+1
  
  x.max <- max(x); x.min <- min(x)
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
      #eval <- unique(x)
      #qseq <- seq(0,1,1/(20+1))
      #eval <- quantile(x, qseq[2:(length(qseq)-1)])
      eval <- seq(x.min, x.max, length.out=30)
    }
    else {
      #eval <- seq(x.min,x.max,length.out=neval)
      #qseq <- seq(0,1,1/(neval+1))
      #eval <- quantile(x, qseq[2:(length(qseq)-1)])
      eval <- seq(x.min, x.max, length.out=neval)
    }
  }
  neval <- length(eval)

  ## Precheck (continued): h/b length must be 1 or neval, otherwise indexing later silently recycles
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
  
  if (vce=="nn") {
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
    if (!is.null(cluster)) cluster <- cluster[order.x]
    weights <- weights[order.x]
  }

  if (!(masspoints %in% c("check", "off"))) {
    stop("masspoints must be one of \"check\" or \"off\".")
  }

  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)
  vce      <- tolower(vce)

  #####################################################   CHECK ERRORS
  exit <- 0
  .err <- function(msg) { warning(msg, call. = FALSE); exit <<- 1 }

  valid_kernels  <- c("gau","gaussian","uni","uniform","tri","triangular","epa","epanechnikov","")
  valid_vce      <- c("nn","hc0","hc1","hc2","hc3","cr1","cr2","cr3","")
  valid_bwselect <- c("mse-dpi","mse-rot","imse-dpi","imse-rot","ce-dpi","ce-rot","manual","")

  if (!(kernel %in% valid_kernels)) .err("kernel incorrectly specified")
  if (!(vce    %in% valid_vce))     .err("vce incorrectly specified")
  ## bwselect may be character(0) (after tolower(NULL)) when h is supplied
  ## manually; in that case it gets set to "Manual" later. Skip check then.
  if (!is.null(bwselect) && length(bwselect) > 0L && nzchar(bwselect) &&
      !(bwselect %in% valid_bwselect)) {
    if (identical(bwselect, "all")) {
      .err("bwselect=\"all\" is only supported by lpbwselect; lprobust requires a single method (e.g. \"mse-dpi\", \"imse-dpi\", \"ce-dpi\").")
    } else {
      .err(paste0("bwselect incorrectly specified (received '", bwselect, "')"))
    }
  }

  if (p < 0 || deriv < 0 || nnmatch <= 0)
    .err("p, deriv must be >=0 and nnmatch must be >0")
  if (deriv > p)
    .err("deriv must be <= p")
  ## integer roundness check
  if (p > 0 && round(p) != p)             .err("p must be an integer")
  if (deriv >= 0 && round(deriv) != deriv) .err("deriv must be an integer")
  if (nnmatch > 0 && round(nnmatch) != nnmatch) .err("nnmatch must be an integer")

  if (exit > 0) stop("nprobust: invalid input (see warnings above).", call. = FALSE)

  # Cluster vce validation: with cluster, only cr1/cr2/cr3 are valid.
  # Without cluster, cr1/cr2/cr3 fall back to hc1/hc2/hc3.
  if (!is.null(cluster)) {
    if (vce %in% c("nn", "")) {
      vce <- "cr1"  # silent default
    } else if (vce %in% c("hc0", "hc1")) {
      warning(paste0("vce='", vce, "' is not a cluster option. Switching to vce='cr1'."), call. = FALSE)
      vce <- "cr1"
    } else if (vce == "hc2") {
      warning("vce='hc2' is not a cluster option. Switching to vce='cr2'.", call. = FALSE)
      vce <- "cr2"
    } else if (vce == "hc3") {
      warning("vce='hc3' is not a cluster option. Switching to vce='cr3'.", call. = FALSE)
      vce <- "cr3"
    }
  } else {
    if (vce == "cr1") {
      warning("vce='cr1' requires a cluster variable. Falling back to vce='hc1'.", call. = FALSE)
      vce <- "hc1"
    } else if (vce == "cr2") {
      warning("vce='cr2' requires a cluster variable. Falling back to vce='hc2'.", call. = FALSE)
      vce <- "hc2"
    } else if (vce == "cr3") {
      warning("vce='cr3' requires a cluster variable. Falling back to vce='hc3'.", call. = FALSE)
      vce <- "hc3"
    }
  }

  # User-facing display label.
  vce_type <- "NN"
  if (vce == "hc0") vce_type <- "HC0"
  if (vce == "hc1") vce_type <- "HC1"
  if (vce == "hc2") vce_type <- "HC2"
  if (vce == "hc3") vce_type <- "HC3"
  if (vce == "cr1") vce_type <- "CR1"
  if (vce == "cr2") vce_type <- "CR2"
  if (vce == "cr3") vce_type <- "CR3"

    if (!is.null(h)) bwselect <- "Manual"
    if (!is.null(h) & rho>0 & is.null(b)) {
        b <- h/rho
    }
    if (!is.null(h) & rho==0 & is.null(b)) {
        stop("When h is provided and rho=0, b must also be provided (b cannot be computed as h/rho).")
    }
    
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

    
  ############################################################################################
    #print("Preparing data.")
    if (is.null(h)) {
        lpbws <- lpbwselect(y=y, x=x,  eval=eval, deriv=deriv, p=p, vce=vce,
                            cluster=cluster, bwselect=bwselect, interior=interior,
                            kernel=kernel, bwcheck=bwcheck, bwregul=bwregul,
                            imsegrid=imsegrid, subset=NULL,
                            weights=weights, masspoints="off")
        h     <- lpbws$bws[,2]
        b     <- lpbws$bws[,3]
        if (rho>0) b <- h/rho
        rho   <- h/b
    }

  # Internal mapping: cr1/cr2/cr3 -> hc1/hc2/hc3. Done AFTER the lpbwselect
  # call so that lpbwselect receives the user-facing cr* label and can do
  # its own (warning-free) normalization without double-firing the
  # cluster-vce warning the user has already seen.
  if (vce == "cr1") vce <- "hc1"
  if (vce == "cr2") vce <- "hc2"
  if (vce == "cr3") vce <- "hc3"
  
    if (length(h)==1 & neval>1) {
      h <- rep(h,neval)
      b <- rep(b,neval)
      rho   <- h/b
    }
  
  dups <- dupsid <- 0
  if (vce=="nn") {
    # x is sorted at this point; rle() gives run-length encoding of duplicates.
    runs <- rle(x)
    dups   <- rep.int(runs$lengths, runs$lengths)
    dupsid <- sequence(runs$lengths)
  }
  
  cov.p = NULL
  ####################################################
  Estimate=matrix(NA,neval,8)
  colnames(Estimate)=c("eval","h","b","N","tau.us","tau.bc","se.us","se.rb")
  
  
  for (i in 1:neval) {
    
    if (!is.null(bwcheck)) {
      bw.min   <- sort(abs(x-eval[i]))[bwcheck]
      #nh       <- sum(abs(x-eval[i]) <= h)
      #nb       <- sum(abs(x-eval[i]) <= b)
      h[i]     <- max(h[i], bw.min)
      b[i]     <- max(b[i], bw.min)
    }
    
    w.h   <- W.fun((x-eval[i])/h[i], kernel)/h[i] * weights
    w.b   <- W.fun((x-eval[i])/b[i], kernel)/b[i] * weights
    ind.h <- w.h>0;  ind.b <- w.b>0
    N.h   <- sum(ind.h);  N.b <- sum(ind.b)
    ind   <- ind.b
    if (h[i]>b[i]) ind <- ind.h

    if (masspoints == "check") {
      n_unique <- length(unique(x[ind.h]))
      if (n_unique < (p + 5)) {
        warning(sprintf(
          "Only %d unique x values within bandwidth at eval=%.4f (p+5=%d); local polynomial may be unreliable. Set masspoints=\"off\" to silence.",
          n_unique, eval[i], p + 5))
      }
    }
    
    #if (N.h.l<5 | N.h.r<5 | N.b.l<5 | N.b.r<5){
    #  stop("Not enough observations to perform calculations")
    #  exit(1)
    #}
    
    eN  <- sum(ind)
    eY  <- y[ind]
    eX  <- x[ind]
    W.h <- w.h[ind]
    W.b <- w.b[ind]
    
    eC = NULL
    if (!is.null(cluster)) eC = cluster[ind]
    
    edups <- edupsid <- 0	
    if (vce=="nn") {
      edups   <- dups[ind]
      edupsid <- dupsid[ind]
    #  for (j in 1:eN) {
    #    edups[j]=sum(eX==eX[j])
    #  }
    #  j=1
    #  while (j<=eN) {
    #    edupsid[j:(j+edups[j]-1)] <- 1:edups[j]
    #    j <- j+edups[j]
    #  }
    }
    
  u   <- (eX-eval[i])/h[i]
  R.q <- matrix(NA,eN,(q+1))
  for (j in 1:(q+1))  R.q[,j] <- (eX-eval[i])^(j-1)
  R.p <- R.q[,1:(p+1)]

  #display("Computing RD estimates.")
  L <- crossprod(R.p*W.h,u^(p+1)) 
  invG.q  <- qrXXinv((sqrt(W.b)*R.q))
  invG.p  <- qrXXinv((sqrt(W.h)*R.p))
  e.p1    <- matrix(0,(q+1),1); e.p1[p+2]=1
  e.v     <- matrix(0,(p+1),1); e.v[deriv+1]=1
  Q.q     <- t(t(R.p*W.h) - h[i]^(p+1)*(L%*%t(e.p1))%*%t(t(invG.q%*%t(R.q))*W.b))
  beta.p  <- invG.p%*%crossprod(R.p*W.h,eY); beta.q <- invG.q%*%crossprod(R.q*W.b,eY); beta.bc <- invG.p%*%crossprod(Q.q,eY) 

  tau.cl <- factorial(deriv)*beta.p[(deriv+1),1]
  tau.bc <- factorial(deriv)*beta.bc[(deriv+1),1]

  hii <- predicts.p <- predicts.q <- 0
  if (vce %in% c("hc0","hc1","hc2","hc3")) {
    predicts.p <- R.p %*% beta.p
    predicts.q <- R.q %*% beta.q
    if (vce %in% c("hc2","hc3") && is.null(eC)) {
      # vectorized diag(R %*% invG %*% t(R*W)) -- row-wise dot product
      hii <- matrix(rowSums((R.p %*% invG.p) * (R.p * W.h)), eN, 1)
    }
  }

  if (is.null(eC)) {
    res.h <- lprobust.res(eX, eY, predicts.p, hii, vce, nnmatch, edups, edupsid, p+1)
    res.b <- if (vce == "nn") res.h else lprobust.res(eX, eY, predicts.q, hii, vce, nnmatch, edups, edupsid, q+1)
    V.Y.cl <- invG.p %*% lprobust.vce(as.matrix(R.p * W.h), res.h, NULL) %*% invG.p
    V.Y.bc <- invG.p %*% lprobust.vce(Q.q, res.b, NULL) %*% invG.p
  } else {
    cr_type <- switch(vce,
                      "hc0" = "CR0",
                      "hc1" = "CR1",
                      "hc2" = "CR2",
                      "hc3" = "CR3",
                      "nn"  = "CR1",
                      "CR1")
    if (vce == "nn") {
      res.h.raw <- lprobust.res(eX, eY, predicts.p, hii, "nn", nnmatch, edups, edupsid, p+1)
      res.b.raw <- res.h.raw
    } else {
      res.h.raw <- matrix(as.numeric(eY) - as.numeric(predicts.p), ncol = 1)
      res.b.raw <- matrix(as.numeric(eY) - as.numeric(predicts.q), ncol = 1)
    }
    sqrtW.h  <- sqrt(W.h)
    X.std.h  <- R.p * sqrtW.h
    r.std.h  <- as.numeric(res.h.raw) * sqrtW.h
    # k_override = q+1 aligns the CR1 df correction with the q-regression
    # that produced res.b.raw. Without it, k=ncol(Q.q)=p+1 was being used.
    meat.cl  <- lprobust.cluster.meat(X.std.h, r.std.h, eC, invG.p, cr_type)
    meat.bc  <- lprobust.cluster.meat(Q.q, as.numeric(res.b.raw), eC, invG.p, cr_type, k_override = q + 1)
    V.Y.cl   <- invG.p %*% meat.cl %*% invG.p
    V.Y.bc   <- invG.p %*% meat.bc %*% invG.p
  }
	se.cl  <- sqrt(factorial(deriv)^2*V.Y.cl[deriv+1,deriv+1])
	se.rb  <- sqrt(factorial(deriv)^2*V.Y.bc[deriv+1,deriv+1])
	
	Estimate[i,] <- c(eval[i], h[i], b[i], eN, tau.cl, tau.bc, se.cl, se.rb) 
  }
  
  cov.us <- NULL
  cov.rb <- NULL

  if (isTRUE(covgrid)) {

    cov.us <- matrix(NA, neval, neval)
    cov.rb <- matrix(NA, neval, neval)

    ## Loop-invariants hoisted out of the (i, j) inner loop.
    e.p1     <- matrix(0, q + 1, 1); e.p1[p + 2]  <- 1
    e.v      <- matrix(0, p + 1, 1); e.v[deriv+1] <- 1
    fact_d   <- factorial(deriv)
    is_hc    <- vce %in% c("hc0", "hc1", "hc2", "hc3")
    is_h2or3 <- vce %in% c("hc2", "hc3")

    ## Diagonal: cov.us[i,i] = se.us[i]^2 (and similarly for cov.rb)
    ## by construction. Avoids neval redundant inner-loop iterations.
    diag(cov.us) <- as.numeric(Estimate[, "se.us"])^2
    diag(cov.rb) <- as.numeric(Estimate[, "se.rb"])^2

    if (neval > 1) for (i in 1:(neval - 1)) {
      for (j in (i + 1):neval) {

        w.h.i   <- W.fun((x-eval[i])/h[i], kernel)/h[i]
        w.b.i   <- W.fun((x-eval[i])/b[i], kernel)/b[i]
        ind.h.i <- w.h.i>0;  ind.b.i <- w.b.i>0
        ind.i   <- ind.h.i | ind.b.i

        w.h.j   <- W.fun((x-eval[j])/h[j], kernel)/h[j]
        w.b.j   <- W.fun((x-eval[j])/b[j], kernel)/b[j]
        ind.h.j <- w.h.j>0;  ind.b.j <- w.b.j>0
        ind.j   <- ind.h.j | ind.b.j

        ind <- ind.i | ind.j

        eN  <- sum(ind)
        eY  <- y[ind]
        eX  <- x[ind]

        W.h.i <- w.h.i[ind]
        W.b.i <- w.b.i[ind]
        W.h.j <- w.h.j[ind]
        W.b.j <- w.b.j[ind]

        eC <- if (!is.null(cluster)) cluster[ind] else NULL

        edups <- edupsid <- 0
        if (vce == "nn") {
          edups   <- dups[ind]
          edupsid <- dupsid[ind]
        }

        u.i   <- (eX-eval[i])/h[i]
        R.q.i <- outer(eX - eval[i], 0:q, `^`)
        R.p.i <- R.q.i[, 1:(p+1), drop = FALSE]

        u.j   <- (eX-eval[j])/h[j]
        R.q.j <- outer(eX - eval[j], 0:q, `^`)
        R.p.j <- R.q.j[, 1:(p+1), drop = FALSE]

        L.i       <- crossprod(R.p.i*W.h.i, u.i^(p+1))
        invG.q.i  <- qrXXinv((sqrt(W.b.i)*R.q.i))
        invG.p.i  <- qrXXinv((sqrt(W.h.i)*R.p.i))
        Q.q.i     <- t(t(R.p.i*W.h.i) - h[i]^(p+1)*(L.i%*%t(e.p1))%*%t(t(invG.q.i%*%t(R.q.i))*W.b.i))
        beta.p.i  <- invG.p.i%*%crossprod(R.p.i*W.h.i, eY)
        beta.q.i  <- invG.q.i%*%crossprod(R.q.i*W.b.i, eY)

        L.j       <- crossprod(R.p.j*W.h.j, u.j^(p+1))
        invG.q.j  <- qrXXinv((sqrt(W.b.j)*R.q.j))
        invG.p.j  <- qrXXinv((sqrt(W.h.j)*R.p.j))
        Q.q.j     <- t(t(R.p.j*W.h.j) - h[j]^(p+1)*(L.j%*%t(e.p1))%*%t(t(invG.q.j%*%t(R.q.j))*W.b.j))
        beta.p.j  <- invG.p.j%*%crossprod(R.p.j*W.h.j, eY)
        beta.q.j  <- invG.q.j%*%crossprod(R.q.j*W.b.j, eY)

        hii.i <- predicts.p.i <- predicts.q.i <- 0
        hii.j <- predicts.p.j <- predicts.q.j <- 0
        if (is_hc) {
          predicts.p.i <- R.p.i %*% beta.p.i
          predicts.q.i <- R.q.i %*% beta.q.i
          predicts.p.j <- R.p.j %*% beta.p.j
          predicts.q.j <- R.q.j %*% beta.q.j
          if (is_h2or3) {
            hii.i <- matrix(rowSums((R.p.i %*% invG.p.i) * (R.p.i * W.h.i)), eN, 1)
            hii.j <- matrix(rowSums((R.p.j %*% invG.p.j) * (R.p.j * W.h.j)), eN, 1)
          }
        }

        res.h.i <- lprobust.res(eX, eY, predicts.p.i, hii.i, vce, nnmatch, edups, edupsid, p+1)
        res.b.i <- if (vce == "nn") res.h.i
                   else lprobust.res(eX, eY, predicts.q.i, hii.i, vce, nnmatch, edups, edupsid, q+1)

        res.h.j <- lprobust.res(eX, eY, predicts.p.j, hii.j, vce, nnmatch, edups, edupsid, p+1)
        res.b.j <- if (vce == "nn") res.h.j
                   else lprobust.res(eX, eY, predicts.q.j, hii.j, vce, nnmatch, edups, edupsid, q+1)

        ## NOTE: V.us.i carries factorial(deriv) (not factorial(deriv)^2) so
        ## that cov(tau_i, tau_j) = (V.us.i %*% t(V.us.j))[d+1,d+1] picks up
        ## factorial(deriv)^2 from the cross-product -- matching the scaling
        ## of se.us = sqrt(factorial(deriv)^2 * V.Y.cl[d+1,d+1]) in the main
        ## loop. The previous implementation had factorial(deriv)^2 here,
        ## which gave factorial(deriv)^4 after the cross-product (a factor
        ## factorial(deriv)^2 too large for deriv >= 2).
        V.us.i <- fact_d * invG.p.i %*% t(c(res.h.i) * R.p.i * W.h.i)
        V.us.j <- fact_d * invG.p.j %*% t(c(res.h.j) * R.p.j * W.h.j)

        V.rb.i <- fact_d * invG.p.i %*% t(c(res.b.i) * Q.q.i)
        V.rb.j <- fact_d * invG.p.j %*% t(c(res.b.j) * Q.q.j)

        cov.us[i, j] <- (V.us.i %*% t(V.us.j))[deriv+1, deriv+1]
        cov.rb[i, j] <- (V.rb.i %*% t(V.rb.j))[deriv+1, deriv+1]

        cov.us[j, i] <- cov.us[i, j]
        cov.rb[j, i] <- cov.rb[i, j]
      }
    }

  }
    
    
  out        <-list(Estimate=Estimate, opt=list(p=p, q=q, deriv=deriv, kernel=kernel.type, n=N, neval=neval, bwselect=bwselect), cov.us=cov.us, cov.rb=cov.rb)
  out$call   <- match.call()
  class(out) <- "lprobust"
  return(out)
}

print.lprobust <- function(x,...){
  cat("Call: lprobust\n\n")

  cat(paste("Sample size (n)                              =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p)    =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated (deriv)        =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Polynomial order for confidence interval (q) =    ", x$opt$q,        "\n", sep=""))
  cat(paste("Kernel function                              =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                             =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  invisible(x)
}

summary.lprobust <- function(object, alpha = 0.05, sep = 5, ...) {
  z    <- qnorm(1 - alpha / 2)
  CI_l <- object$Estimate[, "tau.bc"] - object$Estimate[, "se.rb"] * z
  CI_r <- object$Estimate[, "tau.bc"] + object$Estimate[, "se.rb"] * z

  out <- list(opt      = object$opt,
              Estimate = object$Estimate,
              alpha    = alpha,
              sep      = sep,
              CI_l     = CI_l,
              CI_r     = CI_r)
  class(out) <- "summary.lprobust"
  out
}

print.summary.lprobust <- function(x, ...) {
  cat("Call: lprobust\n\n")

  cat(paste("Sample size (n)                              =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p)    =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated (deriv)        =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Polynomial order for confidence interval (q) =    ", x$opt$q,        "\n", sep=""))
  cat(paste("Kernel function                              =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                             =    ", x$opt$bwselect, "\n", sep=""))
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
  cat(format("h"               , width=10, justify="right"))
  cat(format("Eff.n"           , width=8 , justify="right"))
  cat(format("Est."            , width=10, justify="right"))
  cat(format("Error"           , width=10, justify="right"))
  cat(format(paste("[ ", floor((1-alpha)*100), "%", " C.I. ]", sep=""), width=25, justify="centre"))
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