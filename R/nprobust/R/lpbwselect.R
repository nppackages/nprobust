lpbwselect = function(y, x, eval=NULL, neval=NULL, p=NULL, deriv=NULL, kernel="epa", 
                      bwselect="mse-dpi", bwcheck=21, bwregul=1, imsegrid=30, vce="nn", cluster = NULL, 
                      nnmatch=3, interior=FALSE, subset=NULL){
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  x     <- x[na.ok]
  y     <- y[na.ok]
  if (!is.null(cluster)) cluster = cluster[na.ok]

  
  x.max <- max(x); x.min <- min(x)
  N <- length(x)
  
  if (!is.null(deriv) & is.null(p)) p <- deriv+1
  if (is.null(p))         p <- 1
  if (is.null(deriv)) deriv <- 0
  q <- p+1
  
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

  if  (bwselect=="imse-dpi" | bwselect=="imse-rot") neval =  eval = 1
  
  even <- (p-deriv)%%2==0
  
  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)
  vce      <- tolower(vce)
  
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"
    
  bws           <- matrix(NA,neval,2)
  colnames(bws) <- c("h","b")
  bws.imse      <- NULL  
  
  if (bwselect=="all") {
    bws           <- matrix(NA,neval,8)
    colnames(bws) <- c("h.mse.dpi","b.mse.dpi", "h.mse.rot","b.mse.rot", "h.ce.dpi","b.ce.dpi", "h.ce.rot","b.ce.rot")
    bws.imse      <- matrix(NA,2,2)
  }

  if  (bwselect=="imse-dpi" | bwselect=="all") {
      est <- lpbwselect.imse.dpi(y=y, x=x, cluster=cluster,               p=p, q=q, deriv=deriv, kernel=kernel, bwcheck=bwcheck, bwregul=bwregul, imsegrid=imsegrid, vce=vce, nnmatch=nnmatch, interior=interior)
      h.imse.dpi   <- est$h
      b.imse.dpi   <- est$b
      bws[1,1:2]   <- c(h.imse.dpi,  b.imse.dpi)
  }
    
  if  (bwselect=="imse-rot" | bwselect=="all") {
      est <- lpbwselect.imse.rot(y=y, x=x, p=p, deriv=deriv, kernel=kernel, imsegrid=imsegrid)
      h.imse.rot   <- est$h
      est <- lpbwselect.imse.rot(y=y, x=x, p=q, deriv=p+1,   kernel=kernel, imsegrid=imsegrid)
      b.imse.rot   <- est$h
      bws[1,1:2]   <- c(h.imse.rot,  b.imse.rot)
  }
    
  if  (bwselect=="all") {
    bws.imse[,1] <- c(h.imse.dpi,  b.imse.dpi)
    bws.imse[,2] <- c(h.imse.rot,  b.imse.rot)
  }
  
  if  (bwselect=="all"  | bwselect=="mse-dpi" | bwselect=="mse-rot" | bwselect=="ce-dpi" | bwselect=="ce-rot") {
   
     for (i in 1:neval) {
     
      if  (bwselect=="mse-dpi" | bwselect=="ce-dpi" | bwselect=="ce-rot" | bwselect=="all") {
        est <- lpbwselect.mse.dpi(y=y, x=x, cluster=cluster, eval=eval[i], p=p, q=q, deriv=deriv, kernel=kernel, 
                                bwcheck=bwcheck, bwregul=bwregul, vce=vce, nnmatch=nnmatch, interior=interior)
        h.mse.dpi  <- est$h
        b.mse.dpi  <- est$b
        bws[i,1:2] <- c(h.mse.dpi,  b.mse.dpi)
      }

      if  (bwselect=="mse-rot" | bwselect=="all") {
        est <- lpbwselect.mse.rot(y=y, x=x, eval=eval[i], p=p, deriv=deriv, kernel=kernel)
        h.mse.rot  <- est$h
        est <- lpbwselect.mse.rot(y=y, x=x, eval=eval[i], p=q, deriv=p+1,   kernel=kernel)
        b.mse.rot  <- est$h
        bws[i,1:2] <- c(h.mse.rot,  b.mse.rot)
      }
      
      h.ce.dpi=b.ce.dpi=h.ce.rot=b.ce.rot=0
      
      if  (bwselect=="ce-dpi" | bwselect=="all") {
        h.ce.dpi=b.ce.dpi=0
        #if (deriv==0) {
          if  (even==TRUE ) { 
            h.ce.dpi <- h.mse.dpi*N^(-((p+2)/((2*p+5)*(p+3))))
            b.ce.dpi <- b.mse.dpi*N^(-((q)/((2*q+3)*(q+3))))
          } else{
            est <- lpbwselect.ce.dpi(y=y, x=x, h=h.mse.dpi, b=b.mse.dpi, eval=eval[i], p=p, q=q, deriv=deriv, rho=1, 
                                   kernel=kernel, vce=vce, nnmatch=nnmatch, interior=interior, bwregul=bwregul)
            h.ce.dpi <- est$h
            b.ce.dpi <- b.mse.dpi*N^(-((q+2)/((2*q+5)*(q+3))))
          }
            bws[i,1:2] <- c(h.ce.dpi,  b.ce.dpi)
        #}
      }
      
      if  (bwselect=="ce-rot" | bwselect=="all") {
        h.ce.rot=b.ce.rot=0
        #if (deriv==0) {
          if  (even==TRUE ) { 
            h.ce.rot <- h.mse.dpi*N^(-((p+2)/((2*p+5)*(p+3))))
            b.ce.rot <- b.mse.dpi*N^(-((q)/((2*q+3)*(q+3))))
          } else{
            h.ce.rot <- h.mse.dpi*N^(-((p)/((2*p+3)*(p+3))))
            b.ce.rot <- b.mse.dpi*N^(-((q+2)/((2*q+5)*(q+3))))
          }
          bws[i,1:2]  <- c(h.ce.rot,  b.ce.rot)
      #}
      }
    
      if (bwselect=="all") bws[i,] <- c(h.mse.dpi,b.mse.dpi, h.mse.rot,b.mse.rot,  h.ce.dpi,b.ce.dpi, h.ce.rot,b.ce.rot)
    }
  }
  
  bws <- cbind(eval, bws)
  out        <- list(bws = bws, bws.imse = bws.imse,
                     opt = list(n=N, neval=neval, p=p, q=q, deriv=deriv, kernel=kernel.type, bwselect=bwselect))
  out$call   <- match.call()
  class(out) <- "lpbwselect"
  return(out)
}

print.lpbwselect <- function(x,...){
  cat("Call: lpbwselect\n\n")

  cat(paste("Sample size (n)                              =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p)    =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated (deriv)        =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Polynomial order for confidence interval (q) =    ", x$opt$q,        "\n", sep=""))
  cat(paste("Kernel function                              =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                             =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")
  #cat("Use summary(...) to show bandwidths.\n")
}

summary.lpbwselect <- function(object,...) {
  x <- object
  args <- list(...)
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }
  
  cat("Call: lpbwselect\n\n")
  
  cat(paste("Sample size (n)                              =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p)    =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated (deriv)        =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Polynomial order for confidence interval (q) =    ", x$opt$q,        "\n", sep=""))
  cat(paste("Kernel function                              =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                             =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")
  
  if (x$opt$bwselect=="all") {
    col1.names = c("","", "MSE-DPI","", "MSE-ROT", "","CE-DPI", "","CE-ROT")
    col2.names = rep(c("h", "b"),4)
  } else {
    col1.names = c("")
    col2.names = c("h", "b")
  }

  ### print output
  if (x$opt$bwselect=="imse-dpi" | x$opt$bwselect=="imse-rot") {
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 15 + 8*ncol(x$bws)), collapse="")); cat("\n")
  }
  if (x$opt$bwselect=="all") {
    cat(format(col1.names  , width=8, justify="right"))
    cat("\n")
  }
  if (x$opt$bwselect!="imse-dpi" & x$opt$bwselect!="imse-rot") cat(format("eval", width=10, justify="right"))
  cat(format(col2.names            , width=8, justify="right"))
  cat("\n")
  
  if (x$opt$bwselect=="imse-dpi" | x$opt$bwselect=="imse-rot") {
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 15 + 8*ncol(x$bws)), collapse="")); cat("\n")
  }
  
  if (x$opt$bwselect=="imse-dpi" | x$opt$bwselect=="imse-rot") {
    cat(format(sprintf("%3.3f", x$bws[2:3])  , width=9, justify="right"))
    cat("\n")
  } else {
    for (j in 1:nrow(x$bws)) {
      cat(format(toString(j), width=4))
      cat(format(sprintf("%3.3f", x$bws[j, "eval"]), width=8, justify="right"))
      cat(format(sprintf("%3.3f", x$bws[j, 2:ncol(x$bws)])  , width=8, justify="right"))
      cat("\n")
      if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
        cat(paste(rep("-", 15 + 8*ncol(x$bws)), collapse="")); cat("\n")
      }
    }
  }
  
  if (x$opt$bwselect=="imse-dpi" | x$opt$bwselect=="imse-rot") {
    cat(paste(rep("=", 15 + 8), collapse="")); cat("\n")
  } else {
    cat(paste(rep("=", 15 + 8*ncol(x$bws)), collapse="")); cat("\n")   
  }
  
  if (x$opt$bwselect=="all") {
    cat("\n")
    cat(paste(rep("=", 15 + 10 + 10), collapse="")); cat("\n")
    cat(format(c("","IMSE-DPI","", "IMSE-ROT") , width=8, justify="right"))
    cat("\n")
    cat(format(c("h", "b", "h", "b")            , width=8, justify="right"))
    cat("\n")
    cat(paste(rep("=", 15 + 10 + 10), collapse="")); cat("\n")
    cat(format(sprintf("%3.3f", x$bws.imse)  , width=8, justify="right"))
    cat("\n")
    cat(paste(rep("=", 15 + 10 + 10), collapse="")); cat("\n")
    cat("\n")
  }
  
  
}
