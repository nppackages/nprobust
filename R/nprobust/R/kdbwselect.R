kdbwselect <- function(x, eval=NULL, neval=NULL, kernel="epa", 
                       bwselect="mse-dpi", bwcheck=21, imsegrid=30, subset = NULL){
  
  p <- 2
  deriv <- 0
  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)
  
  if (!is.null(subset)) x <- x[subset]
  na.ok <- complete.cases(x) 
  x <- x[na.ok]

  x.min <- min(x);  x.max <- max(x)
  N <- length(x)

  if (is.null(eval)) {
    if (is.null(neval)) {
      #eval <- unique(x)
      qseq <- seq(0.1,0.9,length.out=30)
      eval <- quantile(x, qseq)
      #eval <- seq(x.min, x.max, length.out=30+2)
      #eval <- eval[2:31]
    }
    else {
      qseq <- seq(0.1,0.9,length.out=neval)
      eval <- quantile(x, qseq)
      #eval <- seq(x.min,x.max,length.out=neval)
      #qseq <- seq(0,1,1/(neval+1))
      #eval <- quantile(x, qseq[2:(length(qseq)-1)])
      #eval <- seq(x.min, x.max, length.out=(neval+2))
      #eval <- eval[2:(neval+1)]
    }
  }
  neval <- length(eval)


  kernel.type <- "Gaussian"  
  C.h<-1.06
  C.b<-1
  if (kernel=="epa") {
    kernel.type <- "Epanechnikov"
    C.h<-2.34
    C.b<-3.49
  }
  if (kernel=="uni") {
    kernel.type <- "Uniform"
    C.h<-1.06
    C.b<-1
  }
  
  #C.fun <- function(v,r) {
  #  R.v <- kd.K.fun(1, v<-v, r<-r, kernel<-"gau")$R.v
  #  k.v <- kd.K.fun(1, v<-v, r<-0, kernel<-"gau")$k.v
  #  C <- 2*((sqrt(pi)*(1+2*d)*factorial(p+d)*R.v)/(2*p*k.v^2*factorial(2*p+2*d)))^(1/(2*p+2*d+1))
  #  return(C)
  #}

  
  bws           <- matrix(NA,neval,2)
  colnames(bws) <- c("h","b")
  bws.imse <- NULL  
  
  if (bwselect=="all") {
    bws           <- matrix(NA,neval,6)
    colnames(bws) <- c("h.mse.dpi","b.mse.dpi",  "h.ce.dpi","b.ce.dpi", "h.ce.rot","b.ce.rot")
    bws.imse      <- matrix(NA,2,2)
  }

  
  ### IMSE-ROT
  h.imse.rot <- sd(x)*C.h*N^(-1/(1+2*p))
  b.imse.rot <- sd(x)*C.b*N^(-1/(1+2*(p+2)+2*p))
  
  if (bwselect=="imse-rot") {
    bws[,1] <- rep(h.imse.rot, neval)
    bws[,2] <- rep(b.imse.rot, neval)
  }
  
  
  ### IMSE-ROT
  if (bwselect=="imse-dpi" | bwselect=="all") {
  
  B.h=V.h=0  
  #qseq.imse <- seq(0,1,1/(imsegrid+1))
  #eval.imse <- quantile(x, qseq[2:(length(qseq)-1)])
  qseq.imse <- seq(0.1,0.9,length.out=imsegrid)
  eval.imse <- quantile(x, qseq.imse)
  
  #eval.imse <- seq(x.min, x.max, length.out=(imsegrid+2))
  #eval.imse <- eval[2:(imsegrid+1)]
  
  for (i in 1:imsegrid) {
    K.b <- kd.K.fun((x-eval.imse[i])/b.imse.rot, v=p+2, r=p,     kernel=kernel)    
    K.h <- kd.K.fun((x-eval.imse[i])/h.imse.rot, v=p,   r=deriv, kernel=kernel)
    f.b <- mean(K.b$Kx)/b.imse.rot^(1+p)
    f.h.rot <- mean(K.h$Kx)/h.imse.rot    
    B.h[i]  <- f.b*K.h$k.v
    V.h[i]  <- f.h.rot*K.h$R.v
  }
  
  h.imse.dpi <- kd.bw.fun(mean(V.h), mean(B.h), N, v=p, r=deriv)
  
  if (bwselect=="imse-dpi") {
    bws[,1] <- rep(h.imse.dpi, neval)
    bws[,2] <- rep(b.imse.rot, neval)
  }

  }
  
  if  (bwselect=="all") {
    bws.imse[,1] <- c(h.imse.dpi,  b.imse.rot)
    bws.imse[,2] <- c(h.imse.rot,  b.imse.rot)
  }
  
  
  if  (bwselect=="all"  | bwselect=="mse-dpi" | bwselect=="ce-dpi" | bwselect=="ce-rot" ) {
  
    B.h=V.h=0 
    
  for (i in 1:neval) {
  
    #bws[i,1:2] <- c(h.imse.rot, b.imse.rot)
    #q.rot <- sd(x)*C.b*N^(-1/(1+2*(p+4)+2*(p+2)))
    
    #V = (mean((K^2)) - mean(K)^2)/(N*h.rot^2)
    #B = mean(rho^(1+p)*M)/h.rot
    
    #K.q = kd.K.fun((x-eval[i])/q.rot, v=p+4, r=p+2,   kernel=kernel)    
    #f.q.rot = mean(K.q$Kx)/q.rot^(1+2*p)
    #f.q.rot <- kd.K.fun(eval[i], v=2, r=4, kernel="gau")$Kx
    
    #V.b <- f.h.rot*kd.K.fun(1, v=p+2, r=p, kernel=kernel)$R.v
    #B.b <- f.q.rot*kd.K.fun(1, v=p+2, r=p, kernel=kernel)$k.v
    #b.mse.dpi <- bw.fun(V.b, B.b, N, v=p+2, r=p)
    
    if (!is.null(bwcheck)) {
      bw.min     <- sort(abs(x-eval[i]))[bwcheck]
      #h.imse.rot <- max(h.imse.rot, bw.min)
      #b.imse.rot <- max(b.imse.rot, bw.min)
    }
    
    K.b <- kd.K.fun((x-eval[i])/b.imse.rot, v=p+2, r=p,     kernel=kernel)    
    K.h <- kd.K.fun((x-eval[i])/h.imse.rot, v=p,   r=deriv, kernel=kernel)
    f.b <- mean(K.b$Kx)/b.imse.rot^(1+p)
    f.h.rot <- mean(K.h$Kx)/h.imse.rot    
    B.h[i]  <- f.b*K.h$k.v
    V.h[i]  <- f.h.rot*K.h$R.v
    h.mse.dpi <- kd.bw.fun(V.h[i], B.h[i], N, v=p, r=deriv)
    
    #f.h.rot=dnorm(eval[i])
    #f.b.rot=(eval[i]^2-1)*dnorm(eval[i])
    #V.h = (f.h.rot*K.h$v.k)/(N*h.rot^(1+2*deriv))
    #B.h = (f.b.rot*K.h$m.k)*b.rot^p
    #u = (x-eval[i])/h.rot
    #K.d = kd.K.fun(u     , r=p,   d=deriv, kernel=kernel)
    #L.r = kd.K.fun(rho[i]*u , r=p+2, d=p,     kernel=kernel)
    #K = K.d$k.x
    #M = rho[i]^(1+p)*L.r$k.x*L.r$m.k
    #f.hat.us = mean(K)/h.rot
    #bias = mean(M)/h.rot
    #var.us = (mean((K^2)) - mean(K)^2)/(N*h.rot^2)
    #bw.fun(var.us, f.hat.bc, r=p, d=deriv)
    #(var.us/(2*p*N*bias^2))^1/(1+2*p)
    b.mse.dpi = b.imse.rot
    
    if (!is.null(bwcheck)) {
      h.mse.dpi <- max(h.mse.dpi, bw.min)
      b.mse.dpi <- max(b.mse.dpi, bw.min)
    }
    
    bws[i,1:2] <- c(h.mse.dpi,    b.mse.dpi)
    
    if (bwselect=="ce-rot" | bwselect=="all") {
      h.ce.rot <- h.mse.dpi*N^(-(p-2)/((1+2*p)*(1+p+2)))
      b.ce.rot <- b.mse.dpi*N^(-(p-2)/((1+2*p)*(1+p+2)))
      
      if (!is.null(bwcheck)) {
        h.ce.rot <- max(h.ce.rot, bw.min)
        b.ce.rot <- max(b.ce.rot, bw.min)
      }
      
      bws[i,1:2] <- c(h.ce.rot, b.ce.rot)
    }
    
    if (bwselect=="ce-dpi"| bwselect=="all") {
      h.ce.dpi <- kd.cer.fun(x, eval[i], h.mse.dpi, h.mse.dpi, p, kernel)
      b.ce.dpi <- b.mse.dpi
      
      if (!is.null(bwcheck)) {
        h.ce.dpi <- max(h.ce.dpi, bw.min)
        b.ce.dpi <- max(b.ce.dpi, bw.min)
      }
      
      bws[i,1:2] <- c(h.ce.dpi,  b.ce.dpi)
    }
    
    if(bwselect=="all") bws[i,1:6] <- c(h.mse.dpi,b.mse.dpi, h.ce.dpi,b.ce.dpi, h.ce.rot,b.ce.rot)
  }
  }
  
  
  
  bws <- cbind(eval,bws)
  out <- list(bws=bws, bws.imse = bws.imse, opt=list(p=p, n=N, neval=neval, kernel=kernel.type, bwselect=bwselect))
  out$call   <- match.call()
  class(out) <- "kdbwselect"
  return(out)
}

print.kdbwselect <- function(x,...){
  cat("Call: kdbwselect\n\n")
  
  cat(paste("Sample size (n)                         =    ", x$opt$n,        "\n", sep=""))
#  cat(paste("Kernel order for point estimation (p)   =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Kernel function                         =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                        =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")
  
  #cat("Use summary(...) to show bandwidths.\n")
}

summary.kdbwselect <- function(object,...) {
  x <- object
  args <- list(...)
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }
  
  cat("Call: kdbwselect\n\n")
  
  cat(paste("Sample size (n)                             =    ", x$opt$n,        "\n", sep=""))
#  cat(paste("Polynomial order for point estimation (p)   =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Kernel function                             =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth selection method                  =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")

  if (x$opt$bwselect=="all") {
    col1.names = c("","", "MSE-DPI", "","CE-DPI","","CE-ROT")
    col2.names = rep(c("h", "b"),3)
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
    cat(format(sprintf("%3.3f", x$bws[1,2:3])  , width=9, justify="right"))
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