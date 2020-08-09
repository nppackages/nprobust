kdrobust <- function(x, eval=NULL, neval=NULL, h=NULL, b=NULL, rho=1, kernel="epa", 
                    bwselect=NULL, bwcheck=21, imsegrid=30, level=95, subset = NULL) {
  
  p <- 2
  deriv <- 0
  if (!is.null(subset)) x <- x[subset]
  na.ok <- complete.cases(x) 
  x <- x[na.ok]
  
  x.min <- min(x);  x.max <- max(x)
  N <- length(x)
  #quant <- -qnorm(abs((1-(level/100))/2))
  
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
  
  if (is.null(h) & is.null(bwselect) & neval==1) bwselect="mse-dpi"
  if (is.null(h) & is.null(bwselect) & neval>1)  bwselect="imse-dpi"  
  
  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)
  
  #####################################################   CHECK ERRORS
  exit<-0
    #if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    #  print("kernel incorrectly specified")
    #  exit = 1
    #}
    
  #if  (bwselect!="imse-rot" & bwselect!="imse-dpi" & bwselect!="mse-dpi" & bwselect!="ce-dpi" & bwselect!="ce-rot" & bwselect!="all" & !is.null(bwselect)){
  #  print("bwselect incorrectly specified")  
  #  exit = 1
  #}
  
    #if (min(eval)<x.min | max(eval)>x.max){
    #  print("evaluation points should be set within the range of x")
    #  exit = 1
    #}
    
    if (p<0 | deriv<0 ){
      print("p should be positive integer")
      exit = 1
    }
    
    
    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit = 1
    }
    
    #if (!is.null(rho)){  
       if (rho<0){
          print("rho should be greater than 0")
          exit = 1
        }
    #}
  
    if (exit>0) stop()
    if (!is.null(h)) bwselect = "Manual"

  kernel.type <- "Gaussian"  
  if (kernel=="epa") kernel.type <- "Epanechnikov"
  if (kernel=="uni") kernel.type <- "Uniform"

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
    
    eN = sum(M>0)
    
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
  
# cat("Use summary(...) to show estimates.\n")
}

summary.kdrobust <- function(object,...) {
  x <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }
  
  cat("Call: kdrobust\n\n")
  
  cat(paste("Sample size (n)                            =     ", x$opt$n,        "\n", sep=""))
  cat(paste("Kernel order for point estimation (p)      =     ", x$opt$p,        "\n", sep=""))
  cat(paste("Kernel function                            =     ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth selection method                 =     ", x$opt$bwselect, "\n", sep=""))
  cat("\n")
  
  ### compute CI
  z <- qnorm(1 - alpha / 2)
  CI_l <- x$Estimate[, "tau.bc"] - x$Estimate[, "se.rb"] * z;
  CI_r <- x$Estimate[, "tau.bc"] + x$Estimate[, "se.rb"] * z;
  
  ### print output
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
    cat(format(paste("[", sprintf("%3.3f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_r[j]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
      cat(paste(rep("-", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
    }
  }
  
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
}