lprobust = function(y, x, eval=NULL, neval=NULL, p=NULL, deriv=NULL, h=NULL, b=NULL, rho=1, 
                    kernel="epa", bwselect=NULL, bwcheck=21, bwregul=1, imsegrid=30, vce="nn", covgrid = FALSE,
                    cluster=NULL, nnmatch=3, level=95, interior = FALSE, subset = NULL) {
  
  if (!is.null(subset)) { 
    x <- x[subset]
    y <- y[subset]
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)

  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  x <- x[na.ok]
  y <- y[na.ok]
  if (!is.null(cluster)) cluster = cluster[na.ok]
  
  if (!is.null(deriv) & is.null(p)) p <- deriv+1
  if (is.null(p))         p <- 1
  if (is.null(deriv)) deriv <- 0
  q <- p+1
  
  x.max <- max(x); x.min <- min(x)
  N <- length(x)
  
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
  
  if (is.null(h) & is.null(bwselect) & neval==1) bwselect="mse-dpi"
  if (is.null(h) & is.null(bwselect) & neval>1)  bwselect="imse-dpi"  
  
  if (vce=="nn") {
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
    if (!is.null(cluster)) cluster = cluster[order.x]
  }

  kernel   <- tolower(kernel)
  bwselect <- tolower(bwselect)
  vce      <- tolower(vce)

  vce_type = "NN"
  if (vce=="hc0")     		vce_type = "HC0"
  if (vce=="hc1")      	  vce_type = "HC1"
  if (vce=="hc2")      	  vce_type = "HC2"
  if (vce=="hc3")      	  vce_type = "HC3"
  if (vce=="cluster")  	  vce_type = "Cluster"
  if (vce=="nncluster") 	vce_type = "NNcluster"
  
  #####################################################   CHECK ERRORS
  exit=0
    if (kernel!="gau" & kernel!="gaussian" & kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit <- 1
    }
    
  #if  (bwselect!="mse-dpi" & bwselect!="mse-rot" & bwselect!="imse-dpi" & bwselect!="imse-rot" & bwselect!="ce-dpi" & bwselect!="ce-rot" & bwselect!=NULL){
  #  print("bwselect incorrectly specified")  
  #  exit <- 1
  #}
  
  if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
    print("vce incorrectly specified")
    exit <- 1
  }
    
   # if (min(eval)<x.min | max(eval)>x.max){
  #    print("eval should be set within the range of x")
  #    exit <- 1
  #  }
    
  if (p<0 | deriv<0 | nnmatch<=0 ){
     print("p,q,deriv and matches should be positive integers")
     exit <- 1
  }
    
  if (deriv>p){
    print("deriv can only be equal or lower p")
    exit <- 1
  }
    
    #p.round <- round(p)/p;  q.round <- round(q)/q;  d.round <- round(deriv+1)/(deriv+1);  m.round <- round(nnmatch)/nnmatch
    
    #if (p.round!=1 | q.round!=1 | d.round!=1 | m.round!=1 ){
    #  print("p,q,deriv and matches should be integer numbers")
    #  exit <- 1
    #}
    
    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit <- 1
    }
    
    #if (!is.null(rho)){  
      if (rho<0){
          print("rho should be greater than 0")
          exit <- 1
      }
    #}
  
    if (exit>0) stop()
    if (!is.null(h)) bwselect <- "Manual"
    #if (is.null(rho)) rho <- 1
    if (!is.null(h) & rho>0 & is.null(b)) {
    #  rho <- 1
        b <- h/rho
    }
    #if (!is.null(h) & rho>0) b <- h/rho
    
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

    
  ############################################################################################
    #print("Preparing data.") 
    if (is.null(h)) {
        lpbws <- lpbwselect(y=y, x=x,  eval=eval, deriv=deriv, p=p, vce=vce, cluster=cluster, bwselect=bwselect, interior=interior, kernel=kernel, bwcheck = bwcheck, bwregul=bwregul, imsegrid=imsegrid, subset=subset)
        h     <- lpbws$bws[,2]
        b     <- lpbws$bws[,3]
        if (rho>0) b <- h/rho
        rho   <- h/b
    }
  
    if (length(h)==1 & neval>1) {
      h <- rep(h,neval)
      b <- rep(b,neval)
      rho   <- h/b
    }
  
  dups <- dupsid <- 0	
  if (vce=="nn") {
    for (j in 1:N) {
      dups[j]=sum(x==x[j])
    }
    j=1
    while (j<=N) {
      dupsid[j:(j+dups[j]-1)] <- 1:dups[j]
      j <- j+dups[j]
    }
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
    
    w.h   <- W.fun((x-eval[i])/h[i], kernel)/h[i]
    w.b   <- W.fun((x-eval[i])/b[i], kernel)/b[i]
    ind.h <- w.h>0;  ind.b <- w.b>0
    N.h   <- sum(ind.h);  N.b <- sum(ind.b)
    ind   <- ind.b
    if (h[i]>b[i]) ind <- ind.h   
    
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

  hii=predicts.p=predicts.q=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts.p <- R.p%*%beta.p
    predicts.q <- R.q%*%beta.q
    if (vce=="hc2" | vce=="hc3") {
      hii <- matrix(NA,eN,1)	
      for (j in 1:eN) hii[j] <- R.p[j,]%*%invG.p%*%(R.p*W.h)[j,]
    }
  }
  						
	res.h <- lprobust.res(eX, eY, predicts.p, hii, vce, nnmatch, edups, edupsid, p+1)
	if (vce=="nn") res.b <- res.h
	else           res.b <- lprobust.res(eX, eY, predicts.q, hii, vce, nnmatch, edups, edupsid, q+1)
  
	V.Y.cl <- invG.p%*%lprobust.vce(as.matrix(R.p*W.h), res.h, eC)%*%invG.p
	V.Y.bc <- invG.p%*%lprobust.vce(Q.q,     res.b, eC)%*%invG.p
	se.cl  <- sqrt(factorial(deriv)^2*V.Y.cl[deriv+1,deriv+1])
	se.rb  <- sqrt(factorial(deriv)^2*V.Y.bc[deriv+1,deriv+1])
	
	Estimate[i,] <- c(eval[i], h[i], b[i], eN, tau.cl, tau.bc, se.cl, se.rb) 
  }
  
  cov.us = cov.rb = matrix(NA,neval,neval)
  
  if (covgrid == TRUE) {
  
    for (i in 1:neval) {
      for (j in i:neval) {
    
        #h.max     <- max(h[i], h[j])
        #b.max     <- max(b[i], b[j])
        
        w.h.i   <- W.fun((x-eval[i])/h[i], kernel)/h[i]
        w.b.i   <- W.fun((x-eval[i])/b[i], kernel)/b[i]
        ind.h.i <- w.h.i>0;  ind.b.i <- w.b.i>0
        N.h.i   <- sum(ind.h.i);  N.b.i <- sum(ind.b.i)
        ind.i   <- ind.b.i
        if (h[i]>b[i]) ind.i <- ind.h.i   
        
        w.h.j   <- W.fun((x-eval[j])/h[j], kernel)/h[j]
        w.b.j   <- W.fun((x-eval[j])/b[j], kernel)/b[j]
        ind.h.j <- w.h.j>0;  ind.b.j <- w.b.j>0
        N.h.j   <- sum(ind.h.j);  N.b.j <- sum(ind.b.j)
        ind.j   <- ind.b.j
        if (h[j]>b[j]) ind.j <- ind.h.j   
        
        ind = ind.i=="TRUE" | ind.j=="TRUE"
        
        eN  <- sum(ind)
        eY  <- y[ind]
        eX  <- x[ind]
        
        W.h.i <- w.h.i[ind]
        W.b.i <- w.b.i[ind]
        
        W.h.j <- w.h.j[ind]
        W.b.j <- w.b.j[ind]
        
        eC = NULL
        if (!is.null(cluster)) eC = cluster[ind]
        
        edups <- edupsid <- 0	
        if (vce=="nn") {
          edups   <- dups[ind]
          edupsid <- dupsid[ind]
        #  for (k in 1:eN) {
        #    edups[k]=sum(eX==eX[k])
        #  }
        #  k=1
        #  while (k<=eN) {
        #    edupsid[k:(k+edups[k]-1)] <- 1:edups[k]
        #    k <- k+edups[k]
        #  }
        }
        
        u.i   <- (eX-eval[i])/h[i]
        R.q.i <- matrix(NA,eN,(q+1))
        for (k in 1:(q+1))  R.q.i[,k] <- (eX-eval[i])^(k-1)
        R.p.i <- R.q.i[,1:(p+1)]
        
        u.j   <- (eX-eval[j])/h[j]
        R.q.j <- matrix(NA,eN,(q+1))
        for (k in 1:(q+1))  R.q.j[,k] <- (eX-eval[j])^(k-1)
        R.p.j <- R.q.j[,1:(p+1)]
        
        e.p1    <- matrix(0,(q+1),1); e.p1[p+2]=1
        e.v     <- matrix(0,(p+1),1); e.v[deriv+1]=1
        
        #display("Computing RD estimates.")
        L.i <- crossprod(R.p.i*W.h.i,u.i^(p+1)) 
        invG.q.i  <- qrXXinv((sqrt(W.b.i)*R.q.i))
        invG.p.i  <- qrXXinv((sqrt(W.h.i)*R.p.i))
        Q.q.i     <- t(t(R.p.i*W.h.i) - h[i]^(p+1)*(L.i%*%t(e.p1))%*%t(t(invG.q.i%*%t(R.q.i))*W.b.i))
        beta.p.i  <- invG.p.i%*%crossprod(R.p.i*W.h.i,eY); beta.q.i <- invG.q.i%*%crossprod(R.q.i*W.b.i,eY); 
        beta.bc.i <- invG.p.i%*%crossprod(Q.q.i,eY) 

        L.j <- crossprod(R.p.j*W.h.j,u.j^(p+1)) 
        invG.q.j  <- qrXXinv((sqrt(W.b.j)*R.q.j))
        invG.p.j  <- qrXXinv((sqrt(W.h.j)*R.p.j))
        Q.q.j     <- t(t(R.p.j*W.h.j) - h[j]^(p+1)*(L.j%*%t(e.p1))%*%t(t(invG.q.j%*%t(R.q.j))*W.b.j))
        beta.p.j  <- invG.p.j%*%crossprod(R.p.j*W.h.j,eY); beta.q.j <- invG.q.j%*%crossprod(R.q.j*W.b.j,eY); 
        beta.bc.j <- invG.p.j%*%crossprod(Q.q.j,eY) 
        
        #tau.cl <- factorial(deriv)*beta.p[(deriv+1),1]
        #tau.bc <- factorial(deriv)*beta.bc[(deriv+1),1]
        
        hii.i=predicts.p.i=predicts.q.i=hii.j=predicts.p.j=predicts.q.j=0
        if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
          predicts.p.i <- R.p.i%*%beta.p.i
          predicts.q.i <- R.q.i%*%beta.q.i
          predicts.p.j <- R.p.j%*%beta.p.j
          predicts.q.j <- R.q.j%*%beta.q.j
          if (vce=="hc2" | vce=="hc3") {
            hii.i = hii.j = matrix(NA,eN,1)	
            for (k in 1:eN) {
              hii.i[k] <- R.p.i[k,]%*%invG.p.i%*%(R.p.i*W.h.i)[k,]
              hii.j[k] <- R.p.j[k,]%*%invG.p.j%*%(R.p.j*W.h.j)[k,]
            }
          }
        }
        
        res.h.i <- lprobust.res(eX, eY, predicts.p.i, hii.i, vce, nnmatch, edups, edupsid, p+1)
        if (vce=="nn") {
          res.b.i <- res.h.i
        }      else           res.b.i <- lprobust.res(eX, eY, predicts.q.i, hii.i, vce, nnmatch, edups, edupsid, q+1)
        
        
        res.h.j <- lprobust.res(eX, eY, predicts.p.j, hii.j, vce, nnmatch, edups, edupsid, p+1)
        if (vce=="nn") {
          res.b.j <- res.h.j
        }       else           res.b.j <- lprobust.res(eX, eY, predicts.q.j, hii.j, vce, nnmatch, edups, edupsid, q+1)
        
        #V.Y.cl <- invG.p%*%lprobust.vce(as.matrix(R.p*W.h), res.h, eC)%*%invG.p
        #V.Y.bc <- invG.p%*%lprobust.vce(Q.q,     res.b, eC)%*%invG.p
        
        V.us.i =  factorial(deriv)^2*invG.p.i%*%t(c(res.h.i)*R.p.i*W.h.i)
        V.us.j =  factorial(deriv)^2*invG.p.j%*%t(c(res.h.j)*R.p.j*W.h.j)
        
        V.rb.i =  factorial(deriv)^2*invG.p.i%*%t(c(res.b.i)*Q.q.i)
        V.rb.j =  factorial(deriv)^2*invG.p.j%*%t(c(res.b.j)*Q.q.j)
        
        cov.us[i,j] = (V.us.i%*%t(V.us.j))[deriv+1,deriv+1]
        cov.rb[i,j] = (V.rb.i%*%t(V.rb.j))[deriv+1,deriv+1]
        
        cov.us[j,i]= cov.us[i,j]
        cov.rb[j,i]= cov.rb[i,j]
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
  
 # cat("Use summary(...) to show estimates.\n")
  
}

summary.lprobust <- function(object,...) {
  x    <- object
  args <- list(...)
  if (is.null(args[['alpha']])) { alpha <- 0.05 } else { alpha <- args[['alpha']] }
  if (is.null(args[['sep']]))   { sep <- 5 } else { sep <- args[['sep']] }
  
  cat("Call: lprobust\n\n")
  
  cat(paste("Sample size (n)                              =    ", x$opt$n,        "\n", sep=""))
  cat(paste("Polynomial order for point estimation (p)    =    ", x$opt$p,        "\n", sep=""))
  cat(paste("Order of derivative estimated (deriv)        =    ", x$opt$deriv,    "\n", sep=""))
  cat(paste("Polynomial order for confidence interval (q) =    ", x$opt$q,        "\n", sep=""))
  cat(paste("Kernel function                              =    ", x$opt$kernel,   "\n", sep=""))
  cat(paste("Bandwidth method                             =    ", x$opt$bwselect, "\n", sep=""))
  cat("\n")
  
  ### compute CI
  z    <- qnorm(1 - alpha / 2)
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
    cat(format(paste("[", sprintf("%3.3f", CI_l[j]), " , ", sep="")  , width=14, justify="right"))
    cat(format(paste(sprintf("%3.3f", CI_r[j]), "]", sep=""), width=11, justify="left"))
    cat("\n")
    if (is.numeric(sep)) if (sep > 0) if (j %% sep == 0) {
      cat(paste(rep("-", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
    }
  }
  
  cat(paste(rep("=", 14 + 10 + 8 + 10 + 10 + 25), collapse="")); cat("\n")
}