W.fun = function(u,kernel){
  if (kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w =          0.5*(abs(u)<=1)
  if (kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w =   dnorm(u)
  return(w)
}

kd.bw.fun = function(V, B, N, v, r) ((1+2*r)*V/(2*v*N*B^2))^(1/(1+2*v+2*r))

kd.K.fun = function(x,v,r,kernel){
  if (v==2) {
    if (kernel=="gau"){
      if (r==0) k = function(u)            dnorm(u)
      if (r==2) k = function(u)     (u^2-1)*dnorm(u)
      if (r==4) k = function(u) (u^4-6*u^2+3)*dnorm(u)
    }
    if (kernel=="uni"){
      if (r==0) k = function(u)  0.5*(abs(u)<=1)
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)  0.75*(1-u^2)*(abs(u)<=1)
    }
  }
  if (v==4) {
    if (kernel=="uni"){
      if (r==0) k = function(u)    (abs(u)<=1)* 3*(-5*u^2+3)/8
      if (r==2) k = function(u)    (abs(u)<=1)*15*(3*u^2-1)/4
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)   (abs(u)<=1)*(15/32)*(7*u^4-10*u^2+3)
      if (r==2) k = function(u)   (abs(u)<=1)*(105/16)*(6*u^2-5*u^4-1)
    }
  }
  if (v==6) {
    if (kernel=="uni"){
      if (r==0) k = function(u)    (abs(u)<=1)*15*(63*u^4-70*u^2+15)/128
      if (r==2) k = function(u)    (abs(u)<=1)*105*(-45*u^4+42*u^2-5)/32
    }
    if (kernel=="epa"){
      if (r==0) k = function(u)   (abs(u)<=1)*(35/256)*(-99*u^6+189*x^4-105*x^2+15)
      if (r==2) k = function(u)   (abs(u)<=1)*(315/64)*(77*x^6-135*x^4+63*x^2-5)
    }
  }
  Kx = k(x)
  k.v = integrate(function(u) (-1)^v*u^(v)*k(u)/factorial(v), -Inf,Inf)$value
  R.v = integrate(function(u) (k(u))^2,   -Inf,Inf)$value
  out=list(Kx=Kx, k.v=k.v, R.v=R.v)
  return(out)
}

kd.cer.fun = function(x,x0,h,b,v,kernel) {
  n = length(x)
  rho = h/b
  range = max(x)-min(x)

  q.rot = sd(x)*n^(-1/(1+2*v+2*(v+2)))
  K.q = kd.K.fun((x-x0)/q.rot, v=2, r=v+2,   kernel="gau")
  f.r.2 = mean(K.q$Kx)/q.rot^(1+2*v)
  #f.q.rot = kd.K.fun(eval[i], v=2, r=4, kernel="gau")$Kx
  #f.r.2 =  mean((u^4 -6*u^2 + 3)*dnorm(u))/h^(v+3)

  #u = (x-x0)/h
  v.K = kd.K.fun(1, v=v, r=0, kernel=kernel)$k.v
  M.fun = function(u) kd.K.fun(u, v=v,   r=0, kernel=kernel)$Kx - rho^(1+v)*kd.K.fun(rho*u, v=v+2, r=v, kernel=kernel)$Kx*v.K
  K.fun = function(u) kd.K.fun(u, v=v,   r=0, kernel=kernel)$Kx
  L.fun = function(u) kd.K.fun(u, v=v+2, r=v, kernel=kernel)$Kx
  v.fun = function(p,kern) integrate(function(u)   (kern(u))^p,-Inf,Inf)$value
  #m.fun = function(m,kern) ((-1)^m/factorial(m))*(integrate(function(u) u^(m)*kern(u),-Inf,Inf)$value)
  m.fun = function(m,kern) integrate(function(u) u^(m)*kern(u),-Inf,Inf)$value
  v.M.2 = v.fun(2,M.fun)
  v.M.3 = v.fun(3,M.fun)
  v.M.4 = v.fun(4,M.fun)
  m.M.4 = m.fun(4,M.fun)
  m.K.4 = m.fun(4,K.fun)
  m.K.2 = m.fun(2,K.fun)
  m.L.2 = m.fun(2,L.fun)
  
  z = qnorm(0.975)
  q1 = v.M.4*(z^2-3)/6 - v.M.3^2*(z^4-4*z^2+15)/9
  q2 = f.r.2^2 * (m.K.4-rho^(-2)*m.K.2*m.L.2/12)^2 * v.M.2
  q3 = f.r.2   * (m.K.4-rho^(-2)*m.K.2*m.L.2/12)   * v.M.3*(2*z^2)/3

  h <- optimize(function(H) {(H^(-1)*q1 - H^(1+2*(v+2))*q2 + H^(v+2)*q3)^2} , interval=c(.Machine$double.eps, range))$minimum*n^(-1/(v+3))
  return(h)
}


qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x)))
}


lp.bw.fun <- function(V, Bsq, p, v, N, kernel) {
  m1 = function(i,j,k)   integrate(function(x) x^i*x^j*k(x),0,Inf)$value
  m2 = function(i,j,k)   integrate(function(x) x^i*x^j*(k(x))^2,0,Inf)$value
  GAMMA = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m1(i,j,k); return(out)}
  NU    = function(p,k) {out=matrix(NA,p+1,1); for (i in 0:p) out[i+1,1]=m1(i,p+1,k); return(out)}
  PSI   = function(p,k) {out=matrix(NA,p+1,p+1); for (i in 0:p) for (j in 0:p) out[i+1,j+1]=m2(i,j,k); return(out)}
  B.lp  = function(p,k) {out=solve(GAMMA(p,k))%*%NU(p,k); out[1]}
  
  C1.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K))
    C1 = (S.inv%*%NU(p0,K))[v+1]
    return(C1)
  }
  
  C2.fun = function(p0,v,K) {
    S.inv = solve(GAMMA(p0,K))
    C2 = (S.inv%*%PSI(p0,K)%*%S.inv)[v+1,v+1]
    return(C2)
  }
  
  k.fun = function(u){
    if (kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
    if (kernel=="uni") w =          0.5*(abs(u)<=1)
    if (kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
    if (kernel=="gau") w =   dnorm(u)
    return(w)
  }
  
  C1.h <- C1.fun(p0=p, v=v, K=k.fun)
  C2.h <- C2.fun(p0=p, v=v, K=k.fun)
  bw   <- ( ((2*v+1)*C2.h*V) / (2*(p+1-v)*C1.h^2*Bsq*N) )^(1/(2*p+3))
  out=list(bw=bw, C1=C1.h, C2=C2.h)
  return(out)
}

lprobust.res = function(X, y, m, hii, vce, matches, dups, dupsid, d) {
  n = length(y)
  res = matrix(NA,n,1)
  if (vce=="nn") {
    for (pos in 1:n) {
      rpos = dups[pos] - dupsid[pos]
      lpos = dupsid[pos] - 1
      while (lpos+rpos < min(c(matches,n-1))) {
        if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
        else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
        else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
        else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
        else {
          rpos = rpos + dups[pos+rpos+1]
          lpos = lpos + dups[pos-lpos-1]
        }
      }
      ind.J = (pos-lpos):min(c(n,(pos+rpos)))
      y.J   = sum(y[ind.J])-y[pos]
      Ji = length(ind.J)-1
      res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] - y.J/Ji)
    }
  }
  else {
    if (vce=="hc0") w = 1
    else if (vce=="hc1") w = sqrt(n/(n-d))
    else if (vce=="hc2") w = sqrt(1/(1-hii))
    else                 w =      1/(1-hii)
    res[,1] = w*(y-m[,1])
  }
  return(res)
}


lprobust.vce = function(RX, res, C) {
  n = length(C)
  k = ncol(RX)
  M = matrix(0,k,k)
  if (is.null(C)) {
    w = 1
    M = crossprod(c(res)*RX)
  }
  else {	
    clusters = unique(C)
    g        = length(clusters)
    w        = ((n-1)/(n-k))*(g/(g-1))
    for (i in 1:g) {
        ind = C==clusters[i]
        Xi  = RX[ind,]
        ri  = res[ind,]
        M   = M + crossprod(t(crossprod(Xi,ri)),t(crossprod(Xi,ri)))
    }
  }
  return(M)
}

lprobust.bw = function(Y, X, cluster, c, o, nu, o.B, h.V, h.B1, h.B2, scale, vce, nnmatch, kernel, dups, dupsid) {
  
  ### Variance
  dC = 0
  eC = NULL
  w = W.fun((X-c)/h.V, kernel)/h.V
  ind.V = w> 0; eY = Y[ind.V];eX = X[ind.V];eW = w[ind.V]
  n.V = sum(ind.V)
  R.V = matrix(NA,n.V,o+1)
  for (j in 1:(o+1)) R.V[,j] = (eX-c)^(j-1)
  invG.V = qrXXinv(R.V*sqrt(eW))
  beta.V = invG.V%*%crossprod(R.V*eW,eY)
  dups.V = dupsid.V = predicts.V = 0

  if (!is.null(cluster)) {
    dC = 1
    eC =  cluster[ind.V] 
  }
  
  if (vce=="nn") {
    dups.V   = dups[ind.V]
    dupsid.V = dupsid[ind.V]
  }

  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts.V=R.V%*%beta.V
    if (vce=="hc2" | vce=="hc3") {
      hii=matrix(NA,n.V,1)
      for (i in 1:n.V) {
        hii[i] = R.V[i,]%*%invG.V%*%(R.V*eW)[i,]
      }
    }
  }
  res.V = lprobust.res(eX, eY, predicts.V, hii, vce, nnmatch, dups.V, dupsid.V, o+1)
  V.V   = (invG.V%*%lprobust.vce(R.V*eW, res.V, eC)%*%invG.V)[nu+1,nu+1]
  
  ### BIAS
  #Hp = diag(c(1,poly(h.V,degree=o,raw=TRUE)))
  Hp = 0
  for (j in 1:(o+1)) Hp[j] = h.V^((j-1))
  v1 = crossprod(R.V*eW,((eX-c)/h.V)^(o+1))
  v2 = crossprod(R.V*eW,((eX-c)/h.V)^(o+2))
  BConst1 = (Hp*(invG.V%*%v1))[nu+1]
  BConst2 = (Hp*(invG.V%*%v2))[nu+1]

  w = W.fun((X-c)/h.B1, kernel)
  ind = w> 0
  n.B = sum(ind)
  eY = Y[ind];eX = X[ind];eW = w[ind]
  
  if (!is.null(cluster)) eC =  cluster[ind] 
  
  R.B1 = matrix(NA,n.B,o.B+1)
  for (j in 1:(o.B+1)) R.B1[,j] = (eX-c)^(j-1)
  invG.B1 = qrXXinv(R.B1*sqrt(eW))
  beta.B1 = invG.B1%*%crossprod(R.B1*eW,eY)

  BWreg=0
  if (scale>0) {
    #e.B = matrix(0,(o.B+1),1); e.B[o+2]=1
    dups.B = dupsid.B = hii = predicts.B = 0
    if (vce=="nn") {
      dups.B   = dups[ind]
      dupsid.B = dupsid[ind]
    }
    if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
      predicts.B = R.B1%*%beta.B1
      if (vce=="hc2" | vce=="hc3") {
        hii=matrix(NA,n.B,1)
        for (i in 1:n.B) {
          hii[i] = R.B1[i,]%*%invG.B1%*%(R.B1*eW)[i,]
        }
      }
    }
    res.B = lprobust.res(eX, eY, predicts.B, hii, vce, nnmatch, dups.B, dupsid.B,o.B+1)
    V.B = (invG.B1%*%lprobust.vce(R.B1*eW, res.B, eC)%*%invG.B1)[o+2,o+2]
    BWreg = 3*BConst1^2*V.B
  }
  
  w = W.fun((X-c)/h.B2, kernel)
  ind = w> 0
  n.B = sum(ind)
  eY = Y[ind];eX = X[ind];eW = w[ind]
  
  R.B2 = matrix(NA,n.B,o.B+2)
  for (j in 1:(o.B+2)) R.B2[,j] = (eX-c)^(j-1)
  invG.B2 = qrXXinv(R.B2*sqrt(eW))
  beta.B2 = invG.B2%*%crossprod(R.B2*eW,eY)
  
  N  <- length(X)
  B1 <- BConst1%*%(beta.B1[o+2,])
  B2 <- BConst2%*%(beta.B2[o+3,])
  V  <- N*h.V^(2*nu+1)*V.V
  R  <- BWreg
  r  <- 1/(2*o+3)
  rB <- 2*(o+1-nu)
  rV <- 2*nu+1
  bw <- ( (rV*V) / (N*rB*(B1^2 + scale*R)) )^r
  
  output = list(V=V, B1=B1, B2=B2, R=R, r=r, rB=rB, rV=rV, bw=bw)
  return(output)
}


lpbwselect.ce.dpi = function(y, x, h, b, eval, p, q, deriv, rho, kernel, vce, nnmatch, interior, bwregul){
  
  rho <- 1
  bwregul <- 0
  
  N = length(x)
  range = max(x)-min(x)
  
  dups <- dupsid <- hii <- predicts <- NULL


  if (vce=="nn") {
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
  }
  
  if (!is.null(rho)){
    b <- h/rho
  } else {
    rho <- h/b
  }

  X.h <- (x-eval)/h;  X.b <- (x-eval)/b
  K.h <- W.fun(X.h,kernel);  L.b <- W.fun(X.b,kernel)
  ind.h <- K.h>0;    ind.b <- L.b>0
  ind <- ind.h
  if (h>b) ind=ind.h
  eN   <- sum(ind)
  eY   <-   y[ind];  eX   <-   x[ind]
  eX.h <- X.h[ind];  eX.b <- X.b[ind]
  eK.h <- K.h[ind];  eL.b <- L.b[ind]

  W.p <- eK.h/h;  W.q <- eL.b/b
  R.p.2 <- matrix(NA,eN,(p+3))
  for (j in 1:(p+3))  R.p.2[,j] <- eX.h^(j-1)
  R.p.1 <- R.p.2[,1:(p+2)]
  R.p   <- R.p.2[,1:(p+1)]
  R.q   <- matrix(NA,eN,(q+1))
  for (j in 1:(q+1))  R.q[,j] <- eX.b^(j-1)

  L.p.1 <- crossprod(R.p*W.p, eX.h^(p+1))/eN
  L.p.2 <- crossprod(R.p*W.p, eX.h^(p+2))/eN
  L.p.3 <- crossprod(R.p*W.p, eX.h^(p+3))/eN
  L.q.1 <- crossprod(R.q*W.q, eX.b^(q+1))/eN
  L.q.2 <- crossprod(R.q*W.q, eX.b^(q+2))/eN
  L.q.3 <- crossprod(R.q*W.q, eX.b^(q+3))/eN
  invG.p <- eN*qrXXinv(R.p*sqrt(W.p))
  invG.q <- eN*qrXXinv(R.q*sqrt(W.q))

  edups <- edupsid <- 0	
  if (vce=="nn") {
    for (j in 1:eN) {
      edups[j]=sum(eX==eX[j])
    }
    j=1
    while (j<=eN) {
      edupsid[j:(j+edups[j]-1)] <- 1:edups[j]
      j <- j+edups[j]
    }
  }
  
  hii=predicts.p=predicts.q=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    H.q <- 0
    for (j in 1:(q+1)) H.q[j] <- b^(-(j-1))
    beta.q <- H.q*invG.q%*%crossprod(R.q*W.q,eY)/eN
    r.q <- matrix(NA,eN,q+1)
    predicts <- 0
    for (j in 1:(q+1))  r.q[,j] <- (eX-eval)^(j-1)
    for (j in 1:eN) predicts[j] <- r.q[j,]%*%beta.q
    if (vce=="hc2" | vce=="hc3") {
      hii <- matrix(NA,eN,1)
      for (j in 1:eN) hii[j] <- (R.p[j,]%*%invG.p%*%(R.p*W.p)[j,])/eN
    }
  }

  res.q <- lprobust.res(eX, eY, as.matrix(predicts), hii, vce, nnmatch, edups, edupsid, q+1)
 
  ### Bias
  k <- p+3
  r.k <- matrix(NA,N,k+3)
  for (j in 1:(k+3))  r.k[,j] <- x^(j-1)
  gamma <- lm(y~r.k-1)
  m.p.3 <- gamma$coeff[p+4]*factorial(p+3) + gamma$coeff[p+5]*factorial(p+4)*eval + gamma$coeff[p+6]*factorial(p+5)*eval^2/2
  #r.k <- r.k[,1:(k+2)]
  gamma <- lm(y~r.k[,1:(k+2)]-1)
  m.p.2 <- gamma$coeff[p+3]*factorial(p+2) + gamma$coeff[p+4]*factorial(p+3)*eval + gamma$coeff[p+5]*factorial(p+4)*eval^2/2

  if (is.na(m.p.3)) m.p.3 <- lprobust(y=y, x=x, h=range, eval=eval, p=p+4, deriv=(p+3), kernel=kernel, vce=vce)$Estimate[5]
  if (is.na(m.p.2)) m.p.2 <- lprobust(y=y, x=x, h=range, eval=eval, p=p+4, deriv=(p+2), kernel=kernel, vce=vce)$Estimate[5]
  
  e.p.1 <- matrix(0,(q+1),1); e.p.1[p+2]=1
  e.0   <- matrix(0,(p+1),1); e.0[1]=1

  q.terms <- lpbwce(y = eY, x = eX,  K=eK.h, L=eL.b, res=res.q, c = eval, p=p, q=q, h=h, b=b, deriv=deriv, fact=factorial(deriv))
  q1.rbc  <- q.terms$q1rbc
  q2.rbc  <- q.terms$q2rbc
  q3.rbc  <- q.terms$q3rbc

  V.reg=0
  if (bwregul>0) {
    o <- p+2
    R.reg <-  matrix(NA,eN,o+1)
    for (j in 1:(o+1)) R.reg[,j] = (eX-eval)^(j-1)
    invG.reg = qrXXinv(R.reg*sqrt(eK.h))
    V.reg = (invG.reg%*%lprobust.vce(R.reg*eK.h, res.q)%*%invG.reg)[o+1,o+1]
    #beta.reg = invG.reg%*%crossprod(R.reg*eK.h,eY)
    #Hp = 0
    #for (j in 1:o) Hp[j] = h^((j-1))
    #v = crossprod(R.reg*eK.h,((eX-eval)/h)^o)
    #BConst = (Hp*(invG.reg%*%v))[1]
  }  
  
  if (interior==TRUE) {
    eta.bc1  <- (t(e.0)%*%invG.p)%*%(   (m.p.2/factorial(p+2))*L.p.2/h     + (m.p.3/factorial(p+3))*L.p.3 )
    eta.bc2  <- rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p)%*%L.p.1%*%t(e.p.1)%*%invG.q%*%( (m.p.2/factorial(p+2))*L.q.1/b     + (m.p.3/factorial(p+3))*L.q.2 )
    eta.bc   <- (eta.bc1 - eta.bc2)
    Reg      <- 3*(t(e.0)%*%invG.p)%*%(L.p.2/h)^2*V.reg
    
    H.bc     <- function(H) {abs(H^(-1)*q1.rbc + H^(1+2*(p+3))*(eta.bc^2 + bwregul*Reg)*q2.rbc + H^(p+3)*(eta.bc + bwregul*Reg)*q3.rbc)}
    h.bc     <- optimize(H.bc , interval=c(.Machine$double.eps, range))
    h.ce.dpi <- h.bc$minimum*N^(-1/(p+4))
  }

  if (interior==FALSE) {
    eta.bc.1 <- (t(e.0)%*%invG.p)%*%( L.p.2 - rho^(-1)*L.p.1%*%t(e.p.1)%*%invG.q%*%L.q.1)*(m.p.2/factorial(p+2))
    eta.bc.2 <- (t(e.0)%*%invG.p)%*%( L.p.3 - rho^(-2)*L.p.1%*%t(e.p.1)%*%invG.q%*%L.q.2)*(m.p.3/factorial(p+3))
    Reg      <- 3*((t(e.0)%*%invG.p)%*%( L.p.2 - rho^(-1)*L.p.1%*%t(e.p.1)%*%invG.q%*%L.q.1))^2*V.reg
    
    phi.z    <- dnorm(1.96)
    E1       <- q1.rbc*phi.z
    E2       <- eta.bc.1
    E3       <- eta.bc.2
    E4       <- q3.rbc*phi.z*eta.bc.1
    E5       <- q3.rbc*phi.z*eta.bc.2
    H.bc <-function(H) {abs(E1/(N*H) + N*H^(2*p+5)*q2.rbc*phi.z*(E2 + H*E3 + bwregul*Reg)^2 + H^(p+2)*(E4 + H*E5 + bwregul*Reg))}
    h.ce.dpi <- optimize(H.bc , interval=c(.Machine$double.eps, range))$minimum
  }
  out <- list(h=h.ce.dpi)
  return(out)
}

lpbwselect.imse.dpi = function(y, x, cluster, p, q, deriv, kernel, bwcheck, bwregul, imsegrid, vce, nnmatch, interior){

  N <- length(x)
  x.max <- max(x);  x.min <- min(x)
  range <- x.max - x.min
  eval <- seq(x.min, x.max, length.out=imsegrid)
  #eval <- quantile(x, probs = seq(0.05, 0.95, 0.025))
  neval = length(eval)
  even <- (p-deriv)%%2==0
  V.h <- B.h <- V.b <- B.b <- 0

  for (i in 1:neval) {
    est <- lpbwselect.mse.dpi(y=y, x=x, cluster=cluster, eval=eval[i], p=p, q=q, deriv=deriv, kernel=kernel,
                              bwcheck=bwcheck, bwregul=bwregul, vce=vce, nnmatch=nnmatch, interior=interior)
    
    V.h[i] <- est$V.h;    B.h[i] <- est$B.h
    V.b[i] <- est$V.b;    B.b[i] <- est$B.b
  }

  if (even==FALSE | interior==TRUE) {
    b.imse.dpi <- (mean(V.b)/(N*mean(B.b)))^(1/(2*q+3))
    h.imse.dpi <- (mean(V.h)/(N*mean(B.h)))^(1/(2*p+3))
  }

  if (even==TRUE & interior==FALSE) {
    b.bw.fun   <- function(H) {abs(H^(2*q+2-2*(p+1))*mean(B.b) + mean(V.b)/(N*H^(1+2*(p+1))))}
    b.imse.dpi <- optimize(b.bw.fun, interval=c(.Machine$double.eps, range))$minimum

    h.bw.fun   <- function(H) {abs(H^(2*p+2-2*deriv)*mean(B.h) + mean(V.h)/(N*H^(1+2*deriv)))}
    h.imse.dpi <- optimize(h.bw.fun, interval=c(.Machine$double.eps, range))$minimum
  }

  out <- list(h=h.imse.dpi, b=b.imse.dpi)
  return(out)
}


lpbwselect.imse.rot = function(y, x, p, deriv, kernel, imsegrid){
  
  even <- (p-deriv)%%2==0
  eval <- quantile(x, probs = seq(0.05, 0.95, 0.025))
  #eval  <- seq(min(x), max(x), length.out=imsegrid)
  neval <- length(eval)
  x.max <- max(x);  x.min <- min(x)
  range <- x.max - x.min
  N     <- length(x)
  V     <- B <- 0

  for (i in 1:neval) {
    est  <- lpbwselect.mse.rot(y=y, x=x, eval=eval[i], p=p, deriv=deriv, kernel=kernel)
    V[i] <- est$V;    B[i] <- est$B
  }
  
  if (even==FALSE) {
    h.imse.rot <- (mean((1+2*deriv)*V)/(N*mean(2*(p+1-deriv)*B^2)))^(1/(2*p+3))
  } else {
    h.bw.fun   <- function(H) {abs(H^(2*p+2-2*deriv)*mean(B^2) + mean(V)/(N*H^(1+2*deriv)))}
    h.imse.rot <- optimize(h.bw.fun, interval=c(.Machine$double.eps, range))$minimum
  }
  out <- list(h=h.imse.rot)
  return(out)
}


lpbwselect.mse.dpi = function(y, x, cluster, eval, p, q, deriv, kernel, bwcheck, bwregul, vce, nnmatch, interior){

  even <- (p-deriv)%%2==0

                     C.c         <- 2.34
  if (kernel=="uni") C.c         <- 1.843
  if (kernel=="tri") C.c         <- 2.576
  if (kernel=="gau") C.c         <- 1.06

  x.iq  <- quantile(x,.75) - quantile(x,.25)
  x.max <- max(x);  x.min <- min(x)
  range <- x.max - x.min
  N     <- length(x)

  bw.max = max(abs(eval-x.min),abs(eval-x.max))
  
  c.bw  <- c(C.c*min(c(sd(x),x.iq/1.349))*N^(-1/5))
  c.bw <- min(c.bw, bw.max)
  
  dups <- dupsid <- hii <- predicts <- NULL
  if (vce=="nn") {
    order.x = order(x)
    x = x[order.x];   y = y[order.x]
    for (j in 1:N) {
      dups[j]=sum(x==x[j])
    }
    j=1
    while (j<=N) {
      dupsid[j:(j+dups[j]-1)] = 1:dups[j]
      j = j+dups[j]
    }
  }

  bw.adj <- 0
  if (!is.null(bwcheck)) {
    bw.min    <- sort(abs(x-eval))[bwcheck]
    #nh        <- sum(abs(x-eval) <= c.bw)
    c.bw      <- max(c.bw, bw.min)
    bw.adj <- 1
  }

   C.d1 <- lprobust.bw(y, x, cluster, c=eval, o=q+1, nu=q+1, o.B=q+2, h.V=c.bw, h.B1=range, h.B2=range, 0, vce, nnmatch, kernel, dups, dupsid)
   if (even==FALSE | interior==TRUE) {
       bw.mp2 <- c(C.d1$bw)
    }
    else{
      bw.fun  <-function(H) {abs(H^(2*(q+1)+2-2*(q+1))*(C.d1$B1 + H*C.d1$B2)^2 + C.d1$V/(N*H^(1+2*(q+1))))}
      bw.mp2 <- optimize(bw.fun, interval=c(.Machine$double.eps, range))$minimum
   }
   

    C.d2 <- lprobust.bw(y, x, cluster, c=eval, o=q+2, nu=q+2, o.B=q+3, h.V=c.bw, h.B1=range, h.B2=range, 0, vce, nnmatch, kernel, dups, dupsid)
    if (even==FALSE | interior==TRUE) {
      bw.mp3 <- c(C.d2$bw)
    }
    else {
      bw.fun  <-function(H) {abs(H^(2*(q+2)+2-2*(q+2))*(C.d2$B1 + H*C.d2$B2)^2 + C.d2$V/(N*H^(1+2*(q+2))))}
      bw.mp3 <- optimize(bw.fun, interval=c(.Machine$double.eps, range))$minimum
    }
    
    # adjust
    bw.mp2 <- min(bw.mp2, bw.max) 
    bw.mp3 <- min(bw.mp3, bw.max) 
    
    if (!is.null(bwcheck)) {
      bw.mp2 <- max(bw.mp2, bw.min)
      bw.mp3 <- max(bw.mp3, bw.min)
    }
    
    
    ### Select preliminar bw b
    C.b <- lprobust.bw(y, x, cluster, c=eval, o=q, nu=p+1, o.B=q+1, h.V=c.bw, h.B1=bw.mp2, h.B2=bw.mp3, bwregul, vce, nnmatch, kernel, dups, dupsid)
    
    if (even==FALSE | interior==TRUE) {
      b.mse.dpi <- c(C.b$bw)
    }
    else {
      b.bw.fun = function(H) {abs(H^(2*q+2-2*(p+1))*(C.b$B1 + H*C.b$B2 + bwregul*C.b$R)^2 + C.b$V/(N*H^(1+2*(p+1))))}
      b.mse.dpi <- optimize(b.bw.fun, interval=c(.Machine$double.eps, range))$minimum
    }
    
    b.mse.dpi <- min(b.mse.dpi, bw.max) 
    if (!is.null(bwcheck)) b.mse.dpi <- max(b.mse.dpi, bw.min)
  
    bw.mp1 = b.mse.dpi
    
    ### Select final bw h
    C.h <- lprobust.bw(y, x, cluster, c=eval, o=p, nu=deriv, o.B=q, h.V=c.bw, h.B1=bw.mp1, h.B2=bw.mp2, bwregul, vce, nnmatch, kernel, dups, dupsid)
    
    if (even==FALSE | interior==TRUE) {
      h.mse.dpi <- c(C.h$bw)
    } else {
      h.bw.fun = function(H) {abs(H^(2*p+2-2*deriv)*(C.h$B1 + H*C.h$B2 + bwregul*C.h$R)^2 + C.h$V/(N*H^(1+2*deriv)))}
      h.mse.dpi <- optimize(h.bw.fun, interval=c(.Machine$double.eps, range))$minimum
    }
    
    h.mse.dpi <- min(h.mse.dpi, bw.max) 
    if (!is.null(bwcheck)) h.mse.dpi <- max(h.mse.dpi, bw.min)
  
    if (even==FALSE | interior==TRUE) {
      V.h <- C.h$rV*C.h$V;    B.h <- C.h$rB*C.h$B1^2
      V.b <- C.b$rV*C.b$V;    B.b <- C.b$rB*C.b$B1^2
    } else {
      V.h <- C.h$V;    B.h <- (C.h$B1 + h.mse.dpi*C.h$B2)^2
      V.b <- C.b$V;    B.b <- (C.b$B1 + b.mse.dpi*C.b$B2)^2
    }

  out <- list(h=h.mse.dpi, b=b.mse.dpi, V.h=V.h, B.h=B.h, V.b=V.b, B.b=B.b)
  return(out)
}

lpbwselect.mse.rot = function(y, x, eval, p, deriv, kernel){

  even <- (p-deriv)%%2==0
  
  C.c         <- 2.34
  if (kernel=="uni") C.c         <- 1.843
  if (kernel=="tri") C.c         <- 2.576

  x.iq  <- quantile(x,.75) - quantile(x,.25)
  x.max <- max(x);  x.min <- min(x)
  range <- x.max - x.min
  N     <- length(x)

  c.bw  <- C.c*min(c(sd(x),x.iq/1.349))*N^(-1/5)
  n.h1  <- sum(abs(x-eval) <= c.bw)
  f.hat <- n.h1/(2*N*c.bw)
  
  k   <- p+3
  r.k <- matrix(NA,N,k+1)
  for (j in 1:(k+1))  r.k[,j] <- x^(j-1)
  gamma.p <- lm(y~r.k-1)
  s2.hat  <- sigma(gamma.p)^2

  #k   <- q+3
  #r.k <- matrix(NA,N,k+1)
  #for (j in 1:(k+1))  r.k[,j] <- x^(j-1)
  #gamma.q <-lm(y~r.k-1)
  #s2.q    <- sigma(gamma.q)^2

  m.p.1 <- lprobust(y=y, x=x, h=range, eval=eval, p=p+3, deriv=(p+1), kernel=kernel, vce="nn")$Estimate[5]
  m.p.2 <- lprobust(y=y, x=x, h=range, eval=eval, p=p+3, deriv=(p+2), kernel=kernel, vce="nn")$Estimate[5]
  
  #m.q.1 <- lprobust(y=y, x=x, h=range, eval=eval, p=q+3, deriv=(q+1), kernel=kernel, vce=vce)$Estimate[5]
  #m.q.2 <- lprobust(y=y, x=x, h=range, eval=eval, p=q+3, deriv=(q+1), kernel=kernel, vce=vce)$Estimate[5]
  
  #2*gamma.p$coeff[p+2] + 3*2*gamma.p$coeff[p+3]*eval + 4*3*gamma.p$coeff[p+4]*eval^2
  #m.p.1     <- gamma.p$coeff[p+1]*factorial(p+1) + gamma.p$coeff[p+2]*factorial(p+2)*eval + gamma.p$coeff[p+3]*factorial(p+3)*eval^2/2
  #m.q.1     <- gamma.q$coeff[q+1]*factorial(q+1) + gamma.q$coeff[q+2]*factorial(q+2)*eval + gamma.q$coeff[q+3]*factorial(q+3)*eval^2/2


  #h.mse.rot <- lp.bw.fun(s2.p/f0.pilot, (m.p.1/factorial(p+1))^2, p, deriv, N, kernel)
  #b.mse.rot <- lp.bw.fun(s2.q/f0.pilot, (m.q.1/factorial(q+1))^2, q, p+1  , N, kernel)

  bw <- lp.bw.fun(s2.hat/f.hat, (m.p.1/factorial(p+1))^2, p, deriv, N, kernel)
  
  B1 = bw$C1*m.p.1/factorial(p+1)
  B2 = bw$C1*m.p.2/factorial(p+2)
  V  = bw$C2*s2.hat/f.hat

  if (even==FALSE) {
    h.mse.rot <- bw$bw
    B <- B1
  } else {
    h.bw.fun  <- function(H) {abs(H^(2*p+2-2*deriv)*(B1 + H*B2)^2 + V/(N*H^(1+2*deriv)))}
    h.mse.rot <- optimize(h.bw.fun, interval=c(.Machine$double.eps, range))$minimum
    B <- B1 + h.mse.rot*B2
  }
  
  out <- list(h=h.mse.rot, V=V, B=B)
  return(out)
}




