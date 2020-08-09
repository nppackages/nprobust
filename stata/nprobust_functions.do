*!version 0.3.1  06-20-2020
 

 
capture mata mata drop nprobust_w_fun()
mata
real matrix nprobust_w_fun(real matrix u, string kernel)
{
  if (kernel=="epa") w = 0.75:*(1:-u:^2):*(abs(u):<=1)
  if (kernel=="uni") w =             0.5:*(abs(u):<=1)
  if (kernel=="tri") w =     (1:-abs(u)):*(abs(u):<=1)
  if (kernel=="gau") w = normalden(u)
  return(w)  
}
mata mosave nprobust_w_fun(), replace
end

capture mata mata drop nprobust_K_fun()
mata
real matrix nprobust_K_fun(real matrix u, string kernel)
{
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = 0.75:*(1:-u:^2):*(abs(u):<=1)
	}
	else {
		kx = 0.5:*(abs(u):<=1)
	}
	return(kx)  
}
mata mosave nprobust_K_fun(), replace
end		 
		
		
capture mata mata drop nprobust_L_fun()
mata
real matrix nprobust_L_fun(real matrix u, string kernel)
{
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = (abs(u):<=1):*(105/16):*(6:*u:^2-5:*u:^4:-1)
	}
	else {
		kx = (abs(u):<=1):*15:*(3:*u:^2:-1)/4
	}
	return(kx)  
}
mata mosave nprobust_L_fun(), replace
end		
		 
mata
void nprobust_K(u, kernel, kx, vK, mK)
         {
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = 0.75:*(1:-u:^2):*(abs(u):<=1)
		mK = 0.5
		vK = 0.5
	}
	else {
		kx = 0.5:*(abs(u):<=1)
		mK = 0.5
		vK = 0.5
	}
         }
		 
void nprobust_L(u, kernel, kx, vK, mK)
         {
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = (abs(u):<=1):*(105/16):*(6:*u:^2-5:*u:^4:-1)
		mK = 0.5
		vK = 0.5
	}
	else {
		kx = (abs(u):<=1):*15:*(3:*u:^2:-1)/4
		mK = 0.5
		vK = 0.5
	}
         }
end
		 





/*

capture mata mata drop nprobust_kd_bw_cer()
mata
real matrix nprobust_kd_bw_cer(real matrix X, real scalar x0, real scalar h, real scalar b, real scalar v, string kernel)
{
h.cer.den = function(x,x0,h,b,v,kernel) {
  n = length(x)
  rho = h/b
  range = max(x)-min(x)
  
  q_rot = sd(x)*n^(-1/(1+2*v+2*(v+2)))
  K_q = kern.fun((x-x0)/q.rot, v=2, r=v+2,   kernel="gau")    
  f_r_2 = mean(K.q$Kx)/q.rot^(1+2*v)
  
  v_K = kern.fun(1, v=v, r=0, kernel=kernel)$k.v
  M_fun = function(u) kern.fun(u, v=v,r=0, kernel=kernel)$Kx - rho^(1+v)*kern.fun(rho*u, v=v+2, r=v, kernel=kernel)$Kx*v.K
  v_fun = function(p,kern) integrate(function(u)   (kern(u))^p,-Inf,Inf)$value
  m_fun = function(m,kern) ((-1)^m/factorial(m))*(integrate(function(u) u^(m)*kern(u),-Inf,Inf)$value)
  v_M_2 = v.fun(2,M.fun)
  v_M_3 = v.fun(3,M.fun)
  v_M_4 = v.fun(4,M.fun)
  m_M_4 = m.fun(4,M.fun)
  
  z = qnorm(0.975)
  q1 = v_M_4*(z^2-3)/6 - v.M.3^2*(z^4-4*z^2+15)/9
  q2 = f_r_2^2 * m_M_4^2 * v_M_2
  q3 = f_r_2   * m_M_4   * v_M_3*(2*z^2)/3
  
  h = optimize(function(H) {(H^(-1)*q1 - H^(1+2*(v+2))*q2 + H^(v+2)*q3)^2} , interval=c(.Machine$double.eps, range))$minimum*n^(-1/(v+3))
  return(h)
}
}
mata mosave nprobust_kd_bw_cer(), replace
end


*/


*********************************
  
capture mata mata drop nprobust_lp_res()
mata
real matrix nprobust_lp_res(real matrix X, real matrix y, real matrix m, real matrix hii, string vce, real scalar matches, dups, dupsid, real scalar d)
{
n = length(y)
res = J(n,1,.)		
if (vce=="nn") {
	for (pos=1; pos<=n; pos++) {
		rpos = dups[pos] - dupsid[pos]
		lpos = dupsid[pos] - 1
		while (lpos+rpos < min((matches,n-1))) {
			if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
			else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
			else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
			else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
			else {
				rpos = rpos + dups[pos+rpos+1]
				lpos = lpos + dups[pos-lpos-1]
			}
		}
		ind_J = (pos-lpos)::(pos+rpos)
		y_J   = sum(y[ind_J])-y[pos]
		Ji = length(ind_J)-1
		res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] :- y_J/Ji)
	}		
}
else if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
 	 if (vce=="hc0") w = 1
	 if (vce=="hc1") w = sqrt(n/(n-d))
	 if (vce=="hc2") w = sqrt(1:/(1:-hii))
	 if (vce=="hc3") w =      1:/(1:-hii)
	 res[,1] = w:*(y-m[,1])	 
}
return(res)
}
mata mosave nprobust_lp_res(), replace
end

*************************************************************************************************************************************************************
capture mata mata drop nprobust_lp_kweight()
mata
real matrix nprobust_lp_kweight(real matrix X, real scalar c, real scalar h, string kernel)
{
u = (X:-c)/h
	if (kernel=="epanechnikov" | kernel=="epa") {
		w = (0.75:*(1:-u:^2):*(abs(u):<=1))/h
	}
	else if (kernel=="uniform" | kernel=="uni") {
		w = (0.5:*(abs(u):<=1))/h
	}
	else if (kernel=="gaussian" | kernel=="gau") {
		w = normalden(u)/h
	}
	else {
		w = ((1:-abs(u)):*(abs(u):<=1))/h
	}	
return(w)	
}
mata mosave nprobust_lp_kweight(), replace
end

*************************************************************************************************************************************************************
capture mata mata drop nprobust_lp_bw()
mata
real matrix nprobust_lp_bw(real matrix Y, real matrix X, real matrix C, real scalar c, real scalar o, real scalar nu, real scalar o_B, real scalar h_V, real scalar h_B1, real scalar h_B2, real scalar scale, string vce, real scalar nnmatch, string kernel, dups, dupsid)
{
	dC = eC = indC = 0
	w = nprobust_lp_kweight(X, c, h_V, kernel)
	ind_V = selectindex(w:> 0); eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
	n_V = length(ind_V)
	D_V = eY
	R_V = J(n_V,o+1,.)
	for (j=1; j<=(o+1); j++) R_V[.,j] = (eX:-c):^(j-1)
	invG_V = cholinv(quadcross(R_V,eW,R_V))
	e_v = J((o+1),1,0); e_v[nu+1]=1
	
	if (rows(C)>1) {
		dC = 1
		eC =  C[ind_V] 
		indC = order(eC,1) 
	}
	
	beta_V = invG_V*quadcross(R_V:*eW,D_V)	
	dups_V=dupsid_V=predicts_V=0
	if (vce=="nn") {
		dups_V   = dups[ind_V]
		dupsid_V = dupsid[ind_V]
	}
	if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
		predicts_V=R_V*beta_V
		if (vce=="hc2" | vce=="hc3") {
			hii=J(n_V,1,.)	
				for (i=1; i<=n_V; i++) {
					hii[i] = R_V[i,]*invG_V*(R_V:*eW)[i,]'
				}
		}
	}	
	res_V = nprobust_lp_res(eX, eY, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1)
	V_V = (invG_V*nprobust_lp_vce(R_V:*eW, res_V, eC, indC)*invG_V)[nu+1,nu+1]
	
	Hp = J(o+1, 1, 1)
	for (j=1; j<=(o+1); j++) Hp[j] = h_V^((j-1))
	
	v1 = quadcross(R_V:*eW,((eX:-c):/h_V):^(o+1))
	v2 = quadcross(R_V:*eW,((eX:-c):/h_V):^(o+2))
	BConst1 = (diag(Hp)*(invG_V*v1))[nu+1]
	BConst2 = (diag(Hp)*(invG_V*v2))[nu+1]
		
	w = nprobust_lp_kweight(X, c, h_B1, kernel)
	ind = selectindex(w:> 0) 
	n_B = length(ind)
	eY = Y[ind];eX = X[ind];eW = w[ind]

	if (dC==1) {
		eC   = C[ind]
		indC = order(eC,1) 
	}	
	
	R_B1 = J(n_B,o_B+1,.)
	for (j=1; j<=(o_B+1); j++) R_B1[.,j] = (eX:-c):^(j-1)
	invG_B1 = cholinv(quadcross(R_B1,eW,R_B1))
	beta_B1 = invG_B1*quadcross(R_B1:*eW,eY)	
	
  
  	BWreg=0
	if (scale>0) {
		e_B = J((o_B+1),1,0); e_B[o+2]=1
		dups_B=dupsid_B=hii=predicts_B=0
		if (vce=="nn") {
			dups_B   = dups[ind]
			dupsid_B = dupsid[ind]
		}
		if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
			predicts_B=R_B1*beta_B1
			if (vce=="hc2" | vce=="hc3") {
				hii=J(n_B,1,.)	
					for (i=1; i<=n_B; i++) {
						hii[i] = R_B1[i,]*invG_B1*(R_B1:*eW)[i,]'
				}
			}
		}	
		res_B = nprobust_lp_res(eX, eY, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B, o_B+1)		        
		V_B = (invG_B1*nprobust_lp_vce(R_B1:*eW, res_B, eC, indC)*invG_B1)[o+2,o+2]
		BWreg = 3*BConst1^2*V_B
	}
	
    w = nprobust_lp_kweight(X, c, h_B2, kernel)
	ind = selectindex(w:> 0) 
	n_B = length(ind)
	eY = Y[ind];eX = X[ind];eW = w[ind]
  
	R_B2 = J(n_B,o_B+2,.)
	for (j=1; j<=(o_B+2); j++) R_B2[.,j] = (eX:-c):^(j-1)
	invG_B2 = cholinv(quadcross(R_B2,eW,R_B2))
	beta_B2 = invG_B2*quadcross(R_B2:*eW,eY)	
	
	N  = length(X)
	B1 =  BConst1*(beta_B1[o+2,]')
	B2 =  BConst2*(beta_B2[o+3,]')
	V  = N*h_V^(2*nu+1)*V_V
	R = BWreg
	B_rate = 2*(o+1-nu)
	V_rate = (2*nu+1)
	h_rate = 1/(2*o+3)
	bw = ((V_rate*V)/(N*B_rate*(B1^2+scale*R)))^h_rate
	
	return(B1, B2, V, R, B_rate, V_rate, h_rate, bw)	
}
mata mosave nprobust_lp_bw(), replace
end


****************************************************
capture mata mata drop nprobust_lp_vce()
mata
real matrix nprobust_lp_vce(real matrix RX, real matrix res, real matrix C, real matrix ind)
{	
	k = cols(RX)
	M = J(k,k,0)
	n = length(C)
	if (n==1) {
		w = 1
		SS = res:^2
		M  = quadcross(RX,SS,RX)
	}
	else {	
		C_o   = C[ind]
		RX_o  = RX[ind,]
		res_o = res[ind,]
		info  = panelsetup(C_o,1)
		g     = rows(info)
		w     = ((n-1)/(n-k))*(g/(g-1))
		for (i=1; i<=g; i++) {
				Xi = panelsubmatrix(RX_o,  i, info)
				ri = panelsubmatrix(res_o, i, info)
				M = M + quadcross(quadcross(Xi,ri)',quadcross(Xi,ri)')
		}
	}
	return(w*M)			
}
mata mosave nprobust_lp_vce(), replace
end


	
	
	
******************************************************************************************************
capture mata mata drop nprobust_lp_q()
mata
real matrix nprobust_lp_q(real matrix y, real matrix x, real matrix K, real matrix L, real matrix res, real scalar c, real scalar h, real scalar b, real scalar p, real scalar q, real scalar deriv)
{

  N = length(y)
  rho = h / b
  
  Wp = K:/h
  Wq = L:/b
  Xh = (x:-c):/h
  Xb = (x:-c):/b
  Rq = J(N,(q+1),.)
  Rp = J(N,(p+1),.)
  for (i=1; i<=(q+1); i++) {
	Rq[.,i] = Xb:^(i-1)
  }
	Rp = Rq[,1::(p+1)]	 
  
  Lp1 = quadcross(Rp:*Wp,Xh:^(p+1))/N
  iGq  = N*cholinv(quadcross(Rq,Wq,Rq))
  iGp  = N*cholinv(quadcross(Rp,Wp,Rp))
  
  ep1 = J(q+1,1,0)
  ep1[p+2]=1
  e0 = J(p+1,1,0)
  e0[deriv+1]=factorial(deriv+1)
        
  lus0 = e0'*iGp*(diag(K)*Rp)'
  lbc0 = lus0 - rho^(p+1)*(e0'*iGp)*Lp1*ep1'*iGq*(diag(L)*Rq)'
  vx = res:^2
  
  sums2=0
  for (i=1; i<=N; i++) {
     sums2 = sums2+ (lbc0[i]^2)*vx[i]
  }
  s2=sums2/(N*h)
  
  Krrp  = J(p+1,p+1,0)
  Krxip = J(1,p+1,0)
  Krxp  = J(1,p+1,0)
  Lrrq  = J(q+1,q+1,0)
  
  for (i=1; i<=N; i++) {
    Rpi = Rp[i,.]
    Rqi = Rq[i,.]
    Krrp  = Krrp +  K[i]*Rpi'*Rpi
    Lrrq  = Lrrq +  L[i]*Rqi'*Rqi
    Krxip = Krxip + K[i]*Rpi*Xh[i]^(p+1) 
	for (j=1; j<=N; j++) {
       if (j != i) Krxp = Krxp + K[i]*Rpi*Xh[j]^(p+1) 
    }
  }
  
  EKrrp = Krrp/N
  EKrxp = Krxp/(N*(N-1))
  EKrxip = Krxip/N
  ELrrq = Lrrq/N

  q1=J(1,1,0)
  q2=J(1,1,0)
  q3=J(1,1,0)
  q4=J(1,1,0)
  q5a=J(1,q+1,0)
  q5b=J(q+1,1,0)
  q6=J(1,1,0)
  q7a=J(1,q+1,0)
  q7b=J(q+1,q+1,0)
  q7c=J(q+1,1,0)
  q8=J(1,1,0)
  q9=J(1,1,0)
  q10=J(1,1,0)
  q11=J(1,1,0)
  q12=J(1,1,0)
  q3a=J(1,1,0)  


  for (i=1; i<=N; i++) {
  
    Rpi = Rp[i,.]
    Rqi = Rq[i,.]
    
    q1 = q1 + (lbc0[i]*res[i])^3
    
    lus1 = factorial(deriv+1)*iGp[deriv+1,.] * (EKrrp - K[i]*Rpi'*Rpi)*iGp*K[i]*Rpi'
    T1   = iGp[deriv+1,.] * ((EKrrp - K[i]*Rpi'*Rpi)*iGp*Lp1*ep1')    *(iGq*L[i]*Rqi')
    T2   = iGp[deriv+1,.] * ((K[i]*Rpi*Xh[i]^(p+1) - EKrxip)')*ep1'*(iGq*L[i]*Rqi')
    T3   = iGp[deriv+1,.] * ((Lp1*ep1'*iGq)*(ELrrq - L[i]*Rqi'*Rqi))  *(iGq*L[i]*Rqi')
    lbc1 = lus1 - factorial(deriv+1)*rho^(p+1)*(T1 + T2 + T3)
    
    q2 = q2 + lbc1*lbc0[i]*res[i]^2
    
    q3 = q3 + (lbc0[i]^4)*((res[i]^4)-vx[i]^2)
    
    q4 = q4 + (lbc0[i]^2)*(Rqi*iGq*L[i]*Rqi')*(res[i]^2)

    q5a = q5a + (lbc0[i]^3)*(Rqi*iGq)*(res[i]^2)
    q5b = q5b + L[i]*Rqi'*lbc0[i]*(res[i]^2)
    
    q7a = q7a + lbc0[i]*(res[i]^2)*L[i]*Rqi*iGq
    q7b = q7b + (lbc0[i]^2)*Rqi'*Rqi*iGq
    q7c = q7c + lbc0[i]*(res[i]^2)*L[i]*Rqi'
    
    q8  = q8 + (lbc0[i]^4)*(res[i]^4)
    q9  = q9 + ((lbc0[i]^2)*vx[i]-h*s2)*(lbc0[i]*res[i])^2
    
    q12 = q12 + ((lbc0[i]^2)*vx[i]-h*s2)^2

    q3a = q3a + (lbc0[i]*res[i])^3

    for (j=1; j<=N; j++) {  
      if (j != i) {
        Rpj = Rp[j,.]
        Rqj = Rq[j,.]
        
        lus1 = iGp[1,.] *  (EKrrp - K[j]*Rpj'*Rpj)*iGp*K[i]*Rpi'
        T1   = iGp[1,.] * ((EKrrp - K[j]*Rpj'*Rpj)*iGp*Lp1*ep1')    *(iGq*L[i]*Rqi')
        T2   = iGp[1,.] * ((K[j]*Rpj*Xh[i]^(p+1) - EKrxp)')*ep1' *(iGq*L[i]*Rqi')
        T3   = iGp[1,.] * ((Lp1*ep1'*iGq)*(ELrrq - L[j]*Rqj'*Rqj))  *(iGq*L[i]*Rqi')
        lbc1 = lus1 - rho^(p+1)*(T1 + T2 + T3)
   
        q10 = q10 + lbc1*lbc0[i]*(lbc0[j]*res[j])^2*vx[i]
        q11 = q11 + lbc1*lbc0[i]*((lbc0[j]^2)*vx[j]-h*s2)*(res[i]^2)
        
        q6 = q6 + (lbc0[i]^2)*(Rqi*iGq*L[j]*Rqj')^2*(res[j]^2) 
       }
    }
  }
  

  Eq1  = (q1/(N*h))^2
  Eq2  = q2/(N*h)
  Eq3  = q3/(N*h)
  Eq4  = q4/(N*h)
  Eq5  = (q5a/(N*h))*(q5b/(N*h))
  Eq6  = q6/(N*(N-1)*(h^2))
  
  Eq7  = (q7a/(N*h))*(q7b/(N*h))*(q7c/(N*h))

  Eq8  = q8/(N*h)
  Eq9  = q9/(N*h)
  Eq10 = q10/(N*(N-1)*(h^2))
  Eq11 = q11/(N*(N-1)*(h^2))
  Eq12 = q12/(N*h)
  
  z  = 1.959964
  pz = 0.05844507
  
  q1bc =   pz*(
       Eq1*((z^3)/3+7*z/4+s2*z*((z^2)-3)/4)/(s2^3)
    +  Eq2*(-z*((z^2)-3)/2)/s2
    +  Eq3*(z*((z^2)-3)/8)/(s2^2) 
    -  Eq4*(z*((z^2)-1)/2)/s2 
    -  Eq5*(z*((z^2)-1))/(s2^2) 
    +  Eq6*(z*((z^2)-1)/4)/s2 
    +  Eq7*(z*((z^2)-1)/2)/(s2^2) 
    +  Eq8*(-z*((z^2)-3)/24)/(s2^2) 
    +  Eq9*(z*((z^2)-1)/4)/(s2^2) 
    +  Eq10*(z*((z^2)-3))/(s2^2) 
    +  Eq11*(-z)/(s2^2) 
    +  Eq12*(-z*((z^2)+1)/8)/(s2^2)
    )
    
    q2bc = -pz*z/(2*s2)
    
    Eq3a = q3a/(N*h)
    q3bc = pz*Eq3a/(s2^2)*(z^3)/3
    
    q1rbc = 2*q1bc/pz
    q2rbc = 2*q2bc/pz
    q3rbc = 2*q3bc/pz
	
	return(q1rbc,q2rbc,q3rbc)	
}
mata mosave nprobust_lp_q(), replace
end


****************************************************
*capture mata mata drop nprobust_lp_cer()
*mata
*void nprobust_lp_cer(real scalar todo, real scalar H,  real vector q_rbc, real scalar eta, real scalar p, val, grad, hess) 
*{	
*		val = (H^(-1)*q_rbc[1] + H^(1+2*(p+3))*eta^2*q_rbc[2] + H^(p+3)*eta*q_rbc[3])^2
*}
*mata mosave nprobust_lp_cer(), replace
*end

capture mata mata drop nprobust_lp_cer()
mata
void nprobust_lp_cer(real scalar todo, real scalar H,  real vector E, real scalar q2_rbc, real scalar p, real scalar N, val, grad, hess) 
{	
	val = abs(H^(-1)*E[1] + N*H^(2*p+5)*q2_rbc*0.05844094*(E[2] + H*E[3])^2 + H^(p+2)*(E[4]+H*E[5]))
}
mata mosave nprobust_lp_cer(), replace
end
	
capture mata mata drop nprobust_lp_mse()
mata
void nprobust_lp_mse(real scalar todo, real scalar h,  real scalar V, real vector B, real scalar p, real scalar v, real scalar n, val, grad, hess) 
{	
	val = abs(h^(2*p+2-v)*(B[1]+h*B[2])^2 + V/(n*h^(1+2*v)))
}
mata mosave nprobust_lp_mse(), replace
end


*** MSE-DPI
capture mata mata drop nprobust_lp_mse_dpi()
mata
real matrix nprobust_lp_mse_dpi(real vector Y, real vector X, real vector C, real scalar eval,  real scalar p, real scalar q, real scalar deriv, real scalar even, string kernel, real scalar c_bw, real scalar bwcheck, real scalar bwregul, string vce, real scalar nnmatch, real vector dups, real vector dupsid, string interior) 
{		
range = max(X)-min(X)
N = length(X)
bw_max = max((abs(eval-min(X)),abs(eval-max(X))))
c_bw = min((c_bw, bw_max))

		if (bwcheck != 0) {
			bw_min = sort(abs(X:-eval), 1)[bwcheck] + 1e-8
			c_bw = max((c_bw, bw_min))
		}
		
    C_d1 = nprobust_lp_bw(Y, X, C, eval, q+1, q+1, q+2, c_bw, range, range, 0, vce, nnmatch, kernel, dups, dupsid)

    if (even==0 | interior!="") {
		bw_mp2  = C_d1[8]
	} else {
		B_d = C_d1[1], C_d1[2]
		V_d  = C_d1[3]
     	S = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_mse())
		optimize_init_params(S, 1e-8)
		optimize_init_argument(S, 1, V_d)
		optimize_init_argument(S, 2, B_d)
		optimize_init_argument(S, 3, q+1)
		optimize_init_argument(S, 4, q+1)
		optimize_init_argument(S, 5, N)
		optimize_init_which(S, "min")	    
		optimize_init_evaluatortype(S, "d0")
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_maxiter(S, 20)
		(void) optimize(S)
		bw_mp2 = optimize_result_params(S)
	}
	
	C_d2 = nprobust_lp_bw(Y, X, C, eval, q+2, q+2, q+3, c_bw, range, range, 0, vce, nnmatch, kernel, dups, dupsid)
    if (even==0 | interior!="") {
		bw_mp3  = C_d2[8]
	} else {
		B_d = C_d2[1], C_d2[2]
		V_d  = C_d2[3]
		S = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_mse())
		optimize_init_params(S, 1e-8)
		optimize_init_argument(S, 1, V_d)
		optimize_init_argument(S, 2, B_d)
		optimize_init_argument(S, 3, q+2)
		optimize_init_argument(S, 4, q+2)
		optimize_init_argument(S, 5, N)
		optimize_init_which(S, "min")	    
		optimize_init_evaluatortype(S, "d0")
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_maxiter(S, 20)
		(void) optimize(S)
		bw_mp3 = optimize_result_params(S)
	}
	
	bw_mp2 = min((bw_mp2, bw_max))
	bw_mp3 = min((bw_mp3, bw_max))	
		
	  	if (bwcheck != 0) {
			bw_mp2 = max((bw_mp2, bw_min))
			bw_mp3 = max((bw_mp3, bw_min))
		}		
		
      C_b  = nprobust_lp_bw(Y, X, C, eval, q, p+1, q+1, c_bw, bw_mp2, bw_mp3, bwregul, vce, nnmatch, kernel, dups, dupsid)
      if (even==0 | interior!="") {
		b_mse_dpi = C_b[8]     
	  } else {
	  	B_b = C_b[1], C_b[2]
		V_b  = C_b[3]
        S = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_mse())
		optimize_init_params(S, 1e-8)
		optimize_init_argument(S, 1, V_b)
		optimize_init_argument(S, 2, B_b)
		optimize_init_argument(S, 3, q)
		optimize_init_argument(S, 4, p+1)
		optimize_init_argument(S, 5, N)
		optimize_init_which(S, "min")	    
		optimize_init_evaluatortype(S, "d0")
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_maxiter(S, 20)
		(void) optimize(S)
		b_mse_dpi = optimize_result_params(S)
	  }
	
	b_mse_dpi = min((b_mse_dpi, bw_max))
		
	   	if (bwcheck != 0) {
			b_mse_dpi = max((b_mse_dpi, bw_min))
		}
		
		bw_mp1 = b_mse_dpi
		
      C_h  = nprobust_lp_bw(Y, X, C, eval, p, deriv, q, c_bw, bw_mp1, bw_mp2, bwregul, vce, nnmatch, kernel, dups, dupsid)
	  if (even==0 | interior!="") {
		h_mse_dpi = C_h[8]     
	  } else {
	  	B_h = C_h[1], C_h[2]
		V_h  = C_h[3]
        S = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_mse())
		optimize_init_params(S, 1e-8)
		optimize_init_argument(S, 1, V_h)
		optimize_init_argument(S, 2, B_h)
		optimize_init_argument(S, 3, p)
		optimize_init_argument(S, 4, deriv)
		optimize_init_argument(S, 5, N)
		optimize_init_which(S, "min")	    
		optimize_init_evaluatortype(S, "d0")
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_maxiter(S, 20)
		(void) optimize(S)
		h_mse_dpi = optimize_result_params(S)
		}
 
	h_mse_dpi = min((h_mse_dpi, bw_max))
 
	  if (bwcheck != 0) {
			h_mse_dpi = max((h_mse_dpi, bw_min))
	  }
		
	  if (even==0 | interior!="") {
			V_h = C_h[6]*C_h[3]
			B_h = C_h[5]*C_h[1]^2
			V_b = C_b[6]*C_b[3]
			B_b = C_b[5]*C_b[1]^2
	  } else {
           V_h = V_h
		   B_h = (B_h[1]+h_mse_dpi*B_h[2])^2
		   V_b = V_b
		   B_b = (B_b[1]+b_mse_dpi*B_b[2])^2
	   }    
	   
		return(h_mse_dpi, b_mse_dpi, V_h, B_h, V_b, B_b)
}
mata mosave nprobust_lp_mse_dpi(), replace
end


capture mata mata drop nprobust_lp_coef() 
mata
real matrix nprobust_lp_coef(real matrix Y, real matrix X, real scalar eval, real scalar deriv, real scalar p, real scalar h, string kernel) 
{
		w_h = nprobust_lp_kweight(X,eval,h,kernel);	
		ind = selectindex(w_h:> 0)
		eN = length(ind)
		eY  = Y[ind];eX  = X[ind];
		W_h = w_h[ind]
		R_p = J(eN,(p+1),.)
			for (j=1; j<=(p+1); j++)  {
				R_p[.,j] = (eX:-eval):^(j-1)
			}		
		invG_p  = cholinv(quadcross(R_p,W_h,R_p))
		beta_p = invG_p*quadcross(R_p:*W_h,eY)
		tau = factorial(deriv)*beta_p[(deriv+1),1]
		return(tau)
}
mata mosave nprobust_lp_coef(), replace
end
	
	


	
	

/*
*** MSE-ROT
capture mata mata drop nprobust_lp_mse_rot()
mata
real matrix nprobust_lp_mse_rot(real vector Y, real vector X, real scalar eval,  real scalar p, real scalar deriv, string kernel, real scalar c_bw, real scalar even) 
{		
		N = length(X)
		n_h1 = sum(abs(X:-eval):<=c_bw)
        f0_pilot=n_h1/(N*c_bw)
		  C.c         = 2.34
  if (kernel=="uni") C.c         = 1.843
  if (kernel=="tri") C.c         = 2.576
	    k = p+3
	  
      rk = J(N,(k+1),.) 
	  for (j=1; j<=(k+1); j++) {
	    rk[.,j] = X:^(j-1)
	  }	
	  iGp  = cholinv(quadcross(rk,rk))
	  gamma_p = iGp*quadcross(rk,Y)
	  s2_hat = mean((rk*gamma_p:-Y):^2)
      
 
	    mp1 = gamma_p[p+1]*factorial(p+1) + gamma_p[p+2]*factorial(p+2)*eval[e] + gamma_p[p+3]*factorial(p+3)*eval[e]^2/2      
             
        h_mse_rot = C_c*((s2_hat/f_hat)/(N*mp1^2))^(1/(2*p+3))
        
		if ("`bwcheck'" != "0") {
			
			h_mse_rot = max((h_mse_rot, bw_min))
		}
		
		B1 = C1*m.p.1/factorial(p+1)
		B2 = C1*m.p.2/factorial(p+2)
		V  = C2*s2_hat/f_hat

	  if (even==FALSE) {
		h.mse.rot = bw$bw
		B = B1
	  } else {
		h.bw.fun  = function(H) {abs(H^(2*p+2-2*deriv)*(B1 + H*B2)^2 + V/(N*H^(1+2*deriv)))}
		h.mse.rot = optimize(h.bw.fun, interval=c(.Machine$double.eps, range))$minimum
		B = B1 + h.mse.rot*B2
	  }
	  
  out = list(h=h.mse.rot, V=V, B=B)
  
        return = h_mse_rot, s2_p/f0_pilot, mp1^2
		
		
}
mata mosave nprobust_lp_mse_rot(), replace
end





/*
capture mata mata drop nprobust_bw_fun()
mata
real matrix nprobust_bw_fun(real scalar V, real scalar B, real scalar v, real scalar r, real scalar N)
{
((1+2*r)*V/(2*v*N*B^2))^1/(1+2*v+2*r)
}
mata mosave nprobust_bw_fun(), replace
end


capture mata mata drop nprobust_lp_res()
mata
real matrix nprobust_lp_res(real matrix X, real matrix y, real matrix m, real matrix hii, string vce, real scalar matches, dups, dupsid, real scalar d)
{
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

bw.lp.fun = function(V,B,p,v,N,kernel) {
  
  k.fun = function(u){
    if (kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
    if (kernel=="uni") w =          0.5*(abs(u)<=1)
    if (kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
    return(w)  
  }
  
  C1.h = C1.fun(p0=p,v=v, K=k.fun)  
  C2.h = C2.fun(p0=p,v=v, K=k.fun)
  bw = ( ((2*v+1)*C2.h*V) / (2*(p+1-v)*C1.h^2*B*N) )^(1/(2*p+3))
  return(bw)
}
}
mata mosave nprobust_lp_kweight(), replace
end


capture mata mata drop nprobust_kd_K()
mata
real matrix nprobust_kd_K(real matrix u, string kernel)
{
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = 0.75:*(1:-u:^2):*(abs(u):<=1)
		mK = 0.5
		vK = 0.5
	}
	else {
		kx = 0.5:*(abs(u):<=1)
		mK = 0.5
		vK = 0.5
	}
return(kx)
return(mK)	
}
mata mosave nprobust_kd_K(), replace
end
capture mata mata drop nprobust_kd_L()
mata
real matrix nprobust_kd_L(real matrix u, string kernel)
{
	if (kernel=="epanechnikov" | kernel=="epa") {
		kx = (abs(u):<=1):*(105/16):*(6:*u:^2-5:*u:^4:-1)
		mK = 0.5
		vK = 0.5
	}
	else {
		kx = (abs(u):<=1):*15:*(3:*u:^2:-1)/4
		mK = 0.5
		vK = 0.5
	}
return(w,mK,vK)	
}
mata mosave nprobust_kd_L(), replace
end
*/




