*!version 1.0.0  2026-05-17
  
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
real matrix nprobust_lp_bw(real matrix Y, real matrix X, real matrix C, real scalar c, real scalar o, real scalar nu, real scalar o_B, real scalar h_V, real scalar h_B1, real scalar h_B2, real scalar scale, string vce, real scalar nnmatch, string kernel, dups, dupsid, real matrix wts)
{
	dC = eC = indC = 0
	w = nprobust_lp_kweight(X, c, h_V, kernel) :* wts
	ind_V = selectindex(w:> 0); eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
	n_V = length(ind_V)
	D_V = eY
	// Q1: Vandermonde via successive multiplication (matches rdrobust).
	R_V = J(n_V,o+1,1)
	if (o >= 1) {
		u_V = eX :- c
		R_V[.,2] = u_V
		for (j=3; j<=(o+1); j++) R_V[.,j] = R_V[.,j-1] :* u_V
	}
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
			// Vectorized: diag(R (R' W R)^{-1} R' W) = rowsum((R invG) :* (R W))
			hii = rowsum((R_V*invG_V) :* (R_V:*eW))
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
		
	w = nprobust_lp_kweight(X, c, h_B1, kernel) :* wts
	ind = selectindex(w:> 0)
	n_B = length(ind)
	eY = Y[ind];eX = X[ind];eW = w[ind]

	if (dC==1) {
		eC   = C[ind]
		indC = order(eC,1)
	}

	// Q1: Vandermonde via successive multiplication.
	R_B1 = J(n_B,o_B+1,1)
	if (o_B >= 0) {  // o_B+1 >= 1 always
		u_B1 = eX :- c
		if (o_B+1 >= 2) R_B1[.,2] = u_B1
		for (j=3; j<=(o_B+1); j++) R_B1[.,j] = R_B1[.,j-1] :* u_B1
	}
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
				hii = rowsum((R_B1*invG_B1) :* (R_B1:*eW))
			}
		}	
		res_B = nprobust_lp_res(eX, eY, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B, o_B+1)		        
		V_B = (invG_B1*nprobust_lp_vce(R_B1:*eW, res_B, eC, indC)*invG_B1)[o+2,o+2]
		BWreg = 3*BConst1^2*V_B
	}
	
	w = nprobust_lp_kweight(X, c, h_B2, kernel) :* wts
	ind = selectindex(w:> 0)
	n_B = length(ind)
	eY = Y[ind];eX = X[ind];eW = w[ind]

	// Q1: Vandermonde via successive multiplication.
	R_B2 = J(n_B,o_B+2,1)
	u_B2 = eX :- c
	R_B2[.,2] = u_B2
	for (j=3; j<=(o_B+2); j++) R_B2[.,j] = R_B2[.,j-1] :* u_B2
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
** Cluster-robust meat for CR0, CR1, CR2, CR3.
** X:       N x k design (standardized for classical; Q for bias-corrected).
** r:       length-N residuals (std sqrt(W)*raw for classical; raw for BC).
** C_o:     cluster IDs for the reduced sample, already sorted.
** ind:     order() index that maps (C, RX, res) into the sorted frame.
** invG:    k x k inverse Gram for the same design.
** cr_type: one of "CR0", "CR1", "CR2", "CR3".
capture mata mata drop nprobust_lp_cluster_meat()
mata
real matrix nprobust_lp_cluster_meat(real matrix X, real matrix r, real matrix C, real matrix ind, real matrix invG, string cr_type, real scalar k_override)
{
	// k_override: when > 0, overrides the (N-1)/(N-k) df correction's k for
	// CR1. Used when X is a "score-like" matrix (e.g. Q_q) whose ncol is
	// smaller than the effective number of regressors that produced r.
	// Pass 0 to use cols(X) as before.
	C_o   = C[ind]
	X_o   = X[ind,]
	r_o   = r[ind]
	info  = panelsetup(C_o,1)
	G     = rows(info)
	N     = length(r_o)
	k     = cols(X_o)
	k_df  = (k_override > 0 ? k_override : k)
	M     = J(k,k,0)

	for (g=1; g<=G; g++) {
		Xg = panelsubmatrix(X_o, g, info)
		rg = panelsubmatrix(r_o, g, info)

		if (cr_type=="CR0" | cr_type=="CR1") {
			r_adj = rg
		}
		else {
			H_gg = Xg*invG*Xg'
			H_gg = 0.5*(H_gg + H_gg')
			I_H  = I(rows(Xg)) - H_gg
			if (cr_type=="CR2") {
				vals = vecs = .
				symeigensystem(I_H, vecs, vals)
				// symeigensystem returns eigenvalues as a row vector;
				// convert to column and floor at 1e-12 for stability.
				vals = vals'
				vals = vals :* (vals :>= 1e-12) :+ 1e-12 :* (vals :< 1e-12)
				r_adj = vecs*((vecs'*rg):/sqrt(vals))
			}
			else { // CR3
				r_adj = lusolve(I_H, rg)
				if (hasmissing(r_adj)) {
					vals = vecs = .
					symeigensystem(I_H, vecs, vals)
					vals = vals'
					vals = vals :* (vals :>= 1e-12) :+ 1e-12 :* (vals :< 1e-12)
					r_adj = vecs*((vecs'*rg):/vals)
				}
			}
		}
		score = Xg'*r_adj
		M = M + score*score'
	}

	mult = 1
	if (cr_type=="CR1") mult = ((N-1)/(N-k_df))*(G/(G-1))
	if (cr_type=="CR3") mult = (G-1)/G
	return(mult*M)
}
mata mosave nprobust_lp_cluster_meat(), replace
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
	N   = length(y)
	rho = h / b
	Wp  = K:/h
	Wq  = L:/b
	Xh  = (x:-c):/h
	Xb  = (x:-c):/b

	Rp = J(N, p+1, 1)
	if (p >= 1) {
		Rp[.,2] = Xh
		for (i=3; i<=(p+1); i++) Rp[.,i] = Rp[.,i-1] :* Xh
	}
	Rq = J(N, q+1, 1)
	if (q >= 1) {
		Rq[.,2] = Xb
		for (i=3; i<=(q+1); i++) Rq[.,i] = Rq[.,i-1] :* Xb
	}

	Lp1 = quadcross(Rp, Wp:*Xh:^(p+1)) / N
	Gp  = quadcross(Rp, Wp, Rp) / N
	Gq  = quadcross(Rq, Wq, Rq) / N
	iGp = luinv(Gp)
	iGq = luinv(Gq)

	ep1 = J(q+1,1,0); ep1[p+2] = 1
	e0  = J(p+1,1,0); e0[deriv+1] = factorial(deriv)

	lus0     = (e0'*iGp) * (Rp:*K)'
	scale_bc = (e0'*iGp) * Lp1
	tail_bc  = (ep1'*iGq) * (Rq:*L)'
	lbc0     = lus0' - rho^(p+1) * scale_bc * tail_bc'

	vx = res:^2
	s2 = sum((lbc0:^2):*vx) / (N*h)

	Krrp   = quadcross(Rp:*K, Rp)
	Lrrq   = quadcross(Rq:*L, Rq)
	Krxip  = quadcross(Rp, K:*Xh:^(p+1))
	sumKRp = quadcross(Rp, K)
	Sxp1   = sum(Xh:^(p+1))
	Krxp   = sumKRp:*Sxp1 - Krxip

	EKrrp      = Krrp / N
	EKrxip_vec = Krxip / N
	EKrxp_vec  = Krxp / (N*(N-1))
	ELrrq      = Lrrq / N

	a_row = factorial(deriv) * iGp[deriv+1,.]
	u_a   = a_row
	v_a   = a_row * EKrrp * iGp

	RpiGp  = Rp * iGp
	RqiGq  = Rq * iGq
	quadRp = rowsum(RpiGp:*Rp)
	quadRq = rowsum(RqiGq:*Rq)

	Rp_va = Rp * v_a'
	Rp_ua = Rp * u_a'

	q1  = sum((lbc0:*res):^3)
	q3  = sum((lbc0:^4):*((res:^4) - (vx:^2)))
	q3a = q1
	q8  = sum((lbc0:*res):^4)
	q9  = sum(((lbc0:^2):*vx :- h*s2) :* ((lbc0:*res):^2))
	q12 = sum(((lbc0:^2):*vx :- h*s2):^2)
	q4  = sum((lbc0:^2):*L:*quadRq:*vx)

	q5a = colsum(RqiGq :* ((lbc0:^3):*vx))
	q5b = colsum(Rq    :* (lbc0:*vx:*L))'
	q7a = colsum(RqiGq :* (lbc0:*vx:*L))
	q7b = quadcross(Rq:*lbc0, Rq:*lbc0) * iGq
	q7c = colsum(Rq    :* (lbc0:*vx:*L))'

	lus1_diag = K:*Rp_va - (K:^2):*Rp_ua:*quadRp

	C1_scalar   = (a_row * EKrrp * iGp * Lp1)[1,1]
	v_iGq_ep1   = iGq * ep1
	C2_vec      = Rq * v_iGq_ep1
	v_lp        = iGp * Lp1
	D2_vec      = Rp * v_lp
	T1_diag     = L:*C2_vec:*(C1_scalar :- K:*Rp_ua:*D2_vec)
	dot_a_krxip = (a_row * EKrxip_vec)[1,1]
	T2_diag     = L:*C2_vec:*(K:*Xh:^(p+1):*Rp_ua :- dot_a_krxip)
	a_Lp1       = (a_row * Lp1)[1,1]
	u_T3        = a_Lp1 * (ep1' * iGq)
	u_T3_ELrrq_iGq = u_T3 * ELrrq * iGq
	Rq_uT3          = Rq * u_T3'
	Rq_uT3_ELrrq    = Rq * u_T3_ELrrq_iGq'
	T3_diag         = L:*Rq_uT3_ELrrq - (L:^2):*Rq_uT3:*quadRq
	lbc1_diag       = lus1_diag - rho^(p+1) * (T1_diag + T2_diag + T3_diag)
	q2              = sum(lbc1_diag:*lbc0:*vx)

	a_q6 = lbc0:^2
	b_q6 = (L:^2):*vx
	A_q6 = quadcross(Rq, a_q6, Rq)
	B_q6 = quadcross(Rq, b_q6, Rq)
	total_q6 = sum(A_q6 :* (iGq * B_q6 * iGq))
	diag_q6  = sum(a_q6:*b_q6:*(quadRq:^2))
	q6 = total_q6 - diag_q6

	dot_a_krxp = (a_row * EKrxp_vec)[1,1]
	Xhpp1      = Xh:^(p+1)
	alpha      = K:*Rp_va
	gamma      = K:*Rp_ua
	delta      = L:*C2_vec
	epsilon    = K:*Rp_ua:*D2_vec
	zeta       = L:*Rq_uT3_ELrrq
	eta_v      = L:*Rq_uT3

	f = lbc0:*vx
	g = (lbc0:^2):*vx
	sum_g   = sum(g)
	sum_df  = sum(delta:*f)
	fg      = f:*g
	sum_dfg = sum(delta:*fg)
	St_A  = sum(alpha:*f) * sum_g
	Sd_A  = sum(alpha:*fg)
	uB    = quadcross(Rp, K:*f)
	vB    = quadcross(Rp, gamma:*g)
	St_B = (uB' * iGp * vB)[1,1]
	Sd_B = sum(K:*gamma:*quadRp:*fg)
	St_T1 = C1_scalar*sum_df*sum_g - sum_df*sum(epsilon:*g)
	Sd_T1 = sum(delta:*fg:*(C1_scalar :- epsilon))
	St_T2 = sum(delta:*Xhpp1:*f)*sum(gamma:*g) - dot_a_krxp*sum_df*sum_g
	Sd_T2 = sum(delta:*Xhpp1:*gamma:*fg) - dot_a_krxp*sum_dfg
	uT3x  = quadcross(Rq, L:*f)
	vT3x  = quadcross(Rq, eta_v:*g)
	St_T3 = sum(zeta:*f)*sum_g - (uT3x' * iGq * vT3x)[1,1]
	Sd_T3 = sum(zeta:*fg) - sum(L:*eta_v:*quadRq:*fg)
	q10 = (St_A - St_B - rho^(p+1)*(St_T1 + St_T2 + St_T3)) - ///
	      (Sd_A - Sd_B - rho^(p+1)*(Sd_T1 + Sd_T2 + Sd_T3))

	g = (lbc0:^2):*vx :- h*s2
	sum_g   = sum(g)
	sum_df  = sum(delta:*f)
	fg      = f:*g
	sum_dfg = sum(delta:*fg)
	St_A  = sum(alpha:*f) * sum_g
	Sd_A  = sum(alpha:*fg)
	uB    = quadcross(Rp, K:*f)
	vB    = quadcross(Rp, gamma:*g)
	St_B = (uB' * iGp * vB)[1,1]
	Sd_B = sum(K:*gamma:*quadRp:*fg)
	St_T1 = C1_scalar*sum_df*sum_g - sum_df*sum(epsilon:*g)
	Sd_T1 = sum(delta:*fg:*(C1_scalar :- epsilon))
	St_T2 = sum(delta:*Xhpp1:*f)*sum(gamma:*g) - dot_a_krxp*sum_df*sum_g
	Sd_T2 = sum(delta:*Xhpp1:*gamma:*fg) - dot_a_krxp*sum_dfg
	uT3x  = quadcross(Rq, L:*f)
	vT3x  = quadcross(Rq, eta_v:*g)
	St_T3 = sum(zeta:*f)*sum_g - (uT3x' * iGq * vT3x)[1,1]
	Sd_T3 = sum(zeta:*fg) - sum(L:*eta_v:*quadRq:*fg)
	q11 = (St_A - St_B - rho^(p+1)*(St_T1 + St_T2 + St_T3)) - ///
	      (Sd_A - Sd_B - rho^(p+1)*(Sd_T1 + Sd_T2 + Sd_T3))

	Eq1  = (q1/(N*h))^2
	Eq2  = q2/(N*h)
	Eq3  = q3/(N*h)
	Eq4  = q4/(N*h)
	Eq5  = (q5a*q5b)/(N*h)^2
	Eq6  = q6/(N*(N-1)*(h^2))
	Eq7  = (q7a*q7b*q7c)/(N*h)^3
	Eq8  = q8/(N*h)
	Eq9  = q9/(N*h)
	Eq10 = q10/(N*(N-1)*(h^2))
	Eq11 = q11/(N*(N-1)*(h^2))
	Eq12 = q12/(N*h)

	z  = 1.959964
	pz = 0.05844507

	q1bc = pz*(
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
	// Bug fix: leading term is E[1]/(N*H), not E[1]/H (mirrors R).
	val = abs(E[1]/(N*H) + N*H^(2*p+5)*q2_rbc*normalden(1.96)*(E[2] + H*E[3])^2 + H^(p+2)*(E[4]+H*E[5]))
}
mata mosave nprobust_lp_cer(), replace
end
	
capture mata mata drop nprobust_lp_mse()
mata
void nprobust_lp_mse(real scalar todo, real scalar h,  real scalar V, real vector B, real scalar p, real scalar v, real scalar n, val, grad, hess)
{
	// R uses H^(2p+2-2*deriv) (npfunctions.R:769, 779, 800, 815, 887);
	// the prior exponent 2*p+2-v was off by a factor 2 on deriv. This
	// only affects the optimize branch (even & !interior cases).
	val = abs(h^(2*p+2-2*v)*(B[1]+h*B[2])^2 + V/(n*h^(1+2*v)))
}
mata mosave nprobust_lp_mse(), replace
end

// IMSE-RoT optimizer: argmin H of |H^(2p+2-2v) mean(B^2) + mean(V)/(N H^(1+2v))|.
// Mirrors the optimize branch of R lpbwselect.imse.rot (npfunctions.R:720-722).
capture mata mata drop nprobust_lp_mse_imse()
mata
void nprobust_lp_mse_imse(real scalar todo, real scalar h, real scalar mV, real scalar mBsq, real scalar p, real scalar v, real scalar n, val, grad, hess)
{
	val = abs(h^(2*p+2-2*v)*mBsq + mV/(n*h^(1+2*v)))
}
mata mosave nprobust_lp_mse_imse(), replace
end

// Golden-section search on [a,b] for the per-point MSE-RoT bandwidth.
// Mirrors R optimize() with interval=c(.Machine$double.eps, range) used
// in lpbwselect.mse.rot (npfunctions.R:887-888). Returns argmin H of
//   abs(H^(2p+2-2v)*(B1 + H*B2)^2 + V/(N*H^(1+2v))).
capture mata mata drop nprobust_lp_brent_mserot()
mata
real scalar nprobust_lp_brent_mserot(real scalar V, real vector B, real scalar p, real scalar v, real scalar N, real scalar lo, real scalar hi)
{
	real scalar cg, a, b, x, w, vv, u, d, e, fx, fw, fv, fu
	real scalar xm, tol, tol1, tol2, pp, qq, rr, iter, eps
	cg  = (3-sqrt(5))/2
	eps = sqrt(2.220446049250313e-16)
	tol = sqrt(eps)
	a = lo; b = hi
	x = a + cg*(b-a)
	w = x
	vv = x
	d = e = 0
	fx = abs(x^(2*p+2-2*v)*(B[1]+x*B[2])^2 + V/(N*x^(1+2*v)))
	fw = fv = fx
	for (iter=1; iter<=500; iter++) {
		xm = 0.5*(a+b)
		tol1 = eps*abs(x) + tol/3
		tol2 = 2*tol1
		if (abs(x-xm) <= tol2 - 0.5*(b-a)) break
		pp = qq = rr = 0
		if (abs(e) > tol1) {
			rr = (x-w)*(fx-fv)
			qq = (x-vv)*(fx-fw)
			pp = (x-vv)*qq - (x-w)*rr
			qq = 2*(qq-rr)
			if (qq > 0) pp = -pp
			else        qq = -qq
			rr = e
			e = d
		}
		if (abs(pp) >= abs(0.5*qq*rr) | pp <= qq*(a-x) | pp >= qq*(b-x)) {
			e = (x < xm ? b-x : a-x)
			d = cg*e
		}
		else {
			d = pp/qq
			u = x + d
			if ((u-a) < tol2 | (b-u) < tol2) d = (xm >= x ? tol1 : -tol1)
		}
		u = (abs(d) >= tol1 ? x+d : x + (d > 0 ? tol1 : -tol1))
		fu = abs(u^(2*p+2-2*v)*(B[1]+u*B[2])^2 + V/(N*u^(1+2*v)))
		if (fu <= fx) {
			if (u < x) b = x
			else       a = x
			vv = w; fv = fw
			w = x;  fw = fx
			x = u;  fx = fu
		}
		else {
			if (u < x) a = u
			else       b = u
			if (fu <= fw | w == x) {
				vv = w; fv = fw
				w = u;  fw = fu
			}
			else if (fu <= fv | vv == x | vv == w) {
				vv = u; fv = fu
			}
		}
	}
	return(x)
}
mata mosave nprobust_lp_brent_mserot(), replace
end

// Golden-section search for CE-DPI interior=FALSE branch. Mirrors R
// optimize() in lpbwselect.ce.dpi npfunctions.R:656-657:
//   abs(E1/(N*H) + N*H^(2p+5)*q2*phi*(E2 + H*E3 + bwregul*Reg)^2
//                + H^(p+2)*(E4 + H*E5 + bwregul*Reg))
// (Reg=0 since R hard-codes bwregul=0; we mirror that.)
capture mata mata drop nprobust_lp_brent_cedpi()
mata
real scalar nprobust_lp_brent_cedpi(real vector E, real scalar q2, real scalar p, real scalar N, real scalar lo, real scalar hi)
{
	real scalar cg, a, b, x, w, vv, u, d, e, fx, fw, fv, fu
	real scalar xm, tol, tol1, tol2, pp, qq, rr, iter, eps, phiz
	cg = (3-sqrt(5))/2
	eps = sqrt(2.220446049250313e-16)
	tol = sqrt(eps)
	phiz = normalden(1.96)  // R: dnorm(1.96) = 0.0584409443334515 (full precision)
	a = lo; b = hi
	x = a + cg*(b-a)
	w = x
	vv = x
	d = e = 0
	fx = abs(E[1]/(N*x) + N*x^(2*p+5)*q2*phiz*(E[2]+x*E[3])^2 + x^(p+2)*(E[4]+x*E[5]))
	fw = fv = fx
	for (iter=1; iter<=500; iter++) {
		xm = 0.5*(a+b)
		tol1 = eps*abs(x) + tol/3
		tol2 = 2*tol1
		if (abs(x-xm) <= tol2 - 0.5*(b-a)) break
		pp = qq = rr = 0
		if (abs(e) > tol1) {
			rr = (x-w)*(fx-fv)
			qq = (x-vv)*(fx-fw)
			pp = (x-vv)*qq - (x-w)*rr
			qq = 2*(qq-rr)
			if (qq > 0) pp = -pp
			else        qq = -qq
			rr = e
			e = d
		}
		if (abs(pp) >= abs(0.5*qq*rr) | pp <= qq*(a-x) | pp >= qq*(b-x)) {
			e = (x < xm ? b-x : a-x)
			d = cg*e
		}
		else {
			d = pp/qq
			u = x + d
			if ((u-a) < tol2 | (b-u) < tol2) d = (xm >= x ? tol1 : -tol1)
		}
		u = (abs(d) >= tol1 ? x+d : x + (d > 0 ? tol1 : -tol1))
		fu = abs(E[1]/(N*u) + N*u^(2*p+5)*q2*phiz*(E[2]+u*E[3])^2 + u^(p+2)*(E[4]+u*E[5]))
		if (fu <= fx) {
			if (u < x) b = x
			else       a = x
			vv = w; fv = fw
			w = x;  fw = fx
			x = u;  fx = fu
		}
		else {
			if (u < x) a = u
			else       b = u
			if (fu <= fw | w == x) {
				vv = w; fv = fw
				w = u;  fw = fu
			}
			else if (fu <= fv | vv == x | vv == w) {
				vv = u; fv = fu
			}
		}
	}
	return(x)
}
mata mosave nprobust_lp_brent_cedpi(), replace
end


*** MSE-DPI
capture mata mata drop nprobust_lp_mse_dpi()
mata
real matrix nprobust_lp_mse_dpi(real vector Y, real vector X, real vector C, real scalar eval,  real scalar p, real scalar q, real scalar deriv, real scalar even, string kernel, real scalar c_bw, real scalar bwcheck, real scalar bwregul, string vce, real scalar nnmatch, real vector dups, real vector dupsid, string interior, real matrix wts)
{
range = max(X)-min(X)
N = length(X)
bw_max = max((abs(eval-min(X)),abs(eval-max(X))))
c_bw = min((c_bw, bw_max))

		if (bwcheck != 0) {
			bw_min = sort(abs(X:-eval), 1)[bwcheck]
			c_bw = max((c_bw, bw_min))
		}

    C_d1 = nprobust_lp_bw(Y, X, C, eval, q+1, q+1, q+2, c_bw, range, range, 0, vce, nnmatch, kernel, dups, dupsid, wts)
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
	
	C_d2 = nprobust_lp_bw(Y, X, C, eval, q+2, q+2, q+3, c_bw, range, range, 0, vce, nnmatch, kernel, dups, dupsid, wts)
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
		
      C_b  = nprobust_lp_bw(Y, X, C, eval, q, p+1, q+1, c_bw, bw_mp2, bw_mp3, bwregul, vce, nnmatch, kernel, dups, dupsid, wts)
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
		
      C_h  = nprobust_lp_bw(Y, X, C, eval, p, deriv, q, c_bw, bw_mp1, bw_mp2, bwregul, vce, nnmatch, kernel, dups, dupsid, wts)
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


// Compute the kernel constants C1 and C2 used by the local-polynomial
// rule-of-thumb bandwidth (mirrors R lp.bw.fun in npfunctions.R:282-314).
//   C1 = (Sinv NU)[v+1]
//   C2 = (Sinv PSI Sinv)[v+1, v+1]
// where Sinv = inv(GAMMA(p,K)), GAMMA[i,j] = m1(i+j,K),
//       NU[i] = m1(i+p+1,K), PSI[i,j] = m2(i+j,K),
// with m1(k,K) = ∫₀^∞ x^k K(x) dx and m2(k,K) = ∫₀^∞ x^k K(x)² dx.
// Closed forms are used for epa, uni, tri; the half-Gaussian moments
// for gau use Γ((k+1)/2)·2^((k-1)/2)/√π and Γ((k+1)/2)/(4π).
capture mata mata drop nprobust_lp_kmom1()
mata
real scalar nprobust_lp_kmom1(real scalar k, string kernel)
{
	if (kernel=="epa") return(0.75 * (1/(k+1) - 1/(k+3)))
	if (kernel=="uni") return(0.5 / (k+1))
	if (kernel=="tri") return(1/(k+1) - 1/(k+2))
	// gaussian: ∫₀^∞ x^k φ(x) dx = 2^((k-1)/2) Γ((k+1)/2) / √π
	return(exp(((k-1)/2)*ln(2) + lngamma((k+1)/2) - 0.5*ln(pi())))
}
mata mosave nprobust_lp_kmom1(), replace
end

capture mata mata drop nprobust_lp_kmom2()
mata
real scalar nprobust_lp_kmom2(real scalar k, string kernel)
{
	if (kernel=="epa") return(0.5625 * (1/(k+1) - 2/(k+3) + 1/(k+5)))
	if (kernel=="uni") return(0.25 / (k+1))
	if (kernel=="tri") return(1/(k+1) - 2/(k+2) + 1/(k+3))
	// gaussian: ∫₀^∞ x^k φ(x)² dx = Γ((k+1)/2) / (4π)·... full form:
	// φ²(x) = (1/(2π)) exp(-x²); ∫₀^∞ x^k exp(-x²) dx = (1/2) Γ((k+1)/2)
	return(exp(lngamma((k+1)/2) - ln(2) - ln(2*pi())))
}
mata mosave nprobust_lp_kmom2(), replace
end

capture mata mata drop nprobust_lp_constants()
mata
real rowvector nprobust_lp_constants(real scalar p, real scalar v, string kernel)
{
	GAMMA = J(p+1, p+1, .)
	NU    = J(p+1, 1, .)
	PSI   = J(p+1, p+1, .)
	for (i=0; i<=p; i++) {
		NU[i+1] = nprobust_lp_kmom1(i+p+1, kernel)
		for (j=0; j<=p; j++) {
			GAMMA[i+1,j+1] = nprobust_lp_kmom1(i+j, kernel)
			PSI[i+1,j+1]   = nprobust_lp_kmom2(i+j, kernel)
		}
	}
	Sinv = invsym(GAMMA)
	C1   = (Sinv*NU)[v+1]
	C2   = (Sinv*PSI*Sinv)[v+1, v+1]
	return((C1, C2))
}
mata mosave nprobust_lp_constants(), replace
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
		// Q1: Vandermonde via successive multiplication.
		R_p = J(eN,(p+1),1)
		if (p >= 1) {
			u_R = eX :- eval
			R_p[.,2] = u_R
			for (j=3; j<=(p+1); j++) R_p[.,j] = R_p[.,j-1] :* u_R
		}
		invG_p  = cholinv(quadcross(R_p,W_h,R_p))
		beta_p = invG_p*quadcross(R_p:*W_h,eY)
		tau = factorial(deriv)*beta_p[(deriv+1),1]
		return(tau)
}
mata mosave nprobust_lp_coef(), replace
end

capture erase lnprobust.mlib
mata: mata mlib create lnprobust, replace
mata: mata mlib add lnprobust nprobust_w_fun() nprobust_K_fun() nprobust_L_fun() nprobust_K() nprobust_L()
mata: mata mlib add lnprobust nprobust_lp_res() nprobust_lp_kweight() nprobust_lp_bw() nprobust_lp_cluster_meat() nprobust_lp_vce()
mata: mata mlib add lnprobust nprobust_lp_q() nprobust_lp_cer() nprobust_lp_mse() nprobust_lp_mse_imse()
mata: mata mlib add lnprobust nprobust_lp_brent_mserot() nprobust_lp_brent_cedpi() nprobust_lp_mse_dpi()
mata: mata mlib add lnprobust nprobust_lp_kmom1() nprobust_lp_kmom2() nprobust_lp_constants() nprobust_lp_coef()
mata: mata mlib index
	
	


	
	

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

      // Q1: Vandermonde via successive multiplication.
      rk = J(N,(k+1),1)
      if (k >= 1) {
        rk[.,2] = X
        for (j=3; j<=(k+1); j++) rk[.,j] = rk[.,j-1] :* X
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



