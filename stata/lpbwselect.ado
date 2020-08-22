*!version 0.3.2  2020-08-22
 
capture program drop lpbwselect
program define lpbwselect, eclass
	syntax anything [if] [in] [, eval(varname) deriv(real 0) neval(real 0) p(real 1) kernel(string) bwselect(string) separator(integer 5) bwregul(real 1) vce(string) genvars interior bwcheck(real 21)]

	marksample touse

	tokenize "`anything'"
	local y `1'
	local x `2'
	local kernel   = lower("`kernel'")
	local bwselect = lower("`bwselect'") 
	
	******************** Set VCE ***************************
	local nnmatch = 3
	tokenize `vce'	
	local w : word count `vce'
	if `w' == 1 {
		local vce_select `"`1'"'
	}
	if `w' == 2 {
		local vce_select `"`1'"'
		if ("`vce_select'"=="nn")      local nnmatch     `"`2'"'
		if ("`vce_select'"=="cluster" | "`vce_select'"=="nncluster") local clustvar `"`2'"'	
	}
	if `w' == 3 {
		local vce_select `"`1'"'
		local clustvar   `"`2'"'
		local nnmatch    `"`3'"'
		if ("`vce_select'"!="cluster" & "`vce_select'"!="nncluster") di as error  "{err}{cmd:vce()} incorrectly specified"  
	}
	if `w' > 3 {
		di as error  "{err}{cmd:vce()} incorrectly specified"  
		exit 125
	}
	
	local vce_type = "NN"
	if ("`vce_select'"=="hc0")     		 local vce_type = "HC0"
	if ("`vce_select'"=="hc1")      	 local vce_type = "HC1"
	if ("`vce_select'"=="hc2")      	 local vce_type = "HC2"
	if ("`vce_select'"=="hc3")      	 local vce_type = "HC3"
	if ("`vce_select'"=="cluster")  	 local vce_type = "Cluster"
	if ("`vce_select'"=="nncluster") 	 local vce_type = "NNcluster"

	if ("`vce_select'"=="cluster" | "`vce_select'"=="nncluster") local cluster = "cluster"
	if ("`vce_select'"=="cluster")       local vce_select = "hc0"
	if ("`vce_select'"=="nncluster")     local vce_select = "nn"
	if ("`vce_select'"=="")              local vce_select = "nn"
	
	
	if ("`deriv'">"0" & "`p'"=="1") local p = (`deriv'+1)
	local q = `p'+1
		
	********************************************************************************
	**** error check: main variable
	cap confirm numeric variable `y'
	if _rc {
		di as err `"`y' is not a numeric variable"'
		exit 198
	}
	cap confirm numeric variable `x'
	if _rc {
		di as err `"`x' is not a numeric variable"'
		exit 198
	}
	else {
		qui su `x' if `x'!=. & `y'!=. & `touse', d
		local x_min = r(min)
		local x_max = r(max)
		local x_sd = r(sd)
		local n = r(N)
		local x_iq = r(p75)-r(p25)
		local range = `x_max' - `x_min'
		if (`n' == 0) {
			di as err `"`x' has length 0"'
			exit 198
		}
	}

	

	********************************************************************************
	**** error check: eval()
	if ("`eval'" == "") {
		if ("`neval'" == "0") {
			*tempvar temp_grid temp_dup
			*qui duplicates tag `x' if `touse', gen(`temp_dup')
			*qui gen `temp_grid' = `x' if `touse' & `temp_dup'==0
			*qui su `temp_grid'
			*local neval = r(N)
			local neval = 30
			tempvar temp_grid	
			*local nquant = `neval'+1
			*pctile `temp_grid' = `x' if `touse', nq(`nquant')	
			qui range `temp_grid' `x_min' `x_max'  `neval'
		}
		else {
				tempvar temp_grid	
				*local nquant = `neval'+1
				*pctile `temp_grid' = `x' if `touse', nq(`nquant')
				qui range `temp_grid' `x_min' `x_max'  `neval'
				
				*su `temp_grid' in 1				
				*tab `temp_grid'
				*local neval = 19
				*egen `temp_grid' = seq() , from(`x_min') to(`x_max')
				*egen d = seq(), f(10) t(12)
				*tempvar temp_grid 
				*su
				*pctile `temp_grid' = `x' 
				*pctile `temp_grid' = `x' if `touse', nq(10)
			}		
		} 
	else {
		cap confirm numeric variable `eval'
		if _rc {
			di as err `"eval(): `eval' is not a numeric variable"'
			exit 198
		}
		else {
			if ("`x'" == "`eval'") {
				qui count if `eval' != . & `touse'
			}
			else {
				qui count if `eval' != .
			}
			
			local neval = r(N)
			
			if (`neval' == 0) {
				di as err `"eval(): `eval' has length 0"'
				exit 198
			}
		}
		local temp_grid "`eval'"
	}

	if ("`bwselect'" == "imse-dpi" | "`bwselect'" == "imse-rot" | "`bwselect'" == "all") {
		tempvar temp_grid_imse	
		local nquant = 30
		*pctile `temp_grid_imse' = `x' if `touse', nq(`nquant') 
		qui range `temp_grid_imse' `x_min' `x_max'  `nquant'
	}

	
	
	********************************************************************************
	**** error check: bwselect()
	if ("`bwselect'" == "") {
		local bwselect = "mse-dpi"
	} 
	else {
		if ("`bwselect'" != "mse-dpi" & "`bwselect'" != "imse-dpi" & "`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot" & "`bwselect'" != "ce-dpi" & "`bwselect'" != "ce-rot" & "`bwselect'" != "all") {
			di as err `"bwselect(): incorrectly specified: options(mse-dpi, imse-dpi, mse-rot, imse-rot)"'
			exit 198
		}
	}
 

	********************************************************************************
	**** error check: kernel()
	if ("`kernel'" == "") {
		local kernel = "epa"
	} 
	else {
		if ("`kernel'" != "triangular" & "`kernel'" != "tri" & "`kernel'" != "uniform" & "`kernel'" != "uni" & "`kernel'" != "epa" & "`kernel'" != "epanechnikov" & "`kernel'" != "gau" & "`kernel'" != "gaussian") {
			di as err `"kernel(): incorrectly specified: options(triangular, uniform, epanechnikov)"'
			exit 198
		}
	}

	if      ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") local kernel_type = "Epanechnikov"
	else if ("`kernel'"=="uniform"      | "`kernel'"=="uni") local kernel_type = "Uniform"
	else if ("`kernel'"=="gaussian"     | "`kernel'"=="gau") local kernel_type = "Gaussian"
	else                                                     local kernel_type = "Triangular"


	********************************************************************************
	**** error check: separator()
	if (`separator' <= 1) {
		local separator = 1
	}

	********************************************************************************
	**** temporaty varaibles for plotting
	if ("`genvars'" != "") local genvars = "lpbwselect"
	
	
	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
		local kernel_type = "Epanechnikov"
		local C_c = 2.34
	}
	else if ("`kernel'"=="uniform" | "`kernel'"=="uni") {
		local kernel_type = "Uniform"
		local C_c = 1.843
	}
	else if ("`kernel'"=="gaussian" | "`kernel'"=="gau") {
		local kernel_type = "Gaussian"
		local C_c = 1.06
	}
	else {
		local kernel_type = "Triangular"
		local C_c = 2.576
	}
	
	/*
	if ("`vce_select'"=="nn") {
		sort `x', stable
		tempvar dups dupsid 
		*ord
		*gen `ord' = _n
		by `x': gen `dups'   = _N if `touse'
		by `x': gen `dupsid' = _n if `touse'
		*sort `ord'		
	}
	*/
	
	local even = !mod(`p'-`deriv',2)
	if (`even'==1 & "`bwselect'"=="ce-dpi") local bwselect="ce-rot"
  

	************
	*** MATA ***
	************
	
	tempvar temp_touse
	qui gen `temp_touse' = `touse'
	
	mata{
	
	Y = st_data(., "`y'", "`temp_touse'") 
	X = st_data(., "`x'", "`temp_touse'") 
	C = 0
	if ("`cluster'"~="") {
		C  = st_data(., "`clustvar'", "`temp_touse'") 
		dC = 1
		indC = order(C,1)
		g = rows(panelsetup(C[indC],1))
	}
				
	if ("`x'" == "`eval'") {
			eval    = st_data(., "`temp_grid'"	, "`temp_touse'") //; grid
	}
	else {	
			eval     = st_data(., "`temp_grid'"	, 0) //; grid
	}

	eval = sort(eval,1)
	N    = length(X)
    dups = dupsid = J(N,1,.) 	
	
	sort_id = 1::N
	
    if ("`vce_select'"=="nn") {
	
		sort_data = X,Y,sort_id
		sort_data = sort(sort_data,1)
		X = sort_data[,1]
		Y = sort_data[,2]
	
	  for (j=1; j<=`n'; j++) {
		*dups[j]=sum(X==X[j])
		dups[j] = length(selectindex(X:==X[j]))
      }
	  j=1
	  while(j<=`n') {
        dupsid[j::(j+dups[j]-1)]= 1::dups[j]
        j = j + dups[j]
      }
	  
	 }
	
	
	***********************************************************************
	c_bw = `C_c'*min((`x_sd',`x_iq'/1.349))*`n'^(-1/5)
		
    if ("`bwselect'"=="imse-rot" | "`bwselect'"=="mse-rot" |  "`bwselect'"=="all") {
      k = `p'+3	  
      rk = J(N,(k+1),.) 
	  for (j=1; j<=(k+1); j++) {
	    rk[.,j] = X:^(j-1)
	  }	
	  iGp     = cholinv(quadcross(rk,rk))
	  gamma_p = iGp*quadcross(rk,Y)
	  s2_p    = mean((rk*gamma_p:-Y):^2)
      
      k  = `q'+3
      rk = J(N,(k+1),.)
	  for (j=1; j<=(k+1); j++) {
	    rk[.,j] = X:^(j-1)
	  }
      iGq     = cholinv(quadcross(rk,rk))
	  gamma_q = iGq*quadcross(rk,Y)
	  s2_q    = mean((rk*gamma_q:-Y):^2)
    }
	
	if ("`bwselect'" == "imse-dpi" | "`bwselect'" == "imse-rot" | "`bwselect'" == "all") {	
	*eval_imse  = st_data(., "`temp_grid_imse'"	, "`temp_touse'") //; grid
	eval_imse  = st_data(., "`temp_grid_imse'"	, 0) //; grid
	eval_imse  = select(eval_imse, rowmissing(eval_imse):==0)
	neval_imse = length(eval_imse)
	bws_C_dpi  = bws_C_rot = J(neval_imse,4,.)
	
	if  ("`bwselect'"=="imse-dpi" |  "`bwselect'"=="all") {	
		for (e=1; e<=neval_imse; e++) {	
			bw_imse = nprobust_lp_mse_dpi(Y, X, C, eval_imse[e], `p', `q', `deriv', `even', "`kernel'", c_bw, `bwcheck', `bwregul', "`vce_select'", `nnmatch', dups, dupsid, "`interior'")         			 
			bws_C_dpi[e,.] = bw_imse[3]  , bw_imse[4] , bw_imse[5], bw_imse[6] 
		}
		    bws_Cmeans = mean(bws_C_dpi)
			b_imse_dpi = (bws_Cmeans[3]/(`n'*bws_Cmeans[4]))^(1/(2*`q'+3))
			h_imse_dpi = (bws_Cmeans[1]/(`n'*bws_Cmeans[2]))^(1/(2*`p'+3))			
	}
	 
	if  ("`bwselect'"=="imse-rot" |  "`bwselect'"=="all") {

		for (e=1; e<=neval_imse; e++) {		
		    *mp1 = gamma_p[`p'+1]*factorial(`p'+1) + gamma_p[`p'+2]*factorial(`p'+2)*eval_imse[e] + gamma_p[`p'+3]*factorial(`p'+3)*eval_imse[e]^2/2      
			*mq1 = gamma_q[`q'+1]*factorial(`q'+1) + gamma_q[`q'+2]*factorial(`q'+2)*eval_imse[e] + gamma_q[`q'+3]*factorial(`q'+3)*eval_imse[e]^2/2  

        	mp1 = nprobust_lp_coef(Y, X, eval_imse[e], `p'+1, `p'+3, `range', "`kernel'")
			mp2 = nprobust_lp_coef(Y, X, eval_imse[e], `p'+2, `p'+3, `range', "`kernel'")  
  
			n_h1 = sum(abs(X:-eval_imse[e]):<=c_bw)
			f0_pilot=n_h1/(N*c_bw)
        
			bws_C_rot[e,.]= s2_p/f0_pilot, mp1^2, s2_q/f0_pilot, mp2^2
      }	  
	    bws_Cmeans = mean(bws_C_rot)
		b_imse_rot = `C_c'*(bws_Cmeans[3]/(`n'*bws_Cmeans[4]))^(1/(2*`q'+3))
		h_imse_rot = `C_c'*(bws_Cmeans[1]/(`n'*bws_Cmeans[2]))^(1/(2*`p'+3))  

	}	
	}
	
	
	Result = J(`neval', 4, .)
	if("`bwselect'"=="all")  Result = J(`neval', 14, .)
	tcols = cols(Result)
	
	Result[., 1] = eval    
	
	*** Start loop over evaluation points	
	for (e=1; e<=`neval'; e++) {
	
		Result[e, tcols] = e		
		
		if ("`bwcheck'" != "0") {
			bw_min = sort(abs(X:-eval[e]), 1)[`bwcheck'] + 1e-8
			c_bw = max((c_bw, bw_min))
		}
		
    *** mse dpi
	bw_mse = nprobust_lp_mse_dpi(Y, X, C, eval[e], `p', `q', `deriv', `even', "`kernel'", c_bw, `bwcheck', `bwregul', "`vce_select'", `nnmatch', dups, dupsid, "`interior'") 
	h_mse_dpi = bw_mse[1]
	b_mse_dpi = bw_mse[2]
	
    *if (!is.null(rho))  b_mse_dpi = h_mse_dpi/rho
    Result[e,2::3] = h_mse_dpi,  b_mse_dpi	   
	  
      if  ("`bwselect'"=="mse-rot" | "`bwselect'"=="all") {
        *mp1 = gamma_p[`p'+1]*factorial(`p'+1) + gamma_p[`p'+2]*factorial(`p'+2)*eval[e] + gamma_p[`p'+3]*factorial(`p'+3)*eval[e]^2/2      
        *mq1 = gamma_q[`q'+1]*factorial(`q'+1) + gamma_q[`q'+2]*factorial(`q'+2)*eval[e] + gamma_q[`q'+3]*factorial(`q'+3)*eval[e]^2/2  
		mp1 = nprobust_lp_coef(Y, X, eval[e], `p'+1, `p'+3, `range', "`kernel'")
		mp2 = nprobust_lp_coef(Y, X, eval[e], `p'+2, `p'+3, `range', "`kernel'")
			
		n_h1 = sum(abs(X:-eval[e]):<=c_bw)
        f0_pilot=n_h1/(2*N*c_bw)
        
        b_mse_rot = `C_c'*((s2_q/f0_pilot)/(N*mp2^2))^(1/(2*`q'+3))
        h_mse_rot = `C_c'*((s2_p/f0_pilot)/(N*mp1^2))^(1/(2*`p'+3))		
        
		if ("`bwcheck'" != "0") {
			b_mse_rot = max((b_mse_rot, bw_min))
			h_mse_rot = max((h_mse_rot, bw_min))
		}		
		Result[e,2::3] = h_mse_rot, b_mse_rot        
      }
      
	  
	if  ("`bwselect'"=="ce-rot" |"`bwselect'"=="ce-dpi" | "`bwselect'"=="all") {
	  h_ce_rot = h_mse_dpi*`n'^(-(`p'/((3+`p')*(3+2*`p'))))
	  b_ce_rot = b_mse_dpi*`n'^(-(`q'/((3+`q')*(3+2*`q'))))
	  Result[e,2::3] = h_ce_rot, b_ce_rot
	}
	
	
      
  if  ("`bwselect'"=="ce-dpi" | "`bwselect'"=="all") { 
  
    h = h_mse_dpi
	b = h_mse_dpi
	*h=b=30
    rho = h/b
    Xh = (X:-eval[e]):/h
	Xb = (X:-eval[e]):/b
	
	Kh = nprobust_w_fun(Xh,"`kernel'")
	Lb = nprobust_w_fun(Xb,"`kernel'")
    indh = selectindex(Kh:> 0)
    indb = selectindex(Lb:> 0)
    ind = indh
    if (h>b) ind = indh   
    eN = length(ind)
    eY   =   Y[ind];eX  =   X[ind]
    eXh = Xh[ind];  eXb = Xb[ind]
    eKh = Kh[ind];  eLb = Lb[ind]
	
    Wp = eKh/h
	Wq = eLb/b
    Rp2 = J(eN,(`p'+3),.)
	  
	for (j=1; j<=(`p'+3); j++) {
	  Rp2[.,j] = eXh:^(j-1)
	  }
    Rp1  = Rp2[,1::(`p'+2)]
    Rp   = Rp2[,1::(`p'+1)]
    Rq   = J(eN,(`q'+1),.)
		
    for (j=1; j<=(`q'+1); j++) {
		Rq[,j] = eXb:^(j-1)
	}    	
	
    Lp1 = quadcross((Rp'*diag(Wp))', eXh:^(`p'+1))/eN 
    Lp2 = quadcross((Rp'*diag(Wp))', eXh:^(`p'+2))/eN 
    Lp3 = quadcross((Rp'*diag(Wp))', eXh:^(`p'+3))/eN
    Lq1 = quadcross((Rq'*diag(Wq))', eXb:^(`q'+1))/eN 
    Lq2 = quadcross((Rq'*diag(Wq))', eXb:^(`q'+2))/eN
    Lq3 = quadcross((Rq'*diag(Wq))', eXb:^(`q'+3))/eN 
    	
	invGp  = eN*cholinv(quadcross(Rp,Wp,Rp))
	invGq  = eN*cholinv(quadcross(Rq,Wq,Rq))
	
	hii = predicts_p = predicts_q = 0
    
    if ("`vce_select'"=="hc0" | "`vce_select'"=="hc1" | "`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
        Hq = J((`q'+1),1,.)
		for (j=1; j<=(`q'+1); j++) {
			Hq[j] = b^(-(j-1))
		}
        beta_q = diag(Hq)*invGq*quadcross((Rq'*diag(Wq))',eY)/eN
        rq = J(eN,`q'+1,.)
        predicts = J(eN,1,.)
		for (j=1; j<=(`q'+1); j++) {
		  rq[,j] = (eX:-eval[e]):^(j-1)
		  }
		for (j=1; j<=eN; j++) {
			predicts[j] = rq[j,]*beta_q
		}
			
        if ("`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
          hii=J(eN,1,.)	
          for (j=1; j<=eN; j++) { 
			hii[j] = (Rp[j,]*invGp*(Rp:*Wp)[j,]')/eN
			}
		}
		}


    res_q = nprobust_lp_res(eX, eY, predicts, hii, "`vce_select'", `nnmatch', dups[ind], dupsid[ind], `q'+1)
    	
	**** Bias 
	k = `p'+3
    rk = J(N,(k+3),.) 
	for (j=1; j<=(k+3); j++) {
	   rk[.,j] = X:^(j-1)
	}	
	iGp  = cholinv(quadcross(rk,rk))
	gamma = iGp*quadcross(rk,Y)		      	  
    mp3 = gamma[`p'+4]*factorial(`p'+3) + gamma[`p'+5]*factorial(`p'+4)*eval[e] + gamma[`p'+6]*factorial(`p'+5)*eval[e]^2/2  
	*gamma = gamma[1::(k+2)]
	rk=rk[,1::(k+2)]
	iGp  = cholinv(quadcross(rk,rk))
	gamma = iGp*quadcross(rk,Y)
	mp2 = gamma[`p'+3]*factorial(`p'+2) + gamma[`p'+4]*factorial(`p'+3)*eval[e] + gamma[`p'+5]*factorial(`p'+4)*eval[e]^2/2 
	    
    *if (is.na(m.p.3)) m.p.3 = lprobust(Y, X, h=range, eval=eval[e], p=p+4, deriv=(p+3), "`kernel'", "`vce_select'")$Estimate[5]
    *if (is.na(m.p.2)) m.p.2 = lprobust(Y, X, h=range, eval=eval[e], p=p+4, deriv=(p+2), "`kernel'", "`vce_select'")$Estimate[5]
    ep1 = J((`q'+1),1,0)
	ep1[`p'+2]=1
    e0   = J((`p'+1),1,0)
	e0[1]=1	
	
    q_rbc = nprobust_lp_q(eY, eX, eKh, eLb, res_q, eval[e], h, b, `p', `q', `deriv')

	  
    if ("`interior'"!="") {
      eta_bc1 = (e0'*invGp)*(   (mp2/factorial(`p'+2))*Lp2/h     + (mp3/factorial(`p'+3))*Lp3 )
      eta_bc2 = rho^(-2)*b^(`q'-`p'-1)*(e0'*invGp)*Lp1*ep1'*invGq*( (mp2/factorial(`p'+2))*Lq1/b     + (mp3/factorial(`p'+3))*Lq2 )
      eta_bc = (eta_bc1-eta_bc2)
	    
		S  = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_cer())
		optimize_init_params(S, h_mse_dpi)
		optimize_init_argument(S, 1, q_rbc)
		optimize_init_argument(S, 2, eta_bc)
		optimize_init_argument(S, 3, `p')
				
		h_ce_dpi = optimize(S)

      *H_bc = function(H) {abs(H^(-1)*q1_rbc + H^(1+2*(p+3))*eta_bc^2*q2_rbc + H^(p+3)*eta_bc*q3_rbc)}
      *h_bc <- optimize(H_bc , interval=c(0, `range'))
      *h_ce_dpi = h_Bc$minimum*N^(-1/(p+4))
    }
    
	
    if ("`interior'"=="") {
		eta_bc_1 = (e0'*invGp)*(Lp2 - rho^(-1)*Lp1*ep1'*invGq*Lq1)*(mp2/factorial(`p'+2))
		eta_bc_2 = (e0'*invGp)*(Lp3 - rho^(-2)*Lp1*ep1'*invGq*Lq2)*(mp3/factorial(`p'+3))
		phiz = 0.05844094
		E = q_rbc[1]*phiz, eta_bc_1, eta_bc_2, q_rbc[3]*phiz*eta_bc_1, q_rbc[3]*phiz*eta_bc_2
		q2_rbc = q_rbc[2]	 
		 
		S = optimize_init()
		optimize_init_evaluator(S, &nprobust_lp_cer())
		optimize_init_params(S, h_mse_dpi)
		optimize_init_argument(S, 1, E)
		optimize_init_argument(S, 2, q2_rbc)
		optimize_init_argument(S, 3, `p')
		optimize_init_argument(S, 4, `n')
		optimize_init_which(S, "min")	    
		optimize_init_evaluatortype(S, "d0")
		optimize_init_tracelevel(S, "none")
		optimize_init_conv_maxiter(S, 20)
		(void) optimize(S)
		h_ce_dpi = optimize_result_params(S)
    }
	
    b_ce_dpi =  b_mse_dpi	
	Result[e,2::3] = h_ce_dpi, b_ce_dpi	
  }
 
	if ("`bwselect'"=="all") Result[e,2::9] = h_mse_dpi, b_mse_dpi, h_mse_rot, b_mse_rot, h_ce_dpi, b_ce_dpi, h_ce_rot, b_ce_rot	
	
	
	 
  }

  if ("`bwselect'"=="imse-dpi" | "`bwselect'"=="all") {
    if ("`bwselect'"=="all") {
      Result[,10] = J(`neval',1,h_imse_dpi)
      Result[,11] = J(`neval',1,b_imse_dpi)
    } else {
      Result[,2]  = J(`neval',1,h_imse_dpi)
      Result[,3]  = J(`neval',1,b_imse_dpi)
    }
  }
    
  if ("`bwselect'"=="imse-rot" | "`bwselect'"=="all") {
    if ("`bwselect'"=="all") {
      Result[,12] = J(`neval',1,h_imse_rot)
      Result[,13] = J(`neval',1,b_imse_rot)
    } else{
      Result[,2] = J(`neval',1,h_imse_rot)
      Result[,3] = J(`neval',1,b_imse_rot)
    }
  }

  
  	genvars = st_local("genvars")
	if (genvars != "") {
		(void) 	st_addvar("double", 	    invtokens((genvars, "eval"), 	"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "eval"), 	"_"), Result[., 1])
		
		if ( "`bwselect'"!="all") {
		(void) 	st_addvar("double", 	    invtokens((genvars, "h"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h"), 		"_"), Result[., 2])
					
		(void) 	st_addvar("double",      	invtokens((genvars, "b"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b"), 		"_"), Result[., 3])
			}
			
		else {
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_mse_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_mse_dpi"), 		"_"), Result[., 2])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_mse_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_mse_dpi"), 		"_"), Result[., 3])
		
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_mse_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_mse_rot"), 		"_"), Result[., 4])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_mse_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_mse_rot"), 		"_"), Result[., 5])
		
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_ce_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_ce_dpi"), 		"_"), Result[., 6])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_ce_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_ce_dpi"), 		"_"), Result[., 7])
		
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_ce_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_ce_rot"), 		"_"), Result[., 8])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_ce_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_ce_rot"), 		"_"), Result[., 9])
		
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_imse_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_imse_dpi"), 		"_"), Result[., 10])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_imse_dpi"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_imse_dpi"), 		"_"), Result[., 11])
		
		(void) 	st_addvar("double", 	    invtokens((genvars, "h_imse_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "h_imse_rot"), 		"_"), Result[., 12])
					
		(void) 	st_addvar("double", 	    invtokens((genvars, "b_imse_rot"), 		"_"))
				st_store(Result[., tcols], 	invtokens((genvars, "b_imse_rot"), 		"_"), Result[., 13])
		}		
	}
	st_matrix("Result", Result[., 1..(tcols-1)])
  }
  
  ********************************************************************************
**** generate variable labels
*if ("`genvars'" != "") {
*   label variable `genvars'_grid 	"lpbwdensity: grid point"
*	label variable `genvars'_bw 	"lpbwdensity: bandwidth"
*	label variable `genvars'_nh 	"lpbwdensity: effective sample size"
*}    
	  **** display
	disp ""
	disp "Bandwidth Selection for Local Polynomial Density Estimation." 
	disp ""

	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %12.0f `n'
	disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %12.0f `p'
	disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %12.0f `deriv'
	disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 12: `kernel_type'}"
	disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 12: `bwselect'}"
	disp ""

	if ( "`bwselect'"!="all") {
		disp in smcl in gr "{hline 35}"
		disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: eval}" _col(14)  "{ralign 10: h}" _col(24) "{ralign 8: b}" _col(32)
		disp in smcl in gr "{hline 35}"
	}
	else {
		disp in smcl in gr "{hline 132}"
		disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: eval}"   ///
		_col(20)  "{ralign 10: mse-dpi}"   ///
		_col(40)  "{ralign 10: mse-rot}"  ///
		_col(60)  "{ralign 10: ce-dpi}"   ///
		_col(80)  "{ralign 10: ce-rot}"    ///
		_col(100)  "{ralign 10: imse-dpi}"  ///
		_col(120) "{ralign 10: imse-rot}"  

		disp in smcl in gr "{ralign 4: }" _col(4)   ///
		_col(15)  "{ralign 8: h}"  _col(24)  "{ralign 8: b}" ///
		_col(32)  "{ralign 8: h}"  _col(44)  "{ralign 8: b}" ///
		_col(52)  "{ralign 8: h}"  _col(64)  "{ralign 8: b}" ///
		_col(72)  "{ralign 8: h}"  _col(84)  "{ralign 8: b}" ///
		_col(92)  "{ralign 8: h}" _col(104) "{ralign 8: b}" ///
		_col(112) "{ralign 8: h}" _col(124) "{ralign 8: b}" 
		disp in smcl in gr "{hline 132}"

	}
*/

	forvalues i = 1(1)`neval' {
	if ( "`bwselect'"!="all") {
	disp in smcl in gr %4.0f `i' _col(4) in ye %10.3f Result[`i',1] _col(14)  %10.3f Result[`i', 2] _col(24) %10.3f Result[`i', 3] _col(32) 
		if (`i' != `neval' & `separator' > 1 & mod(`i', `separator') == 0) {
			disp in smcl in gr "{hline 35}"
		}
	} 
	else {
	disp in smcl in gr %4.0f `i'  _col(4) in ye %10.4f Result[`i',1] ///
	_col(14)  %9.3f Result[`i', 2]  _col(24)  %9.3f Result[`i', 3] ///
	_col(34)  %9.3f Result[`i', 4]  _col(44)  %9.3f Result[`i', 5] ///
	_col(54)  %9.3f Result[`i', 6]  _col(64)  %9.3f Result[`i', 7] ///
	_col(74)  %9.3f Result[`i', 8]  _col(84)  %9.3f Result[`i', 9] ///
	_col(94)  %9.3f Result[`i', 10] _col(104) %9.3f Result[`i', 11] ///
	_col(114) %9.3f Result[`i', 12] _col(124) %9.3f Result[`i', 13] 
	}

	}
	if ( "`bwselect'"!="all") {
	disp in smcl in gr "{hline 35}"
	}
	else {
	disp in smcl in gr "{hline 132}"
	}
  	
	ereturn clear
	ereturn scalar N = `n'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	
	ereturn local kernel = "`kernel_type'"
	ereturn local bwselect = "`bwselect'"
	ereturn local vce_select = "`vce_type'"
	ereturn local yvar "`y'"
	ereturn local xvar "`x'"
	ereturn local cmd "lpbwselect"
	
	*ereturn matrix bws = bws
	
	matrix colnames Result = eval h b
	ereturn matrix result = Result
	
	mata mata clear 

end


