*!version 0.3.1  06-20-2020

capture program drop lprobust 
program define lprobust, eclass
	syntax anything [if] [in] [, eval(varname) neval(real 0) deriv(real 0) p(real 1) h(string) b(string) rho(real 1) kernel(string) bwselect(string) bwcheck(real 21) vce(string) level(real 95) separator(integer 5) bwregul(real 1) interior genvars plot covgrid graph_options(string)]
	
	*disp in yellow "Preparing data." 
	marksample touse
	
	*preserve
	*qui keep if `touse'
	tokenize "`anything'"
	local y `1'
	local x `2'
	local kernel   = lower("`kernel'")
	local bwselect = lower("`bwselect'")
	local vce      = lower("`vce'")
	
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
		
	if ("`deriv'">"0" & "`p'"=="1") local p = `deriv'+1
	local q = `p'+1

	*******************************************************************************
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
		local n = r(N)
		local x_iq = r(p75)-r(p25)
		local range = `x_max' - `x_min'
		if (`n' == 0) {
			di as err `"`x' has length 0"'
			exit 198
		}
	}

	********************************************************************************
	**** error check: grid()
	if ("`eval'" == "") {
		if ("`neval'" == "0") {
			*tempvar temp_grid temp_dup
			*qui duplicates tag `x' if `touse', gen(`temp_dup')
			*qui gen `temp_grid' = `x' if `temp_dup'==0
			*qui su `temp_grid'
			*local neval = r(N)
			local neval = 30
			tempvar temp_grid	
			local nquant = `neval'+1
			pctile `temp_grid' = `x' if `touse', nq(`nquant')
			}
			else {
				tempvar temp_grid	
				local nquant = `neval'+1
				pctile `temp_grid' = `x' if `touse', nq(`nquant')
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
			//disp `ng'
			if (`neval' == 0) {
				di as err `"eval(): `eval' has length 0"'
				exit 198
			}
		}
		local temp_grid "`eval'"		
	}

	********************************************************************************
	**** error check: vce()
	if ("`vce_select'"~="nn" & "`vce_select'"~="" & "`vce_select'"~="cluster" & "`vce_select'"~="nncluster" & "`vce_select'"~="hc1" & "`vce_select'"~="hc2" & "`vce_select'"~="hc3" & "`vce_select'"~="hc0"){ 
		 di as error  "{err}{cmd:vce()} incorrectly specified"  
		 exit 7
	}
		
			
	********************************************************************************
	**** error check: bwselect()
	if ("`bwselect'" == "") {
		if (`neval'==1) local bwselect =  "mse-dpi"
		if (`neval'>1)  local bwselect = "imse-dpi"
	} 
	else {
		if ("`bwselect'" != "mse-dpi" & "`bwselect'" != "imse-dpi" & "`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot" & "`bwselect'" != "ce-dpi" & "`bwselect'" != "ce-rot" & "`bwselect'" != "all") {
			di as err `"bwselect(): incorrectly specified: options(mse-dpi, imse-dpi, mse-rot, imse-rot, ce-dpi, ce-rot)"'
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
	**** error check: level()
	if (`level' <= 0 | `level' >= 100) {
		di as err `"level(): incorrectly specified: should be between 0 and 100"'
		exit 198
	}

	********************************************************************************
	**** error check: bw()
	if ("`h'" == "") {
		tempvar temp_h
		tempvar temp_b		
		*su `temp_grid'
		qui lpbwselect `y' `x' if `touse', eval(`temp_grid') deriv(`deriv') p(`p')  kernel(`kernel') vce(`vce')  bwregul(`bwregul') bwcheck(`bwcheck') `interior' bwselect(`bwselect') genvars 		
		qui gen `temp_h' = lpbwselect_h
		qui gen `temp_b' = lpbwselect_b
		if (`rho'>0)  {
				qui replace `temp_b' = lpbwselect_h/`rho'
		}
		qui capture drop lpbwselect_*
	}
	else {
		cap confirm numeric variable `h'
		if _rc {
			// check if variable exists
			capture confirm variable `h'
			if (!_rc) {
				di as err `"h(): `h' is not a numeric variable"'
				exit 198
			} 
			else {
				tempvar temp_h			
				qui gen `temp_h' = "`h'" if `temp_grid' != .
				capture destring `temp_h', replace
				cap confirm numeric variable `temp_h'
				if _rc {
					di as err `"h(): `h' is not numeric"'
					exit 198
				}
			}
		}
		else {
			if ("`x'" == "`grid'") {
				qui count if `h' != . & `touse'
			}
			else {
				qui count if `h' != .
			}
			local nb = r(N)
			if (`nb' != `neval') {
				di as err `"h(): `h' has different length as grid()"'
				exit 198
			}
			local temp_h "`h'"
			local bwselect = "Manual"
		}
	}


	********************************************************************************
	**** error check: separator()
	if (`separator' <= 1) {
		local separator = 1
	}

	********************************************************************************
	**** temporaty varaibles for plotting
	if ("`plot'" != "" & "`genvars'" == "") {
		tempvar plot_eval
		tempvar plot_gx_us
		tempvar plot_cil_rb
		tempvar plot_cir_rb
	}

	if ("`genvars'" != "") local genvars = "lprobust"

	************
	*** MATA ***
	************
	
	tempvar temp_touse
	qui gen `temp_touse' = `touse'

	mata{
	
		Y = st_data(., "`y'", "`temp_touse'") //; x
		X = st_data(., "`x'", "`temp_touse'") //; x
				
		if ("`x'" == "`eval'") {
			eval    = st_data(., "`temp_grid'"	, "`temp_touse'") //; grid
			h       = st_data(., "`temp_h'"	, "`temp_touse'") //; bw
			b       = st_data(., "`temp_b'"	, "`temp_touse'") //; bw
		}
		else {	
			eval     = st_data(., "`temp_grid'"	, 0) //; grid
			h        = st_data(., "`temp_h'"	, 0) //; bw
			b        = st_data(., "`temp_h'"	, 0) //; bw
		}
		
		C = 0
		if ("`cluster'"~="") {
			C = st_data(., "`clustvar'", "`temp_touse'") 
			dC=1
		}
				
			eval = sort(eval,1)
					
			N = length(X)
			dups = dupsid = J(N,1,.) 
			sort_id = 1::N
			
			if ("`vce_select'"=="nn") {
			*dups      = st_data(.,("`dups'"), 0); dupsid    = st_data(.,("`dupsid'"), 0)
			sort_data = X,Y,sort_id
			if ("`cluster'"~="") {
				sort_data = X,Y,C,sort_id
			}
			
			sort_data = sort(sort_data,1)
			X = sort_data[,1]
			Y = sort_data[,2]
			sort_id = sort_data[,3]
			if ("`cluster'"~="") {
				C = sort_data[,3]
				sort_id = sort_data[,4]
			}
			
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
			
	
		*sort_data = X,Y,sort_id
		*sort_data = sort(sort_data,3)
		*X = sort_data[,1]
		*Y = sort_data[,2]
		
		Result = J(`neval', 11, .)
		

		*** Start loop over evaluation points
		for (c=1; c<=`neval'; c++) {
	
			w_h = nprobust_lp_kweight(X,eval[c],h[c],"`kernel'");	
			w_b = nprobust_lp_kweight(X,eval[c],b[c],"`kernel'");	
			ind_h = selectindex(w_h:> 0);ind_b = selectindex(w_b:> 0);
			N_h = length(ind_h);	N_b = length(ind_b)
			ind = ind_b
			if (h>b) {
				ind = ind_h   
			}
			eN = length(ind)
			eY  = Y[ind];
			eX  = X[ind];
			W_h = w_h[ind];
			W_b = w_b[ind];
		
			eC=indC=0
			if ("`cluster'"~="") {
				eC  = C[ind]
				indC = order(eC,1)
				g = rows(panelsetup(eC[indC],1))
			}
		
			*edups = edups = 0	
			edups   = dups[ind];edupsid = dupsid[ind]
			
			u = (eX:-eval[c])/h[c]
			R_q = J(eN,(`q'+1),.)
				for (j=1; j<=(`q'+1); j++)  {
					R_q[.,j] = (eX:-eval[c]):^(j-1)
				}
			R_p = R_q[,1::(`p'+1)]
					
			*************************************************************************
			*** display("Computing lp estimates.")
			*************************************************************************
			L = quadcross(R_p:*W_h,u:^(`p'+1))
			invG_q  = cholinv(quadcross(R_q,W_b,R_q))
			invG_p  = cholinv(quadcross(R_p,W_h,R_p))
			
			if (rank(invG_p)==. | rank(invG_q)==.){
				display("{err}Invertibility problem: check variability of running variable around cutoff")
				exit(1)
			}
			
			e_p1 = J((`q'+1),1,0); e_p1[`p'+2]=1
			e_v  = J((`p'+1),1,0); e_v[`deriv'+1]=1
	
			Q_q = ((R_p:*W_h)' - h[c]^(`p'+1)*(L*e_p1')*((invG_q*R_q')':*W_b)')'
			D = eY
			beta_p = invG_p*quadcross(R_p:*W_h,D); beta_q = invG_q*quadcross(R_q:*W_b,D); beta_bc = invG_p*quadcross(Q_q,D) 
				
			tau_cl = factorial(`deriv')*beta_p[(`deriv'+1),1]
			tau_bc = factorial(`deriv')*beta_bc[(`deriv'+1),1]
			s_Y = 1
			bias = tau_cl-tau_bc
		
			****************************************************************************
			*display("Computing variance-covariance matrix.")
			****************************************************************************
			hii=predicts_p=predicts_q=0
			if ("`vce_select'"=="hc0" | "`vce_select'"=="hc1" | "`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
				predicts_p=R_p*beta_p
				predicts_q=R_q*beta_q
				if ("`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
					hii=J(eN,1,.)	
						for (i=1; i<=eN; i++) {
							hii[i] = R_p[i,]*invG_p*(R_p:*W_h)[i,]'
					}
				}
			}
			
			res_h = nprobust_lp_res(eX, eY, predicts_p, hii, "`vce_select'", `nnmatch', edups, edupsid, `p'+1)
			if ("`vce_select'"=="nn") {
					res_b = res_h
			}
			else {
					res_b = nprobust_lp_res(eX, eY, predicts_q, hii, "`vce_select'", `nnmatch', edups, edupsid, `q'+1)				
			}

			V_Y_cl = invG_p*nprobust_lp_vce(R_p:*W_h,  res_h, eC, indC)*invG_p
			V_Y_bc = invG_p*nprobust_lp_vce(Q_q,       res_b, eC, indC)*invG_p
			
			V_tau_cl = factorial(`deriv')^2*(V_Y_cl)[`deriv'+1,`deriv'+1]
			V_tau_rb = factorial(`deriv')^2*(V_Y_bc)[`deriv'+1,`deriv'+1]
			se_tau_cl = sqrt(V_tau_cl);	se_tau_rb = sqrt(V_tau_rb)
			quant = -invnormal(abs((1-(`level'/100))/2))
			
			Result[c, 1]  = eval[c]
			Result[c, 2]  = h[c]
			Result[c, 3]  = b[c]
			Result[c, 4]  = eN
			Result[c, 5]  = tau_cl
			Result[c, 6]  = tau_bc
			Result[c, 7]  = se_tau_cl
			Result[c, 8]  = se_tau_rb
			Result[c, 9]  = tau_bc - quant*se_tau_rb
			Result[c, 10] = tau_bc + quant*se_tau_rb
			Result[c, 11] = c
			*** End loop over evaluation points
		}	
		
		
		 cov_us = cov_rb = J(`neval',`neval',.)
 
  if ("`covgrid'" != "") {
  
	for (i=1; i<=`neval'; i++) {
			for (j=i; j<=`neval'; j++) {  
    
		w_h_i = nprobust_lp_kweight(X,eval[i],h[i],"`kernel'");	
		w_b_i = nprobust_lp_kweight(X,eval[i],b[i],"`kernel'");		
		ind_h_i = selectindex(w_h_i:> 0);  ind_b_i = selectindex(w_b_i:> 0)
		N_h_i   = length(ind_h_i);  N_b_i = length(ind_b_i)
        ind_i   = ind_b_i
		if (h[i]>b[i]) ind_i = ind_h_i   
        
		w_h_j = nprobust_lp_kweight(X,eval[j],h[j],"`kernel'");	
		w_b_j = nprobust_lp_kweight(X,eval[j],b[j],"`kernel'");	
        ind_h_j = selectindex(w_h_j:>0); ind_b_j = selectindex(ind_b_j = w_b_j:>0)
        N_h_j   = length(ind_h_j);  N_b_j = length(ind_b_j)
        ind_j   = ind_b_j
        if (h[j]>b[j]) ind_j = ind_h_j   
        		
		ind = uniqrows(ind_i \ind_j)
        *ind = ind_i
				
        eN  = length(ind)
        eY  = Y[ind];     eX  = X[ind]
        
        W_h_i = w_h_i[ind];        W_b_i = w_b_i[ind]
        W_h_j = w_h_j[ind];        W_b_j = w_b_j[ind]
            			
		eC=indC=0
		if ("`cluster'"~="") {
			eC  = C[ind]
			indC = order(eC,1)
			g = rows(panelsetup(eC[indC],1))
		}
        
        edups = edupsid = 0	
        if ("`vce_select'"=="nn") {
          edups   = dups[ind]
          edupsid = dupsid[ind]
        }
        		
        u_i   = (eX:-eval[i])/h[i]
		u_j   = (eX:-eval[j])/h[j]
        R_q_i = R_q_j = J(eN,(`q'+1),.)
		for (k=1; k<=`q'+1; k++) {
			R_q_i[,k] = (eX:-eval[i]):^(k-1)
			R_q_j[,k] = (eX:-eval[j]):^(k-1)
        }
		R_p_i = R_q_i[,1::(`p'+1)]
		R_p_j = R_q_j[,1::(`p'+1)]
			
		L_i = quadcross(R_p_i:*W_h_i,u_i:^(`p'+1))        
        invG_p_i = cholinv(quadcross(R_p_i,W_h_i,R_p_i))
		invG_q_i = cholinv(quadcross(R_q_i,W_b_i,R_q_i))
		Q_q_i = ((R_p_i:*W_h_i)' - h[i]^(`p'+1)*(L_i*e_p1')*((invG_q_i*R_q_i')':*W_b_i)')'
		beta_p_i = invG_p_i*quadcross(R_p_i:*W_h_i,eY); beta_q_i = invG_q_i*quadcross(R_q_i:*W_b_i,eY); beta_bc_i = invG_p_i*quadcross(Q_q_i,eY) 

		L_j = quadcross(R_p_j:*W_h_j,u_j:^(`p'+1))
		invG_p_j = cholinv(quadcross(R_p_j,W_h_j,R_p_j))
		invG_q_j = cholinv(quadcross(R_q_j,W_b_j,R_q_j))
		Q_q_j = ((R_p_j:*W_h_j)' - h[j]^(`p'+1)*(L_j*e_p1')*((invG_q_j*R_q_j')':*W_b_j)')'		       
		beta_p_j = invG_p_j*quadcross(R_p_j:*W_h_j,eY); beta_q_j = invG_q_j*quadcross(R_q_j:*W_b_j,eY); beta_bc_j = invG_p_j*quadcross(Q_q_j,eY) 
       
		
		hii_i=predicts_p_i=predicts_q_i=0
		hii_j=predicts_p_j=predicts_q_j=0
		if ("`vce_select'"=="hc0" | "`vce_select'"=="hc1" | "`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
			predicts_p_i = R_p_i*beta_p_i;	predicts_q_j = R_q_i*beta_q_i
			predicts_p_j = R_p_j*beta_p_j;  predicts_q_j = R_q_j*beta_q_j
			if ("`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
				hii_i=hii_j=J(eN,1,.)	
					for (k=1; k<=eN; k++) {
						hii_i[k] = R_p_i[k,]*invG_p_i*(R_p_i:*W_h_i)[k,]'
						hii_j[k] = R_p_j[k,]*invG_p_j*(R_p_j:*W_h_j)[k,]'
				}
			}
		}		
			
		res_h_i = nprobust_lp_res(eX, eY, predicts_p_i, hii_i, "`vce_select'", `nnmatch', edups, edupsid, `p'+1)
		res_h_j = nprobust_lp_res(eX, eY, predicts_p_j, hii_j, "`vce_select'", `nnmatch', edups, edupsid, `p'+1)
        if ("`vce_select'"=="nn") {
               res_b_i = res_h_i
			   res_b_j = res_h_j
        } else {
			res_b_i = nprobust_lp_res(eX, eY, predicts_q_i, hii_i, "`vce_select'", `nnmatch', edups, edupsid, `q'+1)
			res_b_j = nprobust_lp_res(eX, eY, predicts_q_j, hii_j, "`vce_select'", `nnmatch', edups, edupsid, `q'+1)
         }
			
        V_us_i =  factorial(`deriv')^2*invG_p_i*((res_h_i):*(R_p_i:*W_h_i))'
        V_us_j =  factorial(`deriv')^2*invG_p_j*((res_h_j):*(R_p_j:*W_h_j))'
        
        V_rb_i =  factorial(`deriv')^2*invG_p_i*(res_b_i:*Q_q_i)'
        V_rb_j =  factorial(`deriv')^2*invG_p_j*(res_b_j:*Q_q_j)'
        
        cov_us[i,j] = (V_us_i*V_us_j')[`deriv'+1,`deriv'+1]
        cov_rb[i,j] = (V_rb_i*V_rb_j')[`deriv'+1,`deriv'+1]
        
        cov_us[j,i]= cov_us[i,j]
        cov_rb[j,i]= cov_rb[i,j]
		 
      }
    }
    
	 st_matrix("cov_us", cov_us)
	 st_matrix("cov_rb", cov_rb)
  }
  
 
				
	
	genvars = st_local("genvars")
	*genvars = "lprobust"
	if (genvars != "") {
		(void) 	st_addvar("double", 	invtokens((genvars, "eval"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "eval"), 		"_"), Result[., 1])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "h"), 			"_"))
				st_store(Result[., 11], invtokens((genvars, "h"), 			"_"), Result[., 2])
				
		(void) 	st_addvar("double", 	invtokens((genvars, "b"), 			"_"))
				st_store(Result[., 11], invtokens((genvars, "b"), 			"_"), Result[., 3])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "nh"), 			"_"))
				st_store(Result[., 11], invtokens((genvars, "nh"), 			"_"), Result[., 4])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "gx_us"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "gx_us"), 		"_"), Result[., 5])
				
		(void) 	st_addvar("double", 	invtokens((genvars, "gx_bc"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "gx_bc"), 		"_"), Result[., 6])

		(void) 	st_addvar("double", 	invtokens((genvars, "se_us"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "se_us"), 		"_"), Result[., 7])
				
		(void) 	st_addvar("double", 	invtokens((genvars, "se_rb"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "se_rb"), 		"_"), Result[., 8])
									
		(void) 	st_addvar("double", 	invtokens((genvars, "CI_l_rb"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "CI_l_rb"), 		"_"), Result[., 9])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "CI_r_rb"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "CI_r_rb"), 		"_"), Result[., 10])
	}
	else if ("`plot'" != "") {
		(void) 	st_addvar("double", 	"`plot_eval'")
				st_store(Result[., 11], "`plot_eval'", Result[., 1])
				
		(void) st_addvar("double", 		"`plot_gx_us'")
				st_store(Result[., 11], "`plot_gx_us'", Result[., 5])
				
		(void) st_addvar("double", 		"`plot_cil_rb'")
				st_store(Result[., 11], "`plot_cil_rb'", Result[., 9])
				
		(void) st_addvar("double", 		"`plot_cir_rb'")
				st_store(Result[., 11], "`plot_cir_rb'", Result[., 10])
	}
	
		st_matrix("Result", Result[., 1..10])
}
	
	
	********************************************************************************
	**** generate variable labels
	if ("`genvars'" != "") {
		label variable `genvars'_eval 	"lprobust: eval point"
		label variable `genvars'_h  	"lprobust: bandwidth"
		label variable `genvars'_nh 	"lprobust: effective sample size"
		label variable `genvars'_gx_us 	"lprobust: point estimate with pol. order `p'"
		label variable `genvars'_se_us 	"lprobust: standard error for f_p"
		*label variable `genvars'_CI_l_us 	"lprobust: left level-`level' CI, conventional"
		*label variable `genvars'_CI_r_us 	"lprobust: right level-`level' CI, conventional"
		label variable `genvars'_gx_bc   	"lprobust: point estimate with pol. order `q'"
		label variable `genvars'_se_rb 	    "lprobust: standard error for f_q"
		label variable `genvars'_CI_l_rb 	"lprobust: left level-`level' CI, robust bias-corrected"
		label variable `genvars'_CI_r_rb 	"lprobust: right level-`level' CI, robust bias-corrected"		
	}

		

	********************************************************************************
	**** display
	disp ""
	disp "Local Polynomial Regression Estimation and Inference." 
	disp ""

	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %15.0f `n'
	disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %15.0f `p'
	disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %15.0f `deriv'
	disp in smcl in gr "{lalign 1: Polynomial order for confidence interval (q=)    }" _col(19) in ye %15.0f `q'
	disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 15: `kernel_type'}"
	disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 15: `bwselect'}"
	*if ("`cluster'"!="")     disp in smcl in gr "{lalign 1: Std. Err. adjusted for clusters in "} "`clustvar'"
	disp ""

	disp in smcl in gr "{hline 72}"
	disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: }" _col(14) "{ralign 10: }" _col(24) "{ralign 8: }" _col(32) "{ralign 10: Point}" _col(42) "{ralign 10: Std.}" _col(52) "{ralign 20: Robust B.C.}" _col(72)
	disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: eval}" _col(14)  "{ralign 10: bw}" _col(24) "{ralign 8: Eff.n}" _col(32) "{ralign 10: Est.}" _col(42) "{ralign 10: Error}" _col(52) "{ralign 20: `level'% Conf. Interval}" _col(72)
	disp in smcl in gr "{hline 72}"

	forvalues i = 1(1)`neval' {
	disp in smcl in gr %4.0f `i' _col(4) in ye %10.4f Result[`i', 1] _col(14)  %10.4f Result[`i', 2] _col(24) %8.0f Result[`i', 4] _col(32) %10.4f Result[`i', 5] _col(42) %10.4f Result[`i', 7] _col(52) ///
		in ye %10.4f Result[`i', 9] _col(62) in ye %10.4f Result[`i', 10] _col(72)
		if (`i' != `neval' & `separator' > 1 & mod(`i', `separator') == 0) {
			disp in smcl in gr "{hline 72}"
		}
	}
	disp in smcl in gr "{hline 72}"


	********************************************************************************
	**** plot
	if ("`plot'" != "" & "`genvars'" != "") {
		local plot_eval	    "`genvars'_eval"
		local plot_gx_us	"`genvars'_gx_us"
		local plot_cil_rb	"`genvars'_CI_l_rb"
		local plot_cir_rb 	"`genvars'_CI_r_rb"
	}


	if ("`plot'" != "") {
		if (`"`graph_options'"'=="" ) local graph_options = `"title("lprobust (p=`p', q=`q', deriv=`deriv')", color(gs0)) xtitle("`x'") ytitle("")"'
		twoway 	(rarea `plot_cil_rb' `plot_cir_rb' `plot_eval', sort color(gs11)) ///
				(line `plot_gx_us' `plot_eval',  lcolor(black) sort lwidth(medthin) lpattern(solid)),  ///
				legend(cols(2) order(2 "point estimate" 1 "`level'% C.I." )) `graph_options'
	}
		

	ereturn clear
	
	ereturn scalar N = `n'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	ereturn scalar level = `level'
	
	ereturn local kernel      = "`kernel_type'"
	ereturn local bwselect   = "`bwselect'"
	ereturn local vce_select = "`vce_type'"
	ereturn local yvar "`y'"
	ereturn local xvar "`x'"
	ereturn local cmd "lprobust"
	
	matrix colnames Result = eval h b nh gx_us gx_bc se_us se_rb CI_l_rb CI_r_rb
	ereturn matrix Result = Result
	if ("`covgrid'" != "") {
		ereturn matrix cov_us = cov_us
		ereturn matrix cov_rb = cov_rb
	}
	
	mata mata clear
 
end
	





