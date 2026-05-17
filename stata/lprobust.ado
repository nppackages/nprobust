*!version 1.0.0  2026-05-17

capture program drop lprobust
program define lprobust, eclass
	version 14.0
	syntax anything [if] [in] [, eval(varname) neval(real 0) deriv(real 0) p(real 1) h(string) b(string) rho(real 1) kernel(string) bwselect(string) bwcheck(real 21) imsegrid(real 30) vce(string) level(real 95) separator(integer 5) bwregul(real 1) interior genvars plot covgrid graph_options(string) weights(varname numeric) masspoints(string)]
	
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
		if ("`vce_select'"=="nn") local nnmatch `"`2'"'
		if inlist("`vce_select'","cluster","nncluster","cr1","cr2","cr3","hc0","hc1","hc2","hc3") local clustvar `"`2'"'
	}
	if `w' == 3 {
		local vce_select `"`1'"'
		local clustvar   `"`2'"'
		local nnmatch    `"`3'"'
		if !inlist("`vce_select'","cluster","nncluster","cr1","cr2","cr3","nn") {
			di as error "{err}{cmd:vce()} incorrectly specified"
			exit 125
		}
	}
	if `w' > 3 {
		di as error "{err}{cmd:vce()} incorrectly specified"
		exit 125
	}

	* Disallow vce(nncluster ...): warn and shift to cr1 (default when clusters).
	if ("`vce_select'"=="nncluster") {
		di as text "Warning: vce(nncluster) is not supported. Switching to vce(cr1) (the default when clusters)."
		local vce_select = "cr1"
	}

	* With a cluster variable, map hc0/hc1/hc2/hc3 to cr1/cr1/cr2/cr3.
	* Per cluster_validation design: hc0+cluster is a silent remap to cr1
	* (the default); hc1/2/3 produce a warning so the user knows their
	* requested HC variant is being upgraded to its cluster analogue.
	if ("`clustvar'"!="") {
		if ("`vce_select'"=="hc0") local vce_select = "cr1"
		if ("`vce_select'"=="hc1") {
			di as text "Warning: vce(hc1 `clustvar') is not a cluster option. Switching to vce(cr1 `clustvar')."
			local vce_select = "cr1"
		}
		if ("`vce_select'"=="hc2") {
			di as text "Warning: vce(hc2 `clustvar') is not a cluster option. Switching to vce(cr2 `clustvar')."
			local vce_select = "cr2"
		}
		if ("`vce_select'"=="hc3") {
			di as text "Warning: vce(hc3 `clustvar') is not a cluster option. Switching to vce(cr3 `clustvar')."
			local vce_select = "cr3"
		}
		* bare vce(cluster clustvar) is equivalent to cr1 (default).
		if ("`vce_select'"=="cluster") local vce_select = "cr1"
	}

	* cr1/cr2/cr3 require a cluster variable; otherwise error.
	if inlist("`vce_select'","cr1","cr2","cr3") {
		if ("`clustvar'"=="") {
			di as error "{err}{cmd:vce(`vce_select' clustervar)} requires a cluster variable"
			exit 125
		}
	}

	* Preserve the post-remap vce option for e(vce_select) before the
	* internal cr* -> hc* mapping below. Mirrors rdrobust convention.
	local vce_raw = "`vce_select'"
	if ("`vce_raw'"=="") local vce_raw = "nn"

	* Display label
	local vce_type = "NN"
	if ("`vce_select'"=="hc0") local vce_type = "HC0"
	if ("`vce_select'"=="hc1") local vce_type = "HC1"
	if ("`vce_select'"=="hc2") local vce_type = "HC2"
	if ("`vce_select'"=="hc3") local vce_type = "HC3"
	if ("`vce_select'"=="cr1") local vce_type = "CR1"
	if ("`vce_select'"=="cr2") local vce_type = "CR2"
	if ("`vce_select'"=="cr3") local vce_type = "CR3"

	* Any non-empty clustvar triggers cluster-robust variance.
	if ("`clustvar'"!="") local cluster = "cluster"

	* Internal mapping: cr1/cr2/cr3 -> hc1/hc2/hc3 so the existing residual
	* paths (hc-only) work unchanged. cr_type downstream reads the internal
	* name and applies the right multiplier in nprobust_lp_cluster_meat().
	if ("`vce_select'"=="cr1") local vce_select = "hc1"
	if ("`vce_select'"=="cr2") local vce_select = "hc2"
	if ("`vce_select'"=="cr3") local vce_select = "hc3"
	if ("`vce_select'"=="")    local vce_select = "nn"
		
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
		if ("`neval'" == "0") local neval = 30
		tempvar temp_grid
		qui gen double `temp_grid' = .
		if (`neval' == 1) {
			qui replace `temp_grid' = `x_min' in 1
		}
		else {
			forvalues i = 1/`neval' {
				local eval_i = `x_min' + (`x_max' - `x_min') * (`i' - 1) / (`neval' - 1)
				qui replace `temp_grid' = `eval_i' in `i'
			}
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

	**** error check: level()
	if (`level' >= 100 | `level' <= 0) {
		 di as error  "{err}{cmd:level()} should be a number in (0, 100)"
		 exit 125
	}

	**** error check: rho()
	if (`rho' < 0) {
		 di as error  "{err}{cmd:rho()} should be non-negative"
		 exit 125
	}

	**** error check: bwcheck()
	if (`bwcheck' < 0 | (`bwcheck' != round(`bwcheck'))) {
		 di as error  "{err}{cmd:bwcheck()} must be a non-negative integer"
		 exit 125
	}

	**** error check: masspoints()
	if ("`masspoints'" != "" & ///
	    !inlist("`masspoints'", "check", "off")) {
		 di as error  "{err}{cmd:masspoints()} must be 'check' or 'off'"
		 exit 125
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
	**** error check: masspoints()
	if ("`masspoints'" == "") local masspoints = "check"
	if ("`masspoints'" != "check" & "`masspoints'" != "off") {
		di as err `"masspoints(): incorrectly specified: options(check, off)"'
		exit 198
	}

	********************************************************************************
	**** weights()
	if ("`weights'" != "") {
		cap confirm numeric variable `weights'
		if _rc {
			di as err `"weights(): `weights' is not a numeric variable"'
			exit 198
		}
	}

	********************************************************************************
	**** error check: bw()
	if ("`h'" == "") {
		tempvar temp_h
		tempvar temp_b		
		*su `temp_grid'
		qui lpbwselect `y' `x' if `touse', eval(`temp_grid') deriv(`deriv') p(`p')  kernel(`kernel') vce(`vce')  bwregul(`bwregul') bwcheck(`bwcheck') imsegrid(`imsegrid') `interior' bwselect(`bwselect') genvars `=cond("`weights'"=="","","weights(`weights')")' masspoints(off)
		qui gen double `temp_h' = lpbwselect_h
		qui gen double `temp_b' = lpbwselect_b
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
				qui gen double `temp_h' = real("`h'") if `temp_grid' != .
				cap confirm numeric variable `temp_h'
				if _rc | missing(`temp_h'[1]) {
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

		* Handle b() option in manual-h mode. The auto-bwselect branch
		* (above) populates temp_b from lpbwselect output; the manual-h
		* branch must do the equivalent or downstream st_data calls hit
		* "varlist required" (rc=3598). Logic mirrors the h() handling:
		*  - b() as scalar    -> create tempvar, fill on temp_grid rows
		*  - b() as varname   -> use the user's variable directly
		*  - b() unspecified  -> default to b = h / rho (matches auto path)
		if ("`b'" != "") {
			cap confirm numeric variable `b'
			if _rc {
				capture confirm variable `b'
				if (!_rc) {
					di as err `"b(): `b' is not a numeric variable"'
					exit 198
				}
				else {
					tempvar temp_b
					qui gen double `temp_b' = real("`b'") if `temp_grid' != .
					cap confirm numeric variable `temp_b'
					if _rc | missing(`temp_b'[1]) {
						di as err `"b(): `b' is not numeric"'
						exit 198
					}
				}
			}
			else {
				if ("`x'" == "`grid'") {
					qui count if `b' != . & `touse'
				}
				else {
					qui count if `b' != .
				}
				local nb = r(N)
				if (`nb' != `neval') {
					di as err `"b(): `b' has different length as grid()"'
					exit 198
				}
				local temp_b "`b'"
			}
		}
		else {
			tempvar temp_b
			if (`rho' > 0) {
				qui gen double `temp_b' = `temp_h' / `rho' if `temp_grid' != .
			}
			else {
				qui gen double `temp_b' = `temp_h' if `temp_grid' != .
			}
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
		tempvar plot_tau_us
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

		Y = st_data(., "`y'", "`temp_touse'") //; y
		X = st_data(., "`x'", "`temp_touse'") //; x
		if ("`weights'"=="") {
			wts = J(rows(X), 1, 1)
		}
		else {
			wts = st_data(., "`weights'", "`temp_touse'")
		}

		if ("`x'" == "`eval'") {
			eval    = st_data(., "`temp_grid'"	, "`temp_touse'") //; grid
			h       = st_data(., "`temp_h'"	, "`temp_touse'") //; bw
			b       = st_data(., "`temp_b'"	, "`temp_touse'") //; bw
		}
		else {
			eval     = st_data(., "`temp_grid'"	, 0) //; grid
			h        = st_data(., "`temp_h'"	, 0) //; bw
			b        = st_data(., "`temp_b'"	, 0) //; bw
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
			// Sort (X, Y, wts[, C]) by X, keeping all vectors aligned.
			sort_data = X, Y, wts, sort_id
			if ("`cluster'"~="") {
				sort_data = X, Y, wts, C, sort_id
			}

			sort_data = sort(sort_data,1)
			X   = sort_data[,1]
			Y   = sort_data[,2]
			wts = sort_data[,3]
			if ("`cluster'"~="") {
				C = sort_data[,4]
				sort_id = sort_data[,5]
			}
			else {
				sort_id = sort_data[,4]
			}

			// O(N) duplicate counter on the already-sorted X via run detection.
			if (`n'==1) {
				dups   = 1
				dupsid = 1
			}
			else {
				bdy     = selectindex(X[2::`n'] :!= X[1::(`n'-1)])
				starts  = 1 \ (bdy:+1)
				ends    = (bdy \ `n')
				lens    = ends :- starts :+ 1
				dups    = J(`n',1,.)
				dupsid  = J(`n',1,.)
				for (j=1; j<=length(lens); j++) {
					dups[starts[j]::ends[j]]   = J(lens[j],1,lens[j])
					dupsid[starts[j]::ends[j]] = 1::lens[j]
				}
			}
			}
			
	
		*sort_data = X,Y,sort_id
		*sort_data = sort(sort_data,3)
		*X = sort_data[,1]
		*Y = sort_data[,2]
		
		Result = J(`neval', 11, .)
		

		*** Start loop over evaluation points
		for (c=1; c<=`neval'; c++) {
	
			w_h = nprobust_lp_kweight(X,eval[c],h[c],"`kernel'"):*wts;
			w_b = nprobust_lp_kweight(X,eval[c],b[c],"`kernel'"):*wts;
			ind_h = selectindex(w_h:> 0); ind_b = selectindex(w_b:> 0);
			N_h = length(ind_h);	N_b = length(ind_b)
			// Bug fix: union of h- and b-windows so obs only reached by the
			// larger bandwidth are included for the bias-corrected step.
			ind = selectindex((w_h:> 0) :| (w_b:> 0))

			if ("`masspoints'"=="check") {
				n_unique = rows(uniqrows(X[ind_h]))
				if (n_unique < `p'+5) {
					st_local("mp_warn_msg", sprintf(                                 ///
						"Only %f unique x values within bandwidth at eval=%f (p+5=%f); local polynomial may be unreliable. Set masspoints(off) to silence.", ///
						n_unique, eval[c], `p'+5))
					displayas("text")
					printf("{txt}Warning: %s{txt}\n", st_local("mp_warn_msg"))
				}
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
				if (("`vce_select'"=="hc2" | "`vce_select'"=="hc3") & "`cluster'"=="") {
					hii = rowsum((R_p*invG_p) :* (R_p:*W_h))
				}
			}

			if ("`cluster'"=="") {
				// Standard HC / NN sandwich (meat uses RX = R_p*W_h).
				res_h = nprobust_lp_res(eX, eY, predicts_p, hii, "`vce_select'", `nnmatch', edups, edupsid, `p'+1)
				if ("`vce_select'"=="nn") {
					res_b = res_h
				}
				else {
					res_b = nprobust_lp_res(eX, eY, predicts_q, hii, "`vce_select'", `nnmatch', edups, edupsid, `q'+1)
				}
				V_Y_cl = invG_p*nprobust_lp_vce(R_p:*W_h, res_h, eC, indC)*invG_p
				V_Y_bc = invG_p*nprobust_lp_vce(Q_q,      res_b, eC, indC)*invG_p
			}
			else {
				// Cluster-robust: vce selects the CR variant.
				//   vce(hc0 clvar) -> CR0     vce(hc1 clvar) -> CR1
				//   vce(hc2 clvar) -> CR2     vce(hc3 clvar) -> CR3
				//   vce(cluster/nncluster clvar) -> CR1 (default)
				cr_type = "CR1"
				if ("`vce_select'"=="hc0") cr_type = "CR0"
				if ("`vce_select'"=="hc1") cr_type = "CR1"
				if ("`vce_select'"=="hc2") cr_type = "CR2"
				if ("`vce_select'"=="hc3") cr_type = "CR3"

				if ("`vce_select'"=="nn") {
					res_h_raw = nprobust_lp_res(eX, eY, predicts_p, hii, "nn", `nnmatch', edups, edupsid, `p'+1)
					res_b_raw = res_h_raw
				}
				else {
					res_h_raw = eY - predicts_p
					res_b_raw = eY - predicts_q
				}

				sqrtW_h = sqrt(W_h)
				X_std_h = R_p :* sqrtW_h
				r_std_h = res_h_raw :* sqrtW_h

				// k_override = q+1 aligns CR1 df correction with the q-regression
				// that produced res_b_raw. Without it, k=cols(Q_q)=p+1 was used.
				meat_cl = nprobust_lp_cluster_meat(X_std_h, r_std_h, eC, indC, invG_p, cr_type, 0)
				meat_bc = nprobust_lp_cluster_meat(Q_q,     res_b_raw, eC, indC, invG_p, cr_type, `q'+1)

				V_Y_cl = invG_p * meat_cl * invG_p
				V_Y_bc = invG_p * meat_bc * invG_p
			}
			
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
    
		w_h_i = nprobust_lp_kweight(X,eval[i],h[i],"`kernel'"):*wts;
		w_b_i = nprobust_lp_kweight(X,eval[i],b[i],"`kernel'"):*wts;
		ind_h_i = selectindex(w_h_i:> 0);  ind_b_i = selectindex(w_b_i:> 0)
		N_h_i   = length(ind_h_i);  N_b_i = length(ind_b_i)
		ind_i   = selectindex((w_h_i:> 0) :| (w_b_i:> 0))

		w_h_j = nprobust_lp_kweight(X,eval[j],h[j],"`kernel'"):*wts;
		w_b_j = nprobust_lp_kweight(X,eval[j],b[j],"`kernel'"):*wts;
		// Bug fix: previous versions had an inline assignment inside selectindex.
		ind_h_j = selectindex(w_h_j:>0)
		ind_b_j = selectindex(w_b_j:>0)
		N_h_j   = length(ind_h_j);  N_b_j = length(ind_b_j)
		ind_j   = selectindex((w_h_j:> 0) :| (w_b_j:> 0))

		ind = uniqrows(ind_i \ ind_j)
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
			// Bug fix: previously `predicts_q_j` was assigned twice;
			// `predicts_q_i` was never set and later read as zero.
			predicts_p_i = R_p_i*beta_p_i;   predicts_q_i = R_q_i*beta_q_i
			predicts_p_j = R_p_j*beta_p_j;   predicts_q_j = R_q_j*beta_q_j
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
			
        // NOTE: V_us_i carries factorial(deriv) (NOT factorial(deriv)^2)
        // so that cov(tau_i, tau_j) = (V_us_i * V_us_j')[d+1,d+1] picks up
        // factorial(deriv)^2 from the cross-product -- matching the scaling
        // of se in the main loop. The previous code had factorial(deriv)^2
        // here, giving factorial(deriv)^4 after the cross-product (a factor
        // factorial(deriv)^2 too large for deriv >= 2).
        V_us_i =  factorial(`deriv')*invG_p_i*((res_h_i):*(R_p_i:*W_h_i))'
        V_us_j =  factorial(`deriv')*invG_p_j*((res_h_j):*(R_p_j:*W_h_j))'

        V_rb_i =  factorial(`deriv')*invG_p_i*(res_b_i:*Q_q_i)'
        V_rb_j =  factorial(`deriv')*invG_p_j*(res_b_j:*Q_q_j)'
        
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
					
		(void) 	st_addvar("double", 	invtokens((genvars, "N"), 			"_"))
				st_store(Result[., 11], invtokens((genvars, "N"), 			"_"), Result[., 4])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "tau_us"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "tau_us"), 		"_"), Result[., 5])
				
		(void) 	st_addvar("double", 	invtokens((genvars, "tau_bc"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "tau_bc"), 		"_"), Result[., 6])

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
				
		(void) st_addvar("double", 		"`plot_tau_us'")
				st_store(Result[., 11], "`plot_tau_us'", Result[., 5])
				
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
		label variable `genvars'_N 	"lprobust: effective sample size"
		label variable `genvars'_tau_us 	"lprobust: point estimate with pol. order `p'"
		label variable `genvars'_se_us 	"lprobust: standard error for f_p"
		*label variable `genvars'_CI_l_us 	"lprobust: left level-`level' CI, conventional"
		*label variable `genvars'_CI_r_us 	"lprobust: right level-`level' CI, conventional"
		label variable `genvars'_tau_bc   	"lprobust: point estimate with pol. order `q'"
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
		local plot_tau_us	"`genvars'_tau_us"
		local plot_cil_rb	"`genvars'_CI_l_rb"
		local plot_cir_rb 	"`genvars'_CI_r_rb"
	}


	if ("`plot'" != "") {
		if (`"`graph_options'"'=="" ) local graph_options = `"title("lprobust (p=`p', q=`q', deriv=`deriv')", color(gs0)) xtitle("`x'") ytitle("")"'
		twoway 	(rarea `plot_cil_rb' `plot_cir_rb' `plot_eval', sort color(gs11)) ///
				(line `plot_tau_us' `plot_eval',  lcolor(black) sort lwidth(medthin) lpattern(solid)),  ///
				legend(cols(2) order(2 "point estimate" 1 "`level'% C.I." )) `graph_options'
	}
		

	ereturn clear
	
	ereturn scalar N = `n'
	ereturn scalar p = `p'
	ereturn scalar q = `q'
	ereturn scalar level = `level'
	
	ereturn local kernel      = "`kernel_type'"
	ereturn local bwselect   = "`bwselect'"
	ereturn local vce_select = "`vce_raw'"
	ereturn local vce_type   = "`vce_type'"
	if ("`clustvar'"!="") ereturn local clustvar "`clustvar'"
	ereturn local yvar "`y'"
	ereturn local xvar "`x'"
	ereturn local cmd "lprobust"
	
	matrix colnames Result = eval h b N tau_us tau_bc se_us se_rb CI_l_rb CI_r_rb
	ereturn matrix Result = Result
	if ("`covgrid'" != "") {
		ereturn matrix cov_us = cov_us
		ereturn matrix cov_rb = cov_rb
	}
	
	mata mata clear
 
end
	





