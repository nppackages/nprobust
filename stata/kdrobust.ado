*!version 1.0.0  2026-05-17

capture program drop kdrobust
program define kdrobust, eclass
	version 14.0
	syntax varlist(max=1) [if] [in] [, eval(varname) neval(real 0) h(string) b(string) rho(real 1) kernel(string) bwselect(string) bwcheck(real 21) imsegrid(real 30) level(real 95) separator(integer 5) genvars plot graph_options(string)]
	
	marksample touse
	
	local kernel   = lower("`kernel'")
	local bwselect = lower("`bwselect'")
	local p = 2
	local deriv = 0	
		
	********************************************************************************
	**** error check: main variable
	local x "`varlist'"
	cap confirm numeric variable `x'
	if _rc {
		di as err `"`eval' is not a numeric variable"'
		exit 198
	}
	else {
		qui su `x' if `x' !=. & `touse', d
		local x_min = r(min)
		local x_max = r(max)
		local x_sd  = r(sd)
		local n = r(N)
		if (`n' == 0) {
			di as err `"`x' has length 0"'
			exit 198
		}
	}

	********************************************************************************
	**** error check: grid()
	if ("`eval'" == "") {
		if (`neval' == 0) local neval = 30
		* Quantile grid evenly spaced in [10, 90] — mirrors R seq(0.1, 0.9, length.out=neval)
		* and R's quantile() type 7 (linear interpolation). Done in Mata to avoid
		* a per-eval-point dataset pass via `replace ... in i`.
		tempvar temp_grid
		qui gen double `temp_grid' = .
		mata: xq = sort(st_data(., "`x'", "`touse'"), 1); ne = strtoreal(st_local("neval")); ///
		      probs = J(ne, 1, 0.1); if (ne > 1) probs = 0.1 :+ (0::ne-1) :* (0.8/(ne-1)); ///
		      pos = 1 :+ (rows(xq)-1) :* probs; lo = floor(pos); hi = ceil(pos); ///
		      g = pos :- lo; qs = xq[lo] :+ g :* (xq[hi] :- xq[lo]); ///
		      st_store((1::ne), "`temp_grid'", qs)
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
	**** error check: bwselect()
	if ("`bwselect'" == "") {
		if (`neval'==1) local bwselect =  "mse-dpi"
		if (`neval'>1)  local bwselect = "imse-dpi"
	} 
	else {
		if ("`bwselect'" != "mse-dpi" & "`bwselect'" != "imse-dpi" & "`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot" & "`bwselect'" != "ce-dpi" & "`bwselect'" != "ce-rot") {
			di as err `"bwselect(): incorrectly specified: options(mse-dpi, imse-dpi, mse-rot, imse-rot, ce-dpi, ce-rot)"'
			exit 198
		}
	}

	********************************************************************************
	**** error check: kernel()
	if ("`kernel'" == "") {
		local kernel = "epanechnikov"
	}
	else {
		// kdrobust bias-correction formulas are only implemented for the
		// Epanechnikov and Uniform kernels. Reject triangular/gaussian with
		// a clear error.
		if ("`kernel'" != "uni" & "`kernel'" != "uniform" & "`kernel'" != "epa" & "`kernel'" != "epanechnikov") {
			di as err `"kernel(): incorrectly specified. Supported kernels for kdrobust: epa, uni."'
			exit 198
		}
	}

	if      ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") local kernel_type = "Epanechnikov"
	else                                                     local kernel_type = "Uniform"

	********************************************************************************
	**** error check: level()
	if (`level' <= 0 | `level' >= 100) {
		di as err `"level(): incorrectly specified: should be between 0 and 100"'
		exit 198
	}

	********************************************************************************
	**** error check: bw()
	if ("`h'" == "") {
		*tempvar temp_bw
		*tempvar temp_lpbwdensity
		*qui lpbwdensity `x' if `touse', grid(`temp_grid') p(`p') v(`v') bwselect(`bwselect') kernel(`kernel') cweights(`temp_cweights') pweights(`temp_pweights') genvars(`temp_lpbwdensity') separator(1)
		qui kdbwselect `x' if `touse', eval(`temp_grid') kernel(`kernel') bwselect(`bwselect') bwcheck(`bwcheck') imsegrid(`imsegrid')
		*qui gen `temp_bw' = `temp_lpbwdensity'_bw
		*qui capture drop `temp_lpbwdensity'_grid `temp_lpbwdensity'_bw `temp_lpbwdensity'_nh
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
				tempvar temp_bw			
				qui gen `temp_bw' = "`h'" if `temp_grid' != .
				capture destring `temp_bw', replace
				cap confirm numeric variable `temp_bw'
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
			local temp_bw "`h'"
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
		tempvar plot_tau_us
		tempvar plot_cil_rb
		tempvar plot_cir_rb
	}

	if ("`genvars'" != "") local genvars = "kdrobust"

	********************************************************************************
	**** MATA
	tempvar temp_touse
	qui gen `temp_touse' = `touse'
			
	mata{
	
		x = st_data(., "`x'", "`temp_touse'") //; x
				
		if ("`x'" == "`eval'") {
			eval    = st_data(., "`temp_grid'"	, "`temp_touse'") //; grid
			*h       = st_data(., "`temp_bw'"	, "`temp_touse'") //; bw
		}
		else {	
			eval     = st_data(., "`temp_grid'"	, 0) //; grid
			*h        = st_data(., "`temp_bw'"	, 0) //; bw
		}
			
		eval = sort(eval,1)

		// Fetch bandwidths: from kdbwselect's e(bws) when h was not given,
		// or from the user-provided temp_bw variable when h is manual.
		if ("`h'"=="") {
			bws = st_matrix("e(bws)")
			if (rows(bws) == 0 | cols(bws) < 2) {
				errprintf("kdrobust: internal error — kdbwselect did not return a valid e(bws) matrix (rows=%g, cols=%g)\n", rows(bws), cols(bws))
				exit(error(198))
			}
			h = bws[,1]
			b = bws[,2]
		}
		else {
			h = st_data(., "`temp_bw'", 0)
			h = select(h, rowmissing(h):==0)
			if (rows(h) < `neval') h = J(`neval',1,h[1])
			b = h
		}

		if (`rho'>0)  {
			b = h/`rho'
		}
		
		*h = J(neval, 1, `h')
		*b = J(neval, 1, `h')
		rho = h:/b
		Result = J(`neval', 11, .)				
		quant = invnormal(1 - (100 - `level') / 200)
		
	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") {
				mK = 0.1
				vK = 0.6
				mL = 0.05555556
				vL = 35
	}
	else {
				mK = 0.1666667
				vK = 0.5
				mL = 0.07142857
				vL = 22.5				
	}
			
		*** Start loop over evaluation points
		for (c=1; c<=`neval'; c++) {
			if (`bwcheck' > 0) {
				bw_min = sort(abs(x :- eval[c]), 1)[`bwcheck']
				h[c] = max((h[c], bw_min))
				b[c] = max((b[c], bw_min))
				rho[c] = h[c]/b[c]
			}
		 
		    u = (x:-eval[c])/h[c]
			*nprobust_K(u        , "`kernel'", Kx, mK, vK)
			
			Kx = nprobust_K_fun(u,        "`kernel'")
			Lx = nprobust_L_fun(rho[c]*u, "`kernel'")
			*nprobust_K(rho[c]*u , "`kernel'", Lx, mL, vL)
					
			// Bug fix: bias-correction uses the k_v moment of the L kernel
			// (mL in Stata's naming), not of K (mK).
			Mx = Kx - rho[c]^(1+`p')*Lx*mL
			f_us = mean(Kx)/h[c]
			f_bc = mean(Mx)/h[c]
			se_us = sqrt((mean((Kx:^2)) - mean(Kx):^2)/(`n'*h[c]^2))
			se_rb = sqrt((mean((Mx:^2)) - mean(Mx):^2)/(`n'*h[c]^2))

			// Bug fix: effective N should reflect the union of h- and b-windows.
			nh = sum(abs(x :- eval[c]) :<= max((h[c], b[c])))
			
			Result[c, 1]  = eval[c]
			Result[c, 2]  = h[c]
			Result[c, 3]  = b[c]
			Result[c, 4]  = nh
			Result[c, 5]  = f_us
			Result[c, 6]  = f_bc
			Result[c, 7]  = se_us
			Result[c, 8]  = se_rb
			Result[c, 9]  = f_bc - quant*se_rb
			Result[c, 10] = f_bc + quant*se_rb
			Result[c, 11] = c
		
		*** End loop over evaluation points
		}		
		
		
	genvars = st_local("genvars")
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
		label variable `genvars'_eval 	"kdrobust: grid point"
		label variable `genvars'_h  	"kdrobust: bandwidth"
		label variable `genvars'_N 	"kdrobust: effective sample size"
		label variable `genvars'_tau_us 	"kdrobust: point estimate with pol. order `p'"
		label variable `genvars'_se_us 	"kdrobust: standard error for f_p"
		*label variable `genvars'_CI_l_us 	"kdrobust: left level-`level' CI, conventional"
		*label variable `genvars'_CI_r_us 	"kdrobust: right level-`level' CI, conventional"
		label variable `genvars'_tau_bc   	"kdrobust: point estimate with pol. order `q'"
		label variable `genvars'_se_rb 	    "kdrobust: standard error for f_q"
		label variable `genvars'_CI_l_rb 	"kdrobust: left level-`level' CI, robust bias-corrected"
		label variable `genvars'_CI_r_rb 	"kdrobust: right level-`level' CI, robust bias-corrected"		
	}

	********************************************************************************
	**** display
	disp ""
	disp "Kernel Density Estimation and Inference." 
	disp ""

	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %15.0f `n'
	*disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %15.0f `p'
	*disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %15.0f `deriv'
	*disp in smcl in gr "{lalign 1: Polynomial order for confidence interval (q=)    }" _col(19) in ye %15.0f `q'
	disp in smcl in gr "{lalign 1: Kernel function                                  }" _col(19) in ye "{ralign 15: `kernel_type'}"
	disp in smcl in gr "{lalign 1: Bandwidth selection method                       }" _col(19) in ye "{ralign 15: `bwselect'}"
	disp ""

	disp in smcl in gr "{hline 72}"
	disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: }" _col(14) "{ralign 10: }" _col(24) "{ralign 8: }" _col(32) "{ralign 10: Point}" _col(42) "{ralign 10: Std.}" _col(52) "{ralign 20: Robust B.C.}" _col(72)
	disp in smcl in gr "{ralign 4: }" _col(4) "{ralign 10: grid}" _col(14)  "{ralign 10: bw}" _col(24) "{ralign 8: Eff.n}" _col(32) "{ralign 10: Est.}" _col(42) "{ralign 10: Error}" _col(52) "{ralign 20: `level'% Conf. Interval}" _col(72)
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
		local plot_tau_us	    "`genvars'_tau_us"
		local plot_cil_rb	"`genvars'_CI_l_rb"
		local plot_cir_rb 	"`genvars'_CI_r_rb"
	}

	if ("`plot'" != "") {
		if (`"`graph_options'"'=="" ) local graph_options = `"title("kdrobust", color(gs0)) xtitle("`x'") ytitle("")"'
		twoway 	(rarea `plot_cil_rb' `plot_cir_rb' `plot_eval', sort color(gs11)) ///
				(line `plot_tau_us' `plot_eval',  lcolor(black) sort lwidth(medthin) lpattern(solid)),  ///
				legend(cols(2) order(2 "point estimate" 1 "`level'% C.I." )) `graph_options'
	}
	

	********************************************************************************
	**** returns
	ereturn clear
	ereturn scalar N = `n'
	ereturn scalar level = `level'
	
	ereturn local kernel = "`kernel_type'"
	ereturn local bwselect = "`bwselect'"
	ereturn local xvar "`x'"
	ereturn local cmd "kdrobust"
	
	matrix colnames Result = eval h b N tau_us tau_bc se_us se_rb CI_l_rb CI_r_rb
	ereturn matrix Result = Result
	
	mata mata clear	
 
end
	





