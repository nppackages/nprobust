*!version 0.3.1  06-20-2020

capture program drop kdrobust 
program define kdrobust, eclass
	syntax varlist(max=1) [if] [in] [, eval(varname) neval(real 0) h(string) b(string) rho(real 1) kernel(string) bwselect(string) level(real 95) separator(integer 5) genvars plot graph_options(string)]
	
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
		if ("`neval'" == "0") {
			tempvar temp_grid temp_dup
			qui duplicates tag `x' if `touse', gen(`temp_dup')
			qui gen `temp_grid' = `x' if `temp_dup'==0
			qui su `temp_grid'
			local neval = r(N)
			}
			else {
				qui pctile `temp_grid' = `x' if `touse', nq(`neval'+1)
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
		if ("`kernel'" != "triangular" & "`kernel'" != "uniform" & "`kernel'" != "epanechnikov") {
			di as err `"kernel(): incorrectly specified: options(triangular, uniform, epanechnikov)"'
			exit 198
		}
	}

	if      ("`kernel'"=="epanechnikov" | "`kernel'"=="epa") local kernel_type = "Epanechnikov"
	else if ("`kernel'"=="uniform"      | "`kernel'"=="uni") local kernel_type = "Uniform"
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
		*tempvar temp_bw
		*tempvar temp_lpbwdensity
		*qui lpbwdensity `x' if `touse', grid(`temp_grid') p(`p') v(`v') bwselect(`bwselect') kernel(`kernel') cweights(`temp_cweights') pweights(`temp_pweights') genvars(`temp_lpbwdensity') separator(1)
		qui kdbwselect `x' if `touse', eval(`temp_grid') kernel(`kernel') bwselect(`bwselect') 
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
		tempvar plot_f_us
		tempvar plot_cil_rb
		tempvar plot_cir_rb
	}

	if ("`genvars'" != "") local genvars = "lprobust"

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
		bws = st_matrix("e(bws)")
		h = bws[,1]
		b = bws[,2]
		
		if (`rho'>0)  {
			b = h/`rho'
		}
		
		*h = J(neval, 1, `h')
		*b = J(neval, 1, `h')
		rho = h:/b
		Result = J(`neval', 11, .)				
		
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
		 
		    u = (x:-eval[c])/h[c]
			*nprobust_K(u        , "`kernel'", Kx, mK, vK)
			
			Kx = nprobust_K_fun(u,        "`kernel'")
			Lx = nprobust_L_fun(rho[c]*u, "`kernel'")
			*nprobust_K(rho[c]*u , "`kernel'", Lx, mL, vL)
					
			Mx = Kx - rho[c]^(1+`p')*Lx*mK
			f_us = mean(Kx)/h[c]
			f_bc = mean(Mx)/h[c]
			se_us = sqrt((mean((Kx:^2)) - mean(Kx):^2)/(`n'*h[c]^2))
			se_rb = sqrt((mean((Mx:^2)) - mean(Mx):^2)/(`n'*h[c]^2))

			nh = sum(abs(x :- eval[c]) :<= h[c])
			
			Result[c, 1]  = eval[c]
			Result[c, 2]  = h[c]
			Result[c, 3]  = b[c]
			Result[c, 4]  = nh
			Result[c, 5]  = f_us
			Result[c, 6]  = f_bc
			Result[c, 7]  = se_us
			Result[c, 8]  = se_rb
			Result[c, 9]  = f_bc - 2*se_rb
			Result[c, 10] = f_bc + 2*se_rb
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
					
		(void) 	st_addvar("double", 	invtokens((genvars, "nh"), 			"_"))
				st_store(Result[., 11], invtokens((genvars, "nh"), 			"_"), Result[., 4])
					
		(void) 	st_addvar("double", 	invtokens((genvars, "f_us"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "f_us"), 		"_"), Result[., 5])
				
		(void) 	st_addvar("double", 	invtokens((genvars, "f_bc"), 		"_"))
				st_store(Result[., 11], invtokens((genvars, "f_bc"), 		"_"), Result[., 6])

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
				
		(void) st_addvar("double", 		"`plot_f_us'")
				st_store(Result[., 11], "`plot_f_us'", Result[., 5])
				
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
		label variable `genvars'_nh 	"kdrobust: effective sample size"
		label variable `genvars'_f_us 	"kdrobust: point estimate with pol. order `p'" 
		label variable `genvars'_se_us 	"kdrobust: standard error for f_p"
		*label variable `genvars'_CI_l_us 	"kdrobust: left level-`level' CI, conventional"
		*label variable `genvars'_CI_r_us 	"kdrobust: right level-`level' CI, conventional"
		label variable `genvars'_f_bc   	"kdrobust: point estimate with pol. order `q'"
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
		local plot_f_us	    "`genvars'_f_us"
		local plot_cil_rb	"`genvars'_CI_l_rb"
		local plot_cir_rb 	"`genvars'_CI_r_rb"
	}

	if ("`plot'" != "") {
		if (`"`graph_options'"'=="" ) local graph_options = `"title("kdrobust", color(gs0)) xtitle("`x'") ytitle("")"'
		twoway 	(rarea `plot_cil_rb' `plot_cir_rb' `plot_eval', sort color(gs11)) ///
				(line `plot_f_us' `plot_eval',  lcolor(black) sort lwidth(medthin) lpattern(solid)),  ///
				legend(cols(2) order(2 "point estimate" 1 "`level'% C.I." )) `graph_options'
	}
	

	********************************************************************************
	**** returns
	ereturn clear
	ereturn scalar N = `n'
	ereturn scalar level = `level'
	
	ereturn local kernel = "`kernel_type'"
	ereturn local bwselect = "`bwselect'"
	ereturn local vce_select = "`vce_type'"
	ereturn local varname "`x'"
	ereturn local cmd "kdrobust"
	
	matrix colnames Result = eval h b nh f_us f_bc se_us se_rb CI_l_rb CI_r_rb
	ereturn matrix Result = Result
	
	mata mata clear	
 
end
	





