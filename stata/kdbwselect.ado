*!version 1.0.0  2026-05-17

capture program drop kdbwselect
program define kdbwselect, eclass
	version 14.0
	syntax varlist(max=1) [if] [in] [, eval(varname) neval(real 0) kernel(string) bwselect(string) bwcheck(real 21) imsegrid(real 30) separator(integer 5)]

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
		local bwselect = "mse-dpi"
	} 
	else {
		if ("`bwselect'" != "mse-dpi" & "`bwselect'" != "imse-dpi" & "`bwselect'" != "mse-rot" & "`bwselect'" != "imse-rot" & "`bwselect'" != "all") {
			di as err `"bwselect(): incorrectly specified: options(mse-dpi, imse-dpi, mse-rot, imse-rot)"'
			exit 198
		}
	}

	********************************************************************************
	**** error check: kernel()
	if ("`kernel'" == "") {
		local kernel = "epanechnikov"
	}
	else {
		if ("`kernel'" != "uni" & "`kernel'" != "uniform" & "`kernel'" != "epa" & "`kernel'" != "epanechnikov") {
			di as err `"kernel(): incorrectly specified. Supported kernels for kdbwselect: epa, uni."'
			exit 198
		}
	}


	if ("`kernel'"=="epanechnikov" | "`kernel'"=="epa" | "`kernel'"=="") {
		local kernel_type = "Epanechnikov"
		local C_c = 2.34
		local C_h = 2.34
		local C_b = 3.49
	}
	else {
		local kernel_type = "Uniform"
		local C_c = 1.06
		local C_h = 1.06
		local C_b = 1
	}
	
	********************************************************************************
	**** error check: separator()
	if (`separator' <= 1) {
		local separator = 1
	}

  
  		
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
	bws = J(`neval', 2, .)
	Vh=Bh= J(`neval', 1, .)
	if("`bwselect'"=="all")  bws = J(`neval', 8, .)

	***********************************************************************

	* h_imse_rot/b_imse_rot are used in every bwselect path below (mse-dpi,
	* imse-dpi, mse-rot, imse-rot, all), so compute unconditionally.
	h_imse_rot = `x_sd'*`C_h'*`n'^(-1/(1+2*`p'))
	b_imse_rot = `x_sd'*`C_b'*`n'^(-1/(1+2*(`p'+2)+2*`p'))
	
		
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
			
	for (e=1; e<=`neval'; e++) {

		// bwcheck: per-eval floor on h/b set to the bwcheck-th smallest
		// |x-eval[i]| (mirrors R kdbwselect.R lines 146-150).
		bw_min = .
		if (`bwcheck' > 0) {
			bw_min = sort(abs(x:-eval[e]), 1)[min((`bwcheck', `n'))]
		}

		uh = (x:-eval[e])/h_imse_rot
		ub = (x:-eval[e])/b_imse_rot

		Kx = nprobust_K_fun(uh, "`kernel'")
		Lx = nprobust_L_fun(ub, "`kernel'")

		f_h = mean(Kx)/h_imse_rot
		f_b = mean(Lx)/b_imse_rot^(1+`p')

		Bh[e] = f_b*mK
		Vh[e] = f_h*vK

		h_mse_dpi = ((1+2*`deriv')*Vh[e]/(2*`p'*`n'*Bh[e]^2))^(1/(1+2*`p'+2*`deriv'))
		b_mse_dpi = b_imse_rot

		if (`bwcheck' > 0) {
			h_mse_dpi = max((h_mse_dpi, bw_min))
			b_mse_dpi = max((b_mse_dpi, bw_min))
		}

		bws[e,1::2] = h_mse_dpi, b_mse_dpi



		if ("`bwselect'"=="ce-rot"  | "`bwselect'"=="all") {
		  h_ce_rot = h_mse_dpi*`n'^(-(`p'-2)/((1+2*`p')*(1+`p'+2)))
		  b_ce_rot = b_mse_dpi*`n'^(-(`p'-2)/((1+2*`p')*(1+`p'+2)))

		  if (`bwcheck' > 0) {
			h_ce_rot = max((h_ce_rot, bw_min))
			b_ce_rot = max((b_ce_rot, bw_min))
		  }

		  bws[e,1::2] = h_ce_rot, b_ce_rot
		}

	  if("`bwselect'"=="all") bws[e,1::4] = h_mse_dpi,b_mse_dpi, h_ce_rot,b_ce_rot
    }




	if ("`bwselect'"=="imse-dpi" | "`bwselect'"=="all") {
	// IMSE averaging must be on a representative grid independent of the
	// user's eval points; mirror R kdbwselect imsegrid quantile grid
	// from p=0.1 to p=0.9.
	imsegrid = `imsegrid'
	qseq_imse = range(0.1, 0.9, 0.8/(imsegrid-1))
	xs = sort(x, 1)
	eval_imse = J(imsegrid, 1, .)
	for (gi=1; gi<=imsegrid; gi++) {
		// R quantile() type 7: linear interpolation on (k-1)/(n-1).
		qpos = qseq_imse[gi]*(`n'-1) + 1
		lo   = floor(qpos)
		hi   = ceil(qpos)
		if (lo==hi) eval_imse[gi] = xs[lo]
		else        eval_imse[gi] = xs[lo] + (qpos-lo)*(xs[hi]-xs[lo])
	}
	Vh_imse = J(imsegrid, 1, .)
	Bh_imse = J(imsegrid, 1, .)
	for (gi=1; gi<=imsegrid; gi++) {
		uh = (x:-eval_imse[gi])/h_imse_rot
		ub = (x:-eval_imse[gi])/b_imse_rot
		Kx = nprobust_K_fun(uh, "`kernel'")
		Lx = nprobust_L_fun(ub, "`kernel'")
		f_h = mean(Kx)/h_imse_rot
		f_b = mean(Lx)/b_imse_rot^(1+`p')
		Bh_imse[gi] = f_b*mK
		Vh_imse[gi] = f_h*vK
	}
    h_imse_dpi = ((1+2*`deriv')*mean(Vh_imse)/(2*`p'*`n'*mean(Bh_imse)^2))^(1/(1+2*`p'+2*`deriv'))
	b_imse_dpi = b_imse_rot
    if ("`bwselect'"=="all") {
      bws[,5] = J(`neval',1,h_imse_dpi)
      bws[,6] = J(`neval',1,b_imse_dpi)
    } else {
      bws[,1]  = J(`neval',1,h_imse_dpi)
      bws[,2]  = J(`neval',1,b_imse_dpi)
    }
  }

  
  if ("`bwselect'"=="imse-rot" | "`bwselect'"=="all") {
    if ("`bwselect'"=="all") {
      bws[,7] = J(`neval',1,h_imse_rot)
      bws[,8] = J(`neval',1,b_imse_rot)
    } else{
      bws[,1] = J(`neval',1,h_imse_rot)
      bws[,2] = J(`neval',1,b_imse_rot)
    }
 }  
  
  


   
  st_matrix("bws", bws)
  st_matrix("eval", eval)
  }
  
 
  
	**** display
	disp ""
	disp "Bandwidth Selection for Local Polynomial Density Estimation." 
	disp ""

	disp in smcl in gr "{lalign 1: Sample size                              (n=)    }" _col(19) in ye %12.0f `n'
	*disp in smcl in gr "{lalign 1: Polynomial order for point estimation    (p=)    }" _col(19) in ye %12.0f `p'
	*disp in smcl in gr "{lalign 1: Order of derivative estimated            (v=)    }" _col(19) in ye %12.0f `deriv'
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
		_col(40)  "{ralign 10: ce-rot}"  ///
		_col(60)  "{ralign 10: imse-dpi}"   ///
		_col(80)  "{ralign 10: imse-rot}"   
		
		disp in smcl in gr "{ralign 4: }" _col(4)   ///
		_col(15)  "{ralign 8: h}"  _col(24)  "{ralign 8: b}" ///
		_col(32)  "{ralign 8: h}"  _col(44)  "{ralign 8: b}" ///
		_col(52)  "{ralign 8: h}"  _col(64)  "{ralign 8: b}" ///
		_col(72)  "{ralign 8: h}"  _col(84)  "{ralign 8: b}" 
		disp in smcl in gr "{hline 132}"
	}


	forvalues i = 1(1)`neval' {
		if ( "`bwselect'"!="all") {
		disp in smcl in gr %4.0f `i' _col(4) in ye %10.4f eval[`i',1] _col(14)  %10.4f bws[`i', 1] _col(24) %10.4f bws[`i', 2] _col(32) 
			if (`i' != `neval' & `separator' > 1 & mod(`i', `separator') == 0) {
				disp in smcl in gr "{hline 35}"
			}
		} 
		else {
			disp in smcl in gr %4.0f `i' _col(4) in ye %10.4f eval[`i',1] ///
			_col(14)  %10.4f bws[`i', 1] _col(24) %10.4f bws[`i', 2] ///
			_col(34)  %10.4f bws[`i', 3] _col(44) %10.4f bws[`i', 4] ///
			_col(54)  %10.4f bws[`i', 5] _col(64) %10.4f bws[`i', 6] ///
			_col(74)  %10.4f bws[`i', 7] _col(84) %10.4f bws[`i', 8] 
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
	ereturn local kernel = "`kernel_type'"
	ereturn local bwselect = "`bwselect'"
	ereturn local xvar "`x'"
	ereturn local cmd "kdbwselect"
	ereturn matrix bws = bws
	ereturn matrix eval = eval
	
	mata mata clear 
	
end


