*!version 1.0.0  2026-05-17

capture program drop lpbwselect
program define lpbwselect, eclass
	version 14.0
	syntax anything [if] [in] [, eval(varname) deriv(real 0) neval(real 0) p(real 1) kernel(string) bwselect(string) separator(integer 5) bwregul(real 1) vce(string) genvars interior bwcheck(real 21) imsegrid(real 30) weights(varname numeric) masspoints(string)]

	if ("`masspoints'" == "") local masspoints = "check"
	if ("`masspoints'" != "check" & "`masspoints'" != "off") {
		di as err `"masspoints(): incorrectly specified: options(check, off)"'
		exit 198
	}
	if ("`weights'" != "") {
		cap confirm numeric variable `weights'
		if _rc {
			di as err `"weights(): `weights' is not a numeric variable"'
			exit 198
		}
	}

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
	* paths (hc-only) work unchanged in the bw pilot.
	if ("`vce_select'"=="cr1") local vce_select = "hc1"
	if ("`vce_select'"=="cr2") local vce_select = "hc2"
	if ("`vce_select'"=="cr3") local vce_select = "hc3"
	if ("`vce_select'"=="")    local vce_select = "nn"


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
		local range = `x_max' - `x_min'
		* Stata's summarize, detail uses a different quantile rule than R's
		* quantile() type 7. Compute IQR in Mata using R's rule so c_bw
		* (the pilot bandwidth) matches R lpbwselect exactly.
		mata: _xs = sort(st_data(., "`x'", "`touse'"), 1); _N = rows(_xs); ///
		      _q25 = 1 + 0.25*(_N-1); _q75 = 1 + 0.75*(_N-1); ///
		      _p25 = _xs[floor(_q25)] + (_q25-floor(_q25))*(_xs[ceil(_q25)]-_xs[floor(_q25)]); ///
		      _p75 = _xs[floor(_q75)] + (_q75-floor(_q75))*(_xs[ceil(_q75)]-_xs[floor(_q75)]); ///
		      st_local("x_iq", strofreal(_p75-_p25, "%21.16e"))
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
		* IMSE-DPI in R uses seq(x_min, x_max, length=imsegrid) (range-based);
		* IMSE-ROT in R uses quantile(x, seq(0.05, 0.95, 0.025)) (quantile).
		* Build the DPI grid here; Mata builds the ROT grid separately below.
		local nquant = `imsegrid'
		qui gen double `temp_grid_imse' = `x_min' + (`x_max' - `x_min')*(_n - 1)/(`nquant' - 1) if _n <= `nquant'
		if (`nquant' == 1) qui replace `temp_grid_imse' = `x_min' in 1
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

	* kernel_type is set together with C_c below; the duplicate mapping
	* that used to live here was removed (H5).

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
	// For the even case we use the closed-form CE-DPI scaling rather than
	// the full DPI computation, matching the R implementation.
	local ce_dpi_even = 0
	if (`even'==1 & "`bwselect'"=="ce-dpi") {
		local ce_dpi_even = 1
	}


	************
	*** MATA ***
	************

	tempvar temp_touse
	qui gen `temp_touse' = `touse'

	mata{

	Y = st_data(., "`y'", "`temp_touse'")
	X = st_data(., "`x'", "`temp_touse'")
	if ("`weights'"=="") {
		wts = J(rows(X), 1, 1)
	}
	else {
		wts = st_data(., "`weights'", "`temp_touse'")
	}
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
		// Sort by (X, sort_id) to make this a STABLE sort. Mata sort()
		// is not stable on a single column, so ties in X were getting
		// reordered relative to R order(), giving different y[pos] in NN
		// residual computation and a ~1% V_V drift at boundary eval points.
		sort_data = sort(sort_data,(1,3))
		X = sort_data[,1]
		Y = sort_data[,2]

	  // O(N) duplicate counter on the already-sorted X via run detection.
	  if (`n'==1) {
	    dups   = 1
	    dupsid = 1
	  }
	  else {
	    bdy    = selectindex(X[2::`n'] :!= X[1::(`n'-1)])
	    starts = 1 \ (bdy:+1)
	    ends   = (bdy \ `n')
	    lens   = ends :- starts :+ 1
	    dups   = J(`n',1,.)
	    dupsid = J(`n',1,.)
	    for (j=1; j<=length(lens); j++) {
	      dups[starts[j]::ends[j]]   = J(lens[j],1,lens[j])
	      dupsid[starts[j]::ends[j]] = 1::lens[j]
	    }
	  }

	 }


	***********************************************************************
	c_bw_init = `C_c'*min((`x_sd',`x_iq'/1.349))*`n'^(-1/5)
	c_bw = c_bw_init

    if ("`bwselect'"=="imse-rot" | "`bwselect'"=="mse-rot" |  "`bwselect'"=="all") {
      k = `p'+3
      rk = J(N,(k+1),.)
	  for (j=1; j<=(k+1); j++) {
	    rk[.,j] = X:^(j-1)
	  }
	  iGp     = cholinv(quadcross(rk,rk))
	  gamma_p = iGp*quadcross(rk,Y)
	  // R uses sigma(lm)^2 = RSS/(N-k-1) (unbiased), not mean(resid^2).
	  s2_p    = sum((rk*gamma_p:-Y):^2) / (N - (k+1))

      k  = `q'+3
      rk = J(N,(k+1),.)
	  for (j=1; j<=(k+1); j++) {
	    rk[.,j] = X:^(j-1)
	  }
      iGq     = cholinv(quadcross(rk,rk))
	  gamma_q = iGq*quadcross(rk,Y)
	  s2_q    = sum((rk*gamma_q:-Y):^2) / (N - (k+1))
    }

	if ("`bwselect'" == "imse-dpi" | "`bwselect'" == "imse-rot" | "`bwselect'" == "all") {
	// Read the range-based IMSE grid that the .ado built above. R's
	// lpbwselect.imse.dpi (npfunctions.R:668) uses
	//   eval <- seq(x.min, x.max, length.out=imsegrid)
	// while R's lpbwselect.imse.rot (npfunctions.R:704) uses
	//   eval <- quantile(x, probs = seq(0.05, 0.95, 0.025))
	// We build the imse-rot grid separately in Mata below.
	eval_imse  = st_data(., "`temp_grid_imse'", 0)
	eval_imse  = select(eval_imse, rowmissing(eval_imse):==0)
	neval_imse = length(eval_imse)
	bws_C_dpi  = bws_C_rot = J(neval_imse,4,.)

	if  ("`bwselect'"=="imse-dpi" |  "`bwselect'"=="all") {
		for (e=1; e<=neval_imse; e++) {
			// Reset c_bw per iteration. nprobust_lp_mse_dpi internally
			// bumps c_bw via bwcheck and the bumped value persists across
			// iterations (Mata real scalar args may pass by reference for
			// matrix-stored scalars), giving wrong c_bw for interior eval
			// points after a boundary one is processed.
			c_bw = c_bw_init
			bw_imse = nprobust_lp_mse_dpi(Y, X, C, eval_imse[e], `p', `q', `deriv', `even', "`kernel'", c_bw, `bwcheck', `bwregul', "`vce_select'", `nnmatch', dups, dupsid, "`interior'", wts)
			bws_C_dpi[e,.] = bw_imse[3]  , bw_imse[4] , bw_imse[5], bw_imse[6]
		}
		    bws_Cmeans = mean(bws_C_dpi)
			b_imse_dpi = (bws_Cmeans[3]/(`n'*bws_Cmeans[4]))^(1/(2*`q'+3))
			h_imse_dpi = (bws_Cmeans[1]/(`n'*bws_Cmeans[2]))^(1/(2*`p'+3))
	}

	if  ("`bwselect'"=="imse-rot" |  "`bwselect'"=="all") {
		// R lpbwselect.imse.rot (npfunctions.R:701-725) uses a 37-point
		// quantile grid (probs = seq(0.05, 0.95, 0.025)) and per-point
		// V, B from lpbwselect.mse.rot, with formula
		//   h = ((1+2*deriv) mean(V) / (N * 2*(p+1-deriv) * mean(B^2)))^(1/(2p+3))
		// when even=FALSE; optimize otherwise.
		// b uses the same function with (p=q, deriv=p+1).
		_xs_rot = sort(X, 1)
		_nrot = 37  // (0.95-0.05)/0.025 + 1 = 37
		eval_rot = J(_nrot, 1, .)
		for (gi=1; gi<=_nrot; gi++) {
			_qp = 0.05 + (gi-1)*0.025
			_qpos = _qp*(N-1) + 1
			_lo = floor(_qpos); _hi = ceil(_qpos)
			if (_lo==_hi) eval_rot[gi] = _xs_rot[_lo]
			else          eval_rot[gi] = _xs_rot[_lo] + (_qpos-_lo)*(_xs_rot[_hi]-_xs_rot[_lo])
		}

		Khr = nprobust_lp_constants(`p', `deriv', "`kernel'")
		Kbr = nprobust_lp_constants(`q', `p'+1,   "`kernel'")
		facph1r = factorial(`p'+1); facph2r = factorial(`p'+2)
		facqh1r = factorial(`q'+1); facqh2r = factorial(`q'+2)
		even_h_rot = mod(`p' - `deriv', 2) == 0
		even_b_rot = mod(`q' - (`p'+1), 2) == 0

		Vh_arr  = J(_nrot, 1, .); Bh_arr  = J(_nrot, 1, .)
		Vb_arr  = J(_nrot, 1, .); Bb_arr  = J(_nrot, 1, .)
		for (e=1; e<=_nrot; e++) {
			mp1_h_r = nprobust_lp_coef(Y, X, eval_rot[e], `p'+1, `p'+3, `range', "`kernel'")
			mp2_h_r = nprobust_lp_coef(Y, X, eval_rot[e], `p'+2, `p'+3, `range', "`kernel'")
			mp1_b_r = nprobust_lp_coef(Y, X, eval_rot[e], `q'+1, `q'+3, `range', "`kernel'")
			mp2_b_r = nprobust_lp_coef(Y, X, eval_rot[e], `q'+2, `q'+3, `range', "`kernel'")
			n_h1r = sum(abs(X:-eval_rot[e]):<=c_bw_init)
			f0r   = n_h1r/(2*N*c_bw_init)

			Vh_arr[e] = Khr[2]*(s2_p/f0r)
			B1h_r = Khr[1]*mp1_h_r/facph1r
			B2h_r = Khr[1]*mp2_h_r/facph2r
			if (even_h_rot == 0) {
				Bh_arr[e] = B1h_r
			} else {
				// Per-point h.mse.rot via bracketed search; B = B1 + h*B2.
				h_pt = nprobust_lp_brent_mserot(Vh_arr[e], (B1h_r, B2h_r), `p', `deriv', N, 1e-12, `range')
				Bh_arr[e] = B1h_r + h_pt*B2h_r
			}

			Vb_arr[e] = Kbr[2]*(s2_q/f0r)
			B1b_r = Kbr[1]*mp1_b_r/facqh1r
			B2b_r = Kbr[1]*mp2_b_r/facqh2r
			if (even_b_rot == 0) {
				Bb_arr[e] = B1b_r
			} else {
				b_pt = nprobust_lp_brent_mserot(Vb_arr[e], (B1b_r, B2b_r), `q', `p'+1, N, 1e-12, `range')
				Bb_arr[e] = B1b_r + b_pt*B2b_r
			}
		}

		mVh = mean(Vh_arr); mBh2 = mean(Bh_arr:^2)
		mVb = mean(Vb_arr); mBb2 = mean(Bb_arr:^2)
		if (even_h_rot == 0) {
			h_imse_rot = ((1+2*`deriv')*mVh / (N*2*(`p'+1-`deriv')*mBh2))^(1/(2*`p'+3))
		} else {
			h_imse_rot = nprobust_lp_brent_mserot(mVh, (sqrt(mBh2), 0), `p', `deriv', N, 1e-12, `range')
		}
		if (even_b_rot == 0) {
			b_imse_rot = ((2*(`p'+1)+1)*mVb / (N*2*(`q'-`p')*mBb2))^(1/(2*`q'+3))
		} else {
			b_imse_rot = nprobust_lp_brent_mserot(mVb, (sqrt(mBb2), 0), `q', `p'+1, N, 1e-12, `range')
		}
	}
	}


	Result = J(`neval', 4, .)
	if("`bwselect'"=="all")  Result = J(`neval', 14, .)
	tcols = cols(Result)

	Result[., 1] = eval

	*** Start loop over evaluation points
	for (e=1; e<=`neval'; e++) {

		Result[e, tcols] = e

		// Reset c_bw per eval point. The previous code reused the
		// monotonically-bumped c_bw across iterations, which made
		// per-eval-point bandwidths drift away from R lpbwselect
		// (where c.bw is recomputed inside lpbwselect.mse.dpi each call).
		c_bw = c_bw_init

		if ("`bwcheck'" != "0") {
			// R uses bw.min = sort(|x-eval|)[bwcheck] without the +1e-8.
			// The fudge factor caused boundary X's at distance exactly bw.min
			// to be included (u<1), expanding eN by the tied count and
			// drifting V_V by ~1% at boundary eval points.
			bw_min = sort(abs(X:-eval[e]), 1)[`bwcheck']
			c_bw = max((c_bw, bw_min))
		}

    *** mse dpi
	bw_mse = nprobust_lp_mse_dpi(Y, X, C, eval[e], `p', `q', `deriv', `even', "`kernel'", c_bw, `bwcheck', `bwregul', "`vce_select'", `nnmatch', dups, dupsid, "`interior'", wts)
	h_mse_dpi = bw_mse[1]
	b_mse_dpi = bw_mse[2]

    *if (!is.null(rho))  b_mse_dpi = h_mse_dpi/rho
    Result[e,2::3] = h_mse_dpi,  b_mse_dpi

      if  ("`bwselect'"=="mse-rot" | "`bwselect'"=="all") {
		// R lpbwselect calls lpbwselect.mse.rot TWICE (lpbwselect.R:122-128):
		//   h: lpbwselect.mse.rot(p=p, deriv=deriv) — m.p.* use poly p+3
		//   b: lpbwselect.mse.rot(p=q, deriv=p+1)  — m.p.* use poly q+3
		// Within lpbwselect.mse.rot, even=(p-deriv)%%2==0 selects between
		// closed-form (FALSE) and optimize (TRUE).
		mp1_h = nprobust_lp_coef(Y, X, eval[e], `p'+1, `p'+3, `range', "`kernel'")
		mp2_h = nprobust_lp_coef(Y, X, eval[e], `p'+2, `p'+3, `range', "`kernel'")
		mp1_b = nprobust_lp_coef(Y, X, eval[e], `q'+1, `q'+3, `range', "`kernel'")
		mp2_b = nprobust_lp_coef(Y, X, eval[e], `q'+2, `q'+3, `range', "`kernel'")

		// R lpbwselect.mse.rot (npfunctions.R:847-849) uses the raw
		// Silverman c.bw without the bwcheck floor. Use c_bw_init here
		// instead of the loop-adjusted c_bw, otherwise boundary eval
		// points (where bwcheck bumps c_bw up) get a smaller f0_pilot
		// and hence a smaller h_mse_rot.
		c_bw_rot = c_bw_init
		n_h1 = sum(abs(X:-eval[e]):<=c_bw_rot)
        f0_pilot=n_h1/(2*N*c_bw_rot)

		Kh = nprobust_lp_constants(`p', `deriv', "`kernel'")
		Kb = nprobust_lp_constants(`q', `p'+1,   "`kernel'")
		facph1 = factorial(`p'+1); facph2 = factorial(`p'+2)
		facqh1 = factorial(`q'+1); facqh2 = factorial(`q'+2)
		Vh = Kh[2]*(s2_p/f0_pilot)
		B1h = Kh[1]*mp1_h/facph1
		B2h = Kh[1]*mp2_h/facph2
		Vb = Kb[2]*(s2_q/f0_pilot)
		B1b = Kb[1]*mp1_b/facqh1
		B2b = Kb[1]*mp2_b/facqh2

		even_h = mod(`p' - `deriv', 2) == 0
		even_b = mod(`q' - (`p'+1), 2) == 0

		if (even_h == 0) {
		    h_mse_rot = ((2*`deriv'+1)*Vh / (2*(`p'+1-`deriv')*B1h^2*N))^(1/(2*`p'+3))
		} else {
		    h_mse_rot = nprobust_lp_brent_mserot(Vh, (B1h, B2h), `p', `deriv', N, 1e-12, `range')
		}

		if (even_b == 0) {
		    b_mse_rot = ((2*(`p'+1)+1)*Vb / (2*(`q'-`p')*B1b^2*N))^(1/(2*`q'+3))
		} else {
		    b_mse_rot = nprobust_lp_brent_mserot(Vb, (B1b, B2b), `q', `p'+1, N, 1e-12, `range')
		}

		if ("`bwcheck'" != "0") {
			b_mse_rot = max((b_mse_rot, bw_min))
			h_mse_rot = max((h_mse_rot, bw_min))
		}
		Result[e,2::3] = h_mse_rot, b_mse_rot
      }


	// CE-ROT: even and odd cases use different power rates (mirror R).
	if ("`bwselect'"=="ce-rot" | "`bwselect'"=="all") {
		if (`even'==1) {
			h_ce_rot = h_mse_dpi*`n'^(-(`p'+2)/((2*`p'+5)*(`p'+3)))
			b_ce_rot = b_mse_dpi*`n'^(-(`q')  /((2*`q'+3)*(`q'+3)))
		}
		else {
			h_ce_rot = h_mse_dpi*`n'^(-(`p')/((2*`p'+3)*(`p'+3)))
			b_ce_rot = b_mse_dpi*`n'^(-(`q'+2)/((2*`q'+5)*(`q'+3)))
		}
		Result[e,2::3] = h_ce_rot, b_ce_rot
	}
	// CE-DPI even: closed-form scaling of h_mse_dpi / b_mse_dpi (mirror R).
	if ("`bwselect'"=="ce-dpi" & `even'==1) {
		h_ce_dpi = h_mse_dpi*`n'^(-(`p'+2)/((2*`p'+5)*(`p'+3)))
		b_ce_dpi = b_mse_dpi*`n'^(-(`q')  /((2*`q'+3)*(`q'+3)))
		Result[e,2::3] = h_ce_dpi, b_ce_dpi
	}



  // Full CE-DPI (odd case) uses lpbwce-equivalent computation.
  // Even case was handled above with the closed-form scaling.
  if (("`bwselect'"=="ce-dpi" | "`bwselect'"=="all") & `even'==0) {

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
		predicts = rq*beta_q

        if ("`vce_select'"=="hc2" | "`vce_select'"=="hc3") {
          hii = rowsum((Rp*invGp) :* (Rp:*Wp)) :/ eN
		}
		}


    res_q = nprobust_lp_res(eX, eY, predicts, hii, "`vce_select'", `nnmatch', dups[ind], dupsid[ind], `q'+1)

	**** Bias
	k = `p'+3
	sx = max(abs(X))
	U = X:/sx
	rk = J(N,(k+3),.)
	for (j=1; j<=(k+3); j++) {
	   rk[.,j] = U:^(j-1)
	}
	gamma = qrsolve(rk,Y)
	beta = gamma :/ (sx:^(0::(k+2)))
    mp3 = beta[`p'+4]*factorial(`p'+3) + beta[`p'+5]*factorial(`p'+4)*eval[e] + beta[`p'+6]*factorial(`p'+5)*eval[e]^2/2
	rk=rk[,1::(k+2)]
	gamma = qrsolve(rk,Y)
	beta = gamma :/ (sx:^(0::(k+1)))
	mp2 = beta[`p'+3]*factorial(`p'+2) + beta[`p'+4]*factorial(`p'+3)*eval[e] + beta[`p'+5]*factorial(`p'+4)*eval[e]^2/2

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
		phiz = normalden(1.96)  // R: dnorm(1.96) = 0.0584409443334515 (full precision)
		E = q_rbc[1]*phiz, eta_bc_1, eta_bc_2, q_rbc[3]*phiz*eta_bc_1, q_rbc[3]*phiz*eta_bc_2
		q2_rbc = q_rbc[2]

		// R uses optimize() with bracket [eps, range]; Mata's optimize()
		// finds different local minima. Replace with golden-section search
		// matching R's behavior (lpbwselect.ce.dpi npfunctions.R:657).
		h_ce_dpi = nprobust_lp_brent_cedpi(E, q2_rbc, `p', `n', 1e-12, `range')
    }

	// Scale b by the same rate as R's odd-case CE-DPI.
	b_ce_dpi =  b_mse_dpi * `n'^(-(`q'+2)/((2*`q'+5)*(`q'+3)))
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
	ereturn local vce_select = "`vce_raw'"
	ereturn local vce_type   = "`vce_type'"
	if ("`clustvar'"!="") ereturn local clustvar "`clustvar'"
	ereturn local yvar "`y'"
	ereturn local xvar "`x'"
	ereturn local cmd "lpbwselect"

	if ("`bwselect'"=="all") {
		matrix colnames Result = eval h_mse_dpi b_mse_dpi h_mse_rot b_mse_rot h_ce_dpi b_ce_dpi h_ce_rot b_ce_rot h_imse_dpi b_imse_dpi h_imse_rot b_imse_rot
	}
	else {
		matrix colnames Result = eval h b
	}
	ereturn matrix bws = Result

	mata mata clear

end


