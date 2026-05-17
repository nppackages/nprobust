{smcl}
{* *!version 1.0.0  2026-05-17}{...}
{viewerjumpto "Syntax" "lprobust##syntax"}{...}
{viewerjumpto "Description" "lprobust##description"}{...}
{viewerjumpto "Options" "lprobust##options"}{...}
{viewerjumpto "Examples" "lprobust##examples"}{...}
{viewerjumpto "Stored results" "lprobust##stored_results"}{...}
{viewerjumpto "References" "lprobust##references"}{...}
{viewerjumpto "Authors" "lprobust##authors"}{...}

{title:Title}

{p 4 8}{cmd:lprobust} {hline 2} Local Polynomial Regression Estimation with Robust Bias-Corrected Confidence Intervals and Inference Procedures.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lprobust} {it:yvar} {it:xvar} {ifin} 
[{cmd:,} 
{cmd:eval(}{it:gridvar}{cmd:)} 
{cmd:neval(}{it:#}{cmd:)}
{cmd:deriv(}{it:#}{cmd:)}
{cmd:p(}{it:#}{cmd:)} 
{cmd:h(}{it:hvar}{cmd:)} 
{cmd:b(}{it:bvar}{cmd:)}
{cmd:rho(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernelfn}{cmd:)}
{cmd:bwselect(}{it:bwmethod}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:imsegrid(}{it:#}{cmd:)}
{cmd:vce(}{it:vcetype [vceopt]}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:bwregul(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
{it:interior}
{it:genvars}
{it:covgrid}
{it:plot}
{cmd:weights(}{it:wvar}{cmd:)}
{cmd:masspoints(}{it:check|off}{cmd:)}
{cmd:graph_options(}{it:gphopts}{cmd:)}
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:lprobust} implements local polynomial regression point estimators with robust bias-corrected confidence intervals and inference procedures developed in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf":Calonico, Cattaneo and Farrell (2018)}. 
See also {browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf":Calonico, Cattaneo and Farrell (2022)} for related optimality results.
It also implements other estimation and inference procedures available in the literature.{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.

{p 4 8} Companion command is: {help lpbwselect:lpbwselect} for data-driven bandwidth selection.{p_end}

{p 4 8}Related software useful for empirical analysis is described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}

{p 4 8}{it:Requires Stata 14 or later.}{p_end}


{marker options}{...}
{title:Options}

{p 4 8}{opt eval}({it:gridvar}) specifies the grid of evaluation points for {it:xvar}.
By default it uses 30 equally spaced points over the support of {it:xvar}.{p_end}

{p 4 8}{opt neval}({it:#}) specifies the number of evaluation points to estimate the regression functions. Default is 30 evaluation points. {p_end}

{p 4 8}{opt deriv}({it:#}) specifies the order of the derivative of the regression functions to be estimated.
Default is {opt deriv(0)}.{p_end}

{p 4 8}{opt p}({it:#}) specifies the order of the local polynomial used to construct the point estimator.
Default is {opt p(1)} (local linear regression).{p_end}

{p 4 8}{opt h}({it:hvar}) specifies the main bandwidth ({it:h}) used to construct the point estimator for each evaluation point. If not specified, it is computed by the companion command {help lpbwselect:lpbwselect}.{p_end}

{p 4 8}{opt b}({it:bvar}) specifies the bias bandwidth ({it:b}) used to construct the bias-correction estimator for each evaluation point. If not specified, it is computed by the companion command {help lpbwselect:lpbwselect}.{p_end}

{p 4 8}{opt rho}({it:#}) specifies the value of {it:rho}, so that the bias bandwidth {it:b} equals {it:b}={it:h}/{it:rho}.
Default is {opt rho(1)} if {it:h} is specified but {it:b} is not.{p_end}

{p 4 8}{opt kernel}({it:kernelfn}) specifies the kernel function used to construct the local-polynomial estimator(s). Options are: {opt tri:angular}, {opt epa:nechnikov}, {opt uni:form} and {opt gau:ssian}.
Default is {opt kernel(epanechnikov)}.{p_end}

{p 4 8}{opt bwselect}({it:bwmethod}) bandwidth selection procedure to be used. By default it computes both h and b, unless rho is specified, in which case it only computes h and sets b=h/rho.
Options are:{p_end}
{p 8 12}{opt mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default choice.{p_end}
{p 8 12}{opt mse-rot} ROT implementation of MSE-optimal bandwidth. {p_end}
{p 8 12}{opt imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth. {p_end}
{p 8 12}{opt imse-rot} ROT implementation of IMSE-optimal bandwidth.{p_end}
{p 8 12}{opt ce-dpi} second generation DPI implementation of CE-optimal bandwidth.{p_end}
{p 8 12}{opt ce-rot} ROT implementation of CE-optimal bandwidth.{p_end}
{p 4 12}Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb.{p_end}
{p 8 12}Default is {opt bwselect(mse-dpi)}. For details on implementation see
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.{p_end}

{p 4 8}{opt bwcheck}({it:#}) specifies an optional positive integer so that the selected bandwidth is enlarged to have at least {it:#} effective observations available for each evaluation point.{p_end}

{p 4 8}{opt imsegrid}({it:#}) number of evaluations points used to compute the IMSE bandwidth selector. Default is 30 points.{p_end}

{p 4 8}{opt vce}({it:vcetype [vceopt1]}) specifies the procedure used to compute the variance-covariance matrix estimator.
Options are:{p_end}
{p 8 12}{opt vce}({it:nn [nnmatch]}) heteroskedasticity-robust nearest neighbor variance estimator with {it:nnmatch} minimum neighbors. Default when no cluster variable is supplied.{p_end}
{p 8 12}{opt vce(hc0)} heteroskedasticity-robust plug-in residuals, no weights.{p_end}
{p 8 12}{opt vce(hc1)} heteroskedasticity-robust plug-in residuals, HC1 weights.{p_end}
{p 8 12}{opt vce(hc2)} heteroskedasticity-robust plug-in residuals, HC2 weights (leverage).{p_end}
{p 8 12}{opt vce(hc3)} heteroskedasticity-robust plug-in residuals, HC3 weights (jackknife).{p_end}
{p 8 12}{cmd:vce(cr1 }{it:clustervar}{cmd:)} cluster-robust CR1 variance (raw residuals, {it:((n-1)/(n-k)) * (G/(G-1))} multiplier). Default when a cluster variable is supplied.{p_end}
{p 8 12}{cmd:vce(cr2 }{it:clustervar}{cmd:)} cluster-robust CR2 variance (Bell-McCaffrey block-adjusted residuals). {p_end}
{p 8 12}{cmd:vce(cr3 }{it:clustervar}{cmd:)} cluster-robust CR3 variance (block jackknife, {it:(G-1)/G} multiplier). {p_end}
{p 8 12}{cmd:vce(cluster }{it:clustervar}{cmd:)} alias for {cmd:vce(cr1 }{it:clustervar}{cmd:)} (CR1). {p_end}
{p 8 12}Legacy: {cmd:vce(hc0/hc1 }{it:clustervar}{cmd:)} is remapped to {cmd:cr1}, {cmd:vce(hc2 }{it:clustervar}{cmd:)} to {cmd:cr2}, {cmd:vce(hc3 }{it:clustervar}{cmd:)} to {cmd:cr3}. {cmd:vce(nncluster ...)} is no longer supported and is remapped to {cmd:cr1} with a warning.{p_end}
{p 8 12}Default is {opt vce(nn 3)} when no cluster, {opt vce(cr1)} when a cluster variable is supplied.{p_end}

{p 4 8}{opt level}({it:#}) specifies confidence level for confidence intervals.
Default is {opt level(95)}.{p_end}

{p 4 8}{opt bwregul}({it:#}) specifies scaling factor for the regularization term added to the denominator of the bandwidth selectors. Setting {opt bwregul(0)} removes the regularization term from the bandwidth selectors.
Default is {opt bwregul(1)}.{p_end}

{p 4 8}{opt separator}({it:#}) draws separator line after every {it:#} variables; default is separator(5).{p_end}

{p 4 8}{opt interior} option to set all evaluation points to be interior points. This option affects only data-driven bandwidth selection via {help lpbwselect:lpbwselect}.{p_end}

{p 4 8}{opt covgrid} option to compute two covariance matrices (cov_us and cov_rb) for classical and robust covariances across point estimators over the grid of evaluation points.{p_end}

{p 4 8}{opt weights}({it:wvar}) optional vector of non-negative observation weights. User weights multiply the kernel weights in estimation, bandwidth selection, and variance computation (weighted least squares interpretation).{p_end}

{p 4 8}{opt masspoints}({it:check|off}) handles evaluation points whose bandwidth window contains few unique values of {it:xvar}. Options:{p_end}
{p 8 12}{opt check} (default) warns when fewer than {it:p+5} unique x values fall inside the bandwidth window.{p_end}
{p 8 12}{opt off} disables the check.{p_end}

{p 4 8}{opt plot} generates the local polynomial regression plot.

{p 4 8}{it:genvars} generates new variables storing the following results.{p_end}
{p 8 12}{opt lprobust_eval} evaluation points.{p_end}
{p 8 12}{opt lprobust_h} bandwidth h.{p_end}
{p 8 12}{opt lprobust_b} bandwidth b.{p_end}
{p 8 12}{opt lprobust_N} effective sample size.{p_end}
{p 8 12}{opt lprobust_tau_us} conventional local polynomial estimate.{p_end}
{p 8 12}{opt lprobust_se_us} conventional standard error for the local polynomial estimator.{p_end}
{p 8 12}{opt lprobust_tau_bc} bias-corrected local polynomial regression estimate.{p_end}
{p 8 12}{opt lprobust_se_rb} robust standard error for the local polynomial estimator.{p_end}
{p 8 12}{opt lprobust_ci_l_rb} lower end value of the robust confidence interval.{p_end}
{p 8 12}{opt lprobust_ci_r_rb} upper end value of the robust confidence interval.{p_end}

{p 4 8}{opt graph_options}({it:gphopts}) specifies graphical options to be passed on to the underlying graph command.{p_end}

   {hline}


{marker examples}{...}
{title:Example: Cholesterol Trial Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use nprobust_data.dta}{p_end}

{p 4 8}Local linear regression with second-generation DPI implementation of MSE-optimal bandwidth{p_end}
{p 8 8}{cmd:. lprobust cholf chol1 if t==0}{p_end}

{p 4 8}Same as above, but generating a plot and the corresponding output variables{p_end}
{p 8 8}{cmd:. lprobust cholf chol1 if t==0, plot genvars}{p_end}



{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:lprobust} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size after listwise deletion{p_end}
{synopt:{cmd:e(p)}}order of the polynomial used for estimation of the regression function{p_end}
{synopt:{cmd:e(q)}}order of the polynomial used for bias correction{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}command name{p_end}
{synopt:{cmd:e(yvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(xvar)}}name of running variable{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection choice{p_end}
{synopt:{cmd:e(kernel)}}kernel choice{p_end}
{synopt:{cmd:e(vce_select)}}vce choice (post-remap, lowercase){p_end}
{synopt:{cmd:e(vce_type)}}vce display label{p_end}
{synopt:{cmd:e(clustvar)}}cluster variable (if specified){p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(Result)}}estimation result{p_end}
{synopt:{cmd:e(cov_us)}}covariance matrix over evaluation grid (point estimator), if {cmd:covgrid} is specified{p_end}
{synopt:{cmd:e(cov_rb)}}covariance matrix over evaluation grid (robust bias-corrected), if {cmd:covgrid} is specified{p_end}
 
 
{marker references}{...}
{title:References}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018.
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf":On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}.
{it:Journal of the American Statistical Association}, 113(522): 767-779.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019.
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}.
{it:Journal of Statistical Software}, 91(8): 1-33. {browse "https://doi.org/10.18637/jss.v091.i08":doi: 10.18637/jss.v091.i08}.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2022.
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf":Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}. {it:Bernoulli}, 28(4): 2998-3022.{p_end}


{marker authors}{...}
{title:Authors}

{p 4 8}Sebastian Calonico, University of California, Davis, CA.
{browse "mailto:scalonico@ucdavis.edu":scalonico@ucdavis.edu}.{p_end}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:matias.d.cattaneo@gmail.com":matias.d.cattaneo@gmail.com}.{p_end}

{p 4 8}Max H. Farrell, University of California, Santa Barbara, CA.
{browse "mailto:mhfarrell@gmail.com":mhfarrell@gmail.com}.{p_end}
