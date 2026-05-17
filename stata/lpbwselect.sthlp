{smcl}
{* *!version 1.0.0  2026-05-17}{...}
{viewerjumpto "Syntax" "lpbwselect##syntax"}{...}
{viewerjumpto "Description" "lpbwselect##description"}{...}
{viewerjumpto "Options" "lpbwselect##options"}{...}
{viewerjumpto "Examples" "lpbwselect##examples"}{...}
{viewerjumpto "Stored results" "lpbwselect##stored_results"}{...}
{viewerjumpto "References" "lpbwselect##references"}{...}
{viewerjumpto "Authors" "lpbwselect##authors"}{...}

{title:Title}

{p 4 8}{cmd:lpbwselect} {hline 2} Bandwidth Selection Procedures for Local Polynomial Regression Estimation and Inference.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:lpbwselect } {it:yvar} {it:xvar} {ifin} 
[{cmd:,} 
{cmd:eval(}{it:gridvar}{cmd:)} 
{cmd:neval(}{it:#}{cmd:)}
{cmd:deriv(}{it:#}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernelfn}{cmd:)}
{cmd:bwselect(}{it:bwmethod}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:imsegrid(}{it:#}{cmd:)}
{cmd:vce(}{it:vcetype [vceopt]}{cmd:)}
{cmd:bwregul(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
{it:interior}
{it:genvars}
{cmd:weights(}{it:wvar}{cmd:)}
{cmd:masspoints(}{it:check|off}{cmd:)}
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:lpbwselect} implements bandwidth selectors for local polynomial regression point estimators and inference procedures developed in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf":Calonico, Cattaneo and Farrell (2018)}. 
See also {browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf":Calonico, Cattaneo and Farrell (2022)} for related optimality results.
It also implements other bandwidth selectors available in the literature.{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.

{p 4 8} Companion command is: {help lprobust:lprobust} for local polynomial point estimation and inference procedures.{p_end}

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

{p 4 8}{opt kernel}({it:kernelfn}) specifies the kernel function used to construct the local-polynomial estimator(s). Options are: {opt tri:angular}, {opt epa:nechnikov}, {opt uni:form} and {opt gau:ssian}.
Default is {opt kernel(epanechnikov)}.{p_end}

{p 4 8}{opt bwselect}({it:bwmethod}) bandwidth selection procedure to be used. 
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
{p 8 12}{cmd:vce(cr1 }{it:clustervar}{cmd:)} cluster-robust CR1 variance with {it:((n-1)/(n-k)) * (G/(G-1))} multiplier. Default when a cluster variable is supplied.{p_end}
{p 8 12}{cmd:vce(cr2 }{it:clustervar}{cmd:)} cluster-robust CR2 variance (Bell-McCaffrey block-adjusted residuals).{p_end}
{p 8 12}{cmd:vce(cr3 }{it:clustervar}{cmd:)} cluster-robust CR3 variance (block jackknife, {it:(G-1)/G} multiplier).{p_end}
{p 8 12}{cmd:vce(cluster }{it:clustervar}{cmd:)} alias for {cmd:vce(cr1 }{it:clustervar}{cmd:)}.{p_end}
{p 8 12}Legacy: {cmd:vce(hc0/hc1 }{it:clustervar}{cmd:)} is remapped to {cmd:cr1}, {cmd:vce(hc2 }{it:clustervar}{cmd:)} to {cmd:cr2}, {cmd:vce(hc3 }{it:clustervar}{cmd:)} to {cmd:cr3}. {cmd:vce(nncluster ...)} is no longer supported and is remapped to {cmd:cr1} with a warning.{p_end}
{p 8 12}Default is {opt vce(nn 3)} when no cluster, {opt vce(cr1)} when a cluster variable is supplied.{p_end}

{p 4 8}{opt bwregul}({it:#}) specifies scaling factor for the regularization term added to the denominator of the bandwidth selectors. Setting {opt bwregul(0)} removes the regularization term from the bandwidth selectors.
Default is {opt bwregul(1)}.{p_end}

{p 4 8}{opt separator}({it:#}) draws separator line after every {it:#} variables; default is separator(5).{p_end}

{p 4 8}{opt interior} option to set all evaluation points to be interior points.{p_end}

{p 4 8}{it:genvars} generates new variables {opt lpbwselect_eval}, {opt lpbwselect_h}, and {opt lpbwselect_b} storing the evaluation points and selected bandwidths (or per-method bandwidths when {opt bwselect(all)}).{p_end}

{p 4 8}{opt weights}({it:wvar}) optional vector of non-negative observation weights (multiplicative with kernel weights in all bandwidth-selection steps).{p_end}

{p 4 8}{opt masspoints}({it:check|off}) how to handle evaluation points with few unique {it:xvar} values within bandwidth. Default is {opt check}; {opt off} disables the warning.{p_end}


    {hline}


{marker examples}{...}
{title:Example: Cholesterol Trial Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use nprobust_data.dta}{p_end}

{p 4 8}Second-generation DPI implementation of MSE-optimal bandwidth{p_end}
{p 8 8}{cmd:. lpbwselect cholf chol1 if t==0}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:lpbwselect} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size after listwise deletion{p_end}
{synopt:{cmd:e(p)}}order of the polynomial used for estimation of the regression function{p_end}
{synopt:{cmd:e(q)}}order of the polynomial used for bias correction{p_end}

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
{synopt:{cmd:e(bws)}}estimation result{p_end}


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

