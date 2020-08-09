{smcl}
{* *! version 0.3.1  06-20-2020}{...}
{viewerjumpto "Syntax" "kdrobust##syntax"}{...}
{viewerjumpto "Description" "kdrobust##description"}{...}
{viewerjumpto "Options" "kdrobust##options"}{...}
{viewerjumpto "Examples" "kdrobust##examples"}{...}
{viewerjumpto "Saved results" "kdrobust##saved_results"}{...}

{title:Title}

{p 4 8}{cmd:kdrobust} {hline 2} Kernel Density Estimation with Robust Bias-Corrected Confidence Intervals and Inference Procedures.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:kdrobust} {it:varname} {ifin} 
[{cmd:,} 
{cmd:eval(}{it:gridvar}{cmd:)} 
{cmd:neval(}{it:#}{cmd:)} 
{cmd:h(}{it:#}{cmd:)} 
{cmd:b(}{it:#}{cmd:)}
{cmd:rho(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernelfn}{cmd:)}
{cmd:bwselect(}{it:bwmethod}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:imsegrid(}{it:#}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
{it:genvars}
{it:plot}
{cmd:graph_options(}{it:gphopts}{cmd:)}
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:kdrobust} implements kernel density point estimators with robust bias-corrected confidence intervals and inference procedures developed in
{browse "https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf":Calonico, Cattaneo and Farrell (2018)}. 
See also {browse "http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_CEopt.pdf":Calonico, Cattaneo and Farrell (2020)} for related optimality results. 
It also implements other estimation and inference procedures available in the literature. See Wand and Jones (1995) and Fan and Gijbels (1996) for background references.{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.

{p 4 8} Companion command is: {help kdbwselect:kdbwselect} for data-driven bandwidth selection.{p_end}

{p 4 8}Related Stata and R packages useful for empirical analysis are described in the following website:{p_end}

{p 8 8}{browse "https://sites.google.com/site/nppackages/":https://sites.google.com/site/nppackages/}{p_end}


{marker options}{...}
{title:Options}

{p 4 8}{opt eval}({it:gridvar}) specifies the grid of evaluation points for {it:xvar}.
By default it uses 30 equally spaced points over to support of {it:xvar}.{p_end}

{p 4 8}{opt neval}({it:#}) specifies the number of evaluation points to estimate the regression functions. Default is 30 evaluation points. {p_end}

{p 4 8}{opt h}({it:hvar}) specifies the main bandwidth ({it:h}) used to construct the point estimator for each evaluation point. If not specified, it is computed by the companion command {help kdbwselect:kdbwselect}.{p_end}

{p 4 8}{opt b}({it:bvar}) specifies the bias bandwidth ({it:b}) used to construct the bias-correction estimator for each evaluation point. If not specified, it is computed by the companion command {help kdbwselect:kdbwselect}.{p_end}

{p 4 8}{opt rho}({it:#}) specifies the value of {it:rho}, so that the bias bandwidth {it:b} equals {it:b}={it:h}/{it:rho}.
Default is {opt rho(1)} if {it:h} is specified but {it:b} is not.{p_end}

{p 4 8}{opt kernel}({it:kernelfn}) specifies the kernel function used to construct the kernel density estimator(s). Options are: {opt epa:nechnikov}, and {opt uni:form}.
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
{browse "http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrrell (2019)}.{p_end}

{p 4 8}{opt bwcheck}({it:#}) specifies an optional positive integer so that the selected bandwidth is enlarged to have at least {it:#} effective observations available for each evaluation point.{p_end}

{p 4 8}{opt imsegrid}({it:#}) number of evaluations points used to compute the IMSE bandwidth selector. Default is 30 points.{p_end}

{p 4 8}{opt level}({it:#}) specifies confidence level for confidence intervals.
Default is {opt level(95)}.{p_end}

{p 4 8}{opt separator}({it:#}) draws separator line after every {it:#} variables; default is separator(5).{p_end}

{p 4 8}{opt plot} generates the local polynomial regression plot.

{p 4 8}{it:genvars} generates new variables storing the following results.{p_end}
{p 8 12}{opt kdrobust_eval} evaluation points.{p_end}
{p 8 12}{opt kdrobust_h} bandwidth h.{p_end}
{p 8 12}{opt kdrobust_b} bandwidth b.{p_end}
{p 8 12}{opt kdrobust_nh} effective sample size.{p_end}
{p 8 12}{opt kdrobust_gx_us} conventional local polynomial estimate.{p_end}
{p 8 12}{opt kdrobust_se_us} conventional standard error for the local polynomial estimator.{p_end}
{p 8 12}{opt kdrobust_gx_bc} bias-corrected local polynomial regression estimate.{p_end}
{p 8 12}{opt kdrobust_se_rb} robust standard error for the local polynomial estimator.{p_end}
{p 8 12}{opt kdrobust_ci_l_rb} lower end value of the robust confidence interval.{p_end}
{p 8 12}{opt kdrobust_ci_r_rb} upper end value of the robust confidence interval.{p_end}

{p 4 8}{opt graph_options}({it:gphopts}) specifies graphical options to be passed on to the underlying graph command.{p_end}

   {hline}


{marker examples}{...}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. sysuse auto}{p_end}

{p 4 8}Kernel density estimates for length{p_end}
{p 8 8}{cmd:. kdrobust length}{p_end}

{p 4 8}Kernel density estimates for length{p_end}
{p 8 8}{cmd:. kdrobust length, plot genvars}{p_end}

{marker saved_results}{...}
{title:Saved results}

{p 4 8}{cmd:kdrobust} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}original number of observations{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(varname)}}name of variable{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection choice{p_end}
{synopt:{cmd:e(kernel)}}kernel choice{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(Result)}}estimation result{p_end}
 
 
{title:References}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018.
{browse "https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf":On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}.
{it:Journal of the American Statistical Association}, 113(522): 767-779.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019.
{browse "http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf":nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}.
{it:Journal of Statistical Software}, 91(8): 1-33. {browse "http://dx.doi.org/10.18637/jss.v091.i08":doi: 10.18637/jss.v091.i08}.{p_end}

{p 4 8}Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020.
{browse "http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2020_CEopt.pdf":Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}, working paper.{p_end}

{p 4 8}Fan, J., and Gijbels, I. 1996. Local Polynomial Modelling and Its Applications, London: Chapman and Hall.{p_end}

{p 4 8}Wand, M., and Jones, M. 1995. Kernel Smoothing, Florida: Chapman & Hall/CRC.{p_end}


{title:Authors}

{p 4 8}Sebastian Calonico, Columbia University, New York, NY.
{browse "mailto:sebastian.calonico@columbia.edu":sebastian.calonico@columbia.edu}.{p_end}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Max H. Farrell, University of Chicago, Chicago, IL.
{browse "mailto:max.farrell@chicagobooth.edu":max.farrell@chicagobooth.edu}.{p_end}

