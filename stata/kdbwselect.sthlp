{smcl}
{* *!version 1.0.0  2026-05-17}{...}
{viewerjumpto "Syntax" "kdbwselect##syntax"}{...}
{viewerjumpto "Description" "kdbwselect##description"}{...}
{viewerjumpto "Options" "kdbwselect##options"}{...}
{viewerjumpto "Examples" "kdbwselect##examples"}{...}
{viewerjumpto "Stored results" "kdbwselect##stored_results"}{...}
{viewerjumpto "References" "kdbwselect##references"}{...}
{viewerjumpto "Authors" "kdbwselect##authors"}{...}

{title:Title}

{p 4 8}{cmd:kdbwselect} {hline 2} Bandwidth Selection Procedures for Kernel Density Estimation and Inference.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:kdbwselect } {it:varname} {ifin} 
[{cmd:,} 
{cmd:eval(}{it:gridvar}{cmd:)} 
{cmd:neval(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernelfn}{cmd:)}
{cmd:bwselect(}{it:bwmethod}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:imsegrid(}{it:#}{cmd:)}
{cmd:separator(}{it:#}{cmd:)}
]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:kdbwselect} implements bandwidth selectors for kernel density point estimators and inference procedures developed in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf":Calonico, Cattaneo and Farrell (2018)}. 
See also {browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf":Calonico, Cattaneo and Farrell (2022)} for related optimality results.
It also implements other bandwidth selectors available in the literature.{p_end}

{p 4 8} A detailed introduction to this command is given in
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.

{p 4 8} Companion command is: {help kdrobust:kdrobust} for kernel density point estimation and inference procedures.{p_end}

{p 4 8}Related software useful for empirical analysis is described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/":https://nppackages.github.io/}{p_end}

{p 4 8}{it:Requires Stata 14 or later.}{p_end}


{marker options}{...}
{title:Options}

{p 4 8}{opt eval}({it:gridvar}) specifies the grid of evaluation points for {it:xvar}.
By default it uses 30 quantile-spaced points (deciles 0.1 through 0.9 in equal steps) over the support of {it:xvar}.{p_end}

{p 4 8}{opt neval}({it:#}) specifies the number of quantile-spaced evaluation points used to estimate the density. Default is 30 evaluation points. {p_end}

{p 4 8}{opt kernel}({it:kernelfn}) specifies the kernel function used. Options are: {opt epa:nechnikov}, and {opt uni:form}.
Default is {opt kernel(epanechnikov)}.{p_end}

{p 4 8}{opt bwselect}({it:bwmethod}) bandwidth selection procedure to be used.
Options are:{p_end}
{p 8 12}{opt mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default choice.{p_end}
{p 8 12}{opt imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth. {p_end}
{p 8 12}{opt imse-rot} ROT implementation of IMSE-optimal bandwidth.{p_end}
{p 8 12}{opt ce-dpi} second generation DPI implementation of CE-optimal bandwidth.{p_end}
{p 8 12}{opt ce-rot} ROT implementation of CE-optimal bandwidth.{p_end}
{p 4 12}Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb.{p_end}
{p 8 12}Default is {opt bwselect(mse-dpi)}. For details on implementation see
{browse "https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf":Calonico, Cattaneo and Farrell (2019)}.{p_end}

{p 4 8}{opt bwcheck}({it:#}) specifies an optional positive integer so that the selected bandwidth is enlarged to have at least {it:#} effective observations available at each evaluation point.{p_end}

{p 4 8}{opt imsegrid}({it:#}) number of evaluations points used to compute the IMSE bandwidth selector. Default is 30 points.{p_end}

{p 4 8}{opt separator}({it:#}) draws separator line after every {it:#} variables; default is separator(5).{p_end}


    {hline}


{marker examples}{...}
{title:Example: Cholesterol Trial Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use nprobust_data.dta}{p_end}

{p 4 8}MSE bandwidth selection procedure{p_end}
{p 8 8}{cmd:. kdbwselect chol1 if t==0}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:kdbwselect} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}sample size after listwise deletion{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}command name{p_end}
{synopt:{cmd:e(xvar)}}name of independent variable{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selection choice{p_end}
{synopt:{cmd:e(kernel)}}kernel choice{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bws)}}estimation result{p_end}
{synopt:{cmd:e(eval)}}evaluation grid{p_end}


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

