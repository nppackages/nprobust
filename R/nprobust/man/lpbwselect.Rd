\name{lpbwselect}
\alias{lpbwselect}
\title{Bandwidth Selection Procedures for Local Polynomial Regression Estimation and Inference}

\description{
\code{\link{lpbwselect}} implements bandwidth selectors for local polynomial regression point estimators and inference procedures developed in Calonico, Cattaneo and Farrell (2018). See also Calonico, Cattaneo and Farrell (2022) for related optimality results. 
It also implements other bandwidth selectors available in the literature.

Companion commands: \code{\link{lprobust}} for local polynomial point estimation and inference procedures.

A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2019). For more details and related software useful for empirical analysis, visit \url{https://nppackages.github.io/}.
}
\usage{
lpbwselect(y, x, eval = NULL, neval = NULL, p = NULL, deriv = NULL,
kernel = "epa", bwselect = "mse-dpi", bwcheck = 21, bwregul = 1,
imsegrid = 30, vce = "nn", cluster = NULL,
nnmatch = 3, interior = FALSE, subset = NULL,
weights = NULL, masspoints = "check", data = NULL)
}
\arguments{
  \item{y}{dependent variable.}
  
  \item{x}{independent variable.}
  
  \item{eval}{vector of evaluation point(s). By default it uses 30 equally spaced points over the support of \code{x}.}

  \item{neval}{number of equally spaced evaluation points on the support of \code{x}. Default is \code{neval=30}.}
  
  \item{p}{polynomial order used to construct point estimator; default is \code{p = 1} (local linear regression).}
  
  \item{deriv}{derivative order of the regression function to be estimated. Default is \code{deriv=0} (regression function).}
  
  \item{kernel}{kernel function used to construct local polynomial estimators. Options are \code{epa} for the epanechnikov kernel, \code{tri} for the triangular kernel, \code{uni} for the uniform kernel and \code{gau} for the gaussian kernel. Default is \code{kernel = epa}.}
  
  \item{bwselect}{bandwidth selection procedure to be used. Options are:
  
  \code{mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default option.
	
	\code{mse-rot} ROT implementation of MSE-optimal bandwidth.
	
	\code{imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{imse-rot} ROT implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{ce-dpi} second generation DPI implementation of CE-optimal bandwidth.
	
	\code{ce-rot} ROT implementation of CE-optimal bandwidth.
	
	\code{all} reports all available bandwidth selection procedures.
	
	Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb. For details on implementation see Calonico, Cattaneo and Farrell (2019).}
	
  \item{bwcheck}{if a positive integer is provided, then the selected bandwidth is enlarged so that at least \code{bwcheck} effective observations are available at each evaluation point. Default is \code{bwcheck = 21}.}
  
  \item{bwregul}{specifies scaling factor for the regularization term added to the denominator of bandwidth selectors. Setting \code{bwregul = 0} removes the regularization term from the bandwidth selectors. Default is \code{bwregul = 1}.}
  
  \item{imsegrid}{number of evaluations points used to compute the IMSE bandwidth selector. Default is \code{imsegrid = 30}.}
  
  \item{vce}{procedure used to compute the variance-covariance matrix estimator. Options are:

  \code{nn} heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} the (minimum) number of neighbors to be used. Default choice (when \code{cluster} is not supplied).

	\code{hc0} heteroskedasticity-robust plug-in residuals variance estimator without weights.

	\code{hc1} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc1} weights.

	\code{hc2} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc2} weights.

	\code{hc3} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc3} weights.

	\code{cr1} cluster-robust variance estimator with \code{((n-1)/(n-k)) * (G/(G-1))} multiplier (Stata-style CR1). Default choice when \code{cluster} is supplied. Requires \code{cluster}.

	\code{cr2} cluster-robust variance estimator with Bell-McCaffrey block-adjusted residuals. Requires \code{cluster}.

	\code{cr3} cluster-robust variance estimator with block jackknife-style adjustment and \code{(G-1)/G} multiplier. Requires \code{cluster}.

	When \code{cluster} is supplied, \code{vce=hc0} or \code{vce=hc1} is remapped to \code{cr1}, \code{vce=hc2} to \code{cr2}, and \code{vce=hc3} to \code{cr3} (with a warning). \code{vce=nn} with \code{cluster} silently defaults to \code{cr1}.}

  \item{cluster}{indicates the cluster ID variable used for cluster-robust variance estimation. When supplied, the default \code{vce} becomes \code{cr1}; \code{cr2} and \code{cr3} are also available.}
  
  \item{nnmatch}{to be combined with for \code{vce=nn} for heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} indicating the minimum number of neighbors to be used. Default is \code{nnmatch=3}}.
  
  \item{interior}{if TRUE, all evaluation points are assumed to be interior points. This option affects only data-driven bandwidth selection via \code{\link{lpbwselect}}. Default is \code{interior = FALSE}.}
  
  \item{subset}{optional rule specifying a subset of observations to be used.}

  \item{weights}{optional vector of non-negative observation weights (multiplicative with kernel weights in all bandwidth-selection steps).}

  \item{masspoints}{how to handle evaluation points with few unique \code{x} values within bandwidth. Options: \code{"check"} (default), \code{"off"}.}

  \item{data}{an optional data frame. When supplied, \code{y}, \code{x}, \code{cluster}, \code{weights}, and \code{subset} may be given as bare variable names referring to columns of \code{data}.}
}
\value{
\item{bws}{A matrix containing \code{eval} (grid points), \code{h} and \code{b} (selected bandwidths) for each evaluation point. When \code{bwselect="all"}, the matrix is widened to include both MSE and CER bandwidths from each of the DPI and ROT estimators.}
\item{bws.imse}{IMSE-optimal bandwidths returned when \code{bwselect="all"}; \code{NULL} otherwise.}
\item{opt}{A list containing options passed to the function.}
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}. Journal of the American Statistical Association, 113(522): 767-779. \doi{10.1080/01621459.2017.1285776}.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019. \href{https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2019_JSS.pdf}{nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}. Journal of Statistical Software, 91(8): 1-33. \doi{10.18637/jss.v091.i08}.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2022. \href{https://nppackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf}{Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}. Bernoulli, 28(4): 2998-3022.
}
\author{
Sebastian Calonico, University of California, Davis, CA.  \email{scalonico@ucdavis.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{matias.d.cattaneo@gmail.com}.

Max H. Farrell, University of California, Santa Barbara, CA. \email{mhfarrell@gmail.com}.
}
\examples{
x   <- runif(500)
y   <- sin(4*x) + rnorm(500)
est <- lpbwselect(y,x)
summary(est)
}
\seealso{
  \code{\link{lprobust}}
}
