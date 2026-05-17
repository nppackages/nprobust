\name{lprobust}
\alias{lprobust}
\title{Local Polynomial Methods with Robust Bias-Corrected Inference}

\description{
\code{\link{lprobust}} implements local polynomial regression point estimators, with robust bias-corrected confidence intervals and inference procedures developed in Calonico, Cattaneo and Farrell (2018). See also Calonico, Cattaneo and Farrell (2022) for related optimality results.
It also implements other estimation and inference procedures available in the literature.

Companion commands: \code{\link{lpbwselect}} for local polynomial data-driven bandwidth selection, and \code{\link{nprobust.plot}} for plotting results.

A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2019). For more details and related software useful for empirical analysis, visit \url{https://nppackages.github.io/}.
}
\usage{
lprobust(y, x, eval = NULL, neval = NULL, p = NULL, deriv = NULL,
h = NULL, b = NULL, rho = 1, kernel = "epa", bwselect = NULL,
bwcheck = 21, bwregul = 1, imsegrid = 30, vce = "nn", covgrid = FALSE,
cluster = NULL, nnmatch = 3, level = 95, interior = FALSE, subset = NULL,
weights = NULL, masspoints = "check", data = NULL)
}
\arguments{
  \item{y}{dependent variable.}
  
  \item{x}{independent variable.}
  
  \item{eval}{vector of evaluation point(s). By default it uses 30 equally spaced points over the support of \code{x}.}

  \item{neval}{number of equally spaced evaluation points on the support of \code{x}. Default is \code{neval=30}.}
  
  \item{p}{polynomial order used to construct point estimator; default is \code{p = 1} (local linear regression).}
  
  \item{deriv}{derivative order of the regression function to be estimated. Default is \code{deriv=0} (regression function).}
  
  \item{h}{main bandwidth used to construct local polynomial point estimator. Can be either scalar (same bandwidth for all evaluation points), or vector of same dimension as \code{eval}. If not specified, bandwidth \code{h} is computed by the companion command \code{\link{lpbwselect}}.}
  
  \item{b}{bias bandwidth used to construct the bias-correction estimator. Can be either scalar (same bandwidth for all evaluation points), or vector of same dimension as \code{eval}. By default it is set equal to \code{h}. If \code{rho} is set to zero, \code{b} is computed by the companion command \code{\link{lpbwselect}}.}
  
  \item{rho}{Sets \code{b=h/rho}. Default is \code{rho = 1}.}
  
  \item{kernel}{kernel function used to construct local polynomial estimators. Options are \code{epa} for the epanechnikov kernel, \code{tri} for the triangular kernel, \code{uni} for the uniform kernel and \code{gau} for the gaussian kernel. Default is \code{kernel = epa}.}
  
  \item{bwselect}{bandwidth selection procedure to be used via \code{\link{lpbwselect}}. By default it computes \code{h} and sets \code{b=h/rho} (with \code{rho=1} by default). It computes both \code{h} and \code{b} if \code{rho} is set equal to zero. Options are:
  
  \code{mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default option if only one evaluation point is chosen.
	
	\code{mse-rot} ROT implementation of MSE-optimal bandwidth.
	
	\code{imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth (computed using a grid of evaluation points). Default option if more than one evaluation point is chosen.
	
	\code{imse-rot} ROT implementation of IMSE-optimal bandwidth (computed using a grid of evaluation points).
	
	\code{ce-dpi} second generation DPI implementation of CE-optimal bandwidth.
	
	\code{ce-rot} ROT implementation of CE-optimal bandwidth.

	Use \code{\link{lpbwselect}} with \code{bwselect="all"} to report all available bandwidth selection procedures.
	
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

  \item{covgrid}{if TRUE, it computes two covariance matrices (cov.us and cov.rb) for classical and robust covariances across point estimators over the grid of evaluation points.}

  \item{cluster}{indicates the cluster ID variable used for cluster-robust variance estimation. When supplied, the default \code{vce} becomes \code{cr1}; \code{cr2} and \code{cr3} are also available.}

  \item{nnmatch}{to be combined with for \code{vce=nn} for heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} indicating the minimum number of neighbors to be used. Default is \code{nnmatch=3}}.
  
  \item{level}{confidence level used for confidence intervals; default is \code{level = 95}.}
  
  \item{interior}{if TRUE, all evaluation points are assumed to be interior points. This option affects only data-driven bandwidth selection via \code{\link{lpbwselect}}. Default is \code{interior = FALSE}.}
  
  \item{subset}{optional rule specifying a subset of observations to be used.}

  \item{weights}{optional vector of non-negative observation weights of the same length as \code{x}. User weights multiply the kernel weights in estimation, bandwidth selection, and variance computation (weighted least squares interpretation).}

  \item{masspoints}{how to handle evaluation points whose bandwidth window contains few unique \code{x} values (the local polynomial may become unreliable in that case). Options: \code{"check"} (default) warns when there are fewer than \code{p+5} unique values inside the main bandwidth; \code{"off"} disables the check.}

  \item{data}{an optional data frame. When supplied, \code{y}, \code{x}, \code{cluster}, \code{weights}, and \code{subset} may be given as bare variable names referring to columns of \code{data}.}
}
\value{
\item{Estimate}{A matrix containing \code{eval} (grid points), \code{h}, \code{b} (bandwidths), \code{N}
  (effective sample sizes), \code{tau.us} (point estimates with p-th order local polynomial),
  \code{tau.bc} (bias-corrected point estimates with (p+1)-th order local polynomial),
  \code{se.us} (standard error corresponding to \code{tau.us}), and \code{se.rb} (robust standard error).}
\item{cov.us}{Conventional-estimator covariance matrix across the evaluation grid (\code{neval} by \code{neval}). Returned only when \code{covgrid=TRUE}; \code{NULL} otherwise.}
\item{cov.rb}{Robust-bias-corrected covariance matrix across the evaluation grid (\code{neval} by \code{neval}). Returned only when \code{covgrid=TRUE}; \code{NULL} otherwise.}
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
est <- lprobust(y,x)
summary(est)
}
\keyword{ LPR }
\keyword{ Robust Estimation }
\seealso{
  \code{\link{lpbwselect}}
}
