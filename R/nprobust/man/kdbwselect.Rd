\name{kdbwselect}
\alias{kdbwselect}
\title{Bandwidth Selection Procedures for Kernel Density Estimation and Inference}

\description{
\code{\link{kdbwselect}} implements bandwidth selectors for kernel density point estimators and inference procedures developed in Calonico, Cattaneo and Farrell (2018). See also Calonico, Cattaneo and Farrell (2019b) for related optimality results.
It also implements other bandwidth selectors available in the literature. See Wand and Jones (1995) for background references.

Companion commands are: \code{\link{kdrobust}} for kernel density point estimation and inference procedures.

A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2019a). For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/nppackages/}.
}
\usage{
kdbwselect(x, eval = NULL, neval = NULL, kernel = "epa", 
bwselect = "mse-dpi", bwcheck=21, imsegrid=30, subset = NULL)
}
\arguments{
  \item{x}{independent variable.}
  
  \item{eval}{vector of evaluation point(s). By default it uses 30 equally spaced points over to support of \code{x}.}
  
  \item{neval}{number of quantile-spaced evaluation points on support of \code{x}. Default is \code{neval=30}.}
  
  \item{kernel}{kernel function used to construct the kernel estimators. Options are \code{epa} for the epanechnikov kernel, and \code{uni} for the uniform kernel. Default is \code{kernel = epa}.}
  
  \item{bwselect}{bandwidth selection procedure to be used. Options are:
  
  \code{mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default option.
	
	\code{imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{imse-rot} ROT implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{ce-dpi} second generation DPI implementation of CE-optimal bandwidth.
	
	\code{ce-rot} ROT implementation of CE-optimal bandwidth.
	
	\code{all} reports all available bandwidth selection procedures.
	
	Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb. For details on implementation see Calonico, Cattaneo and Farrell (2019a).}
	
  \item{bwcheck}{if a positive integer is provided, then the selected bandwidth is enlarged so that at least \code{bwcheck} effective observations are available at each evaluation point. Default is \code{bwcheck = 15}.}
  
  \item{imsegrid}{number of evaluations points used to compute the IMSE bandwidth selector. Default is \code{imsegrid = 30}.}
  
  \item{subset}{optional rule specifying a subset of observations to be used.}
}
\value{
\item{Estimate}{A matrix containing \code{eval} (grid points), \code{h} and \code{b} (bandwidths).}
\item{opt}{A list containing options passed to the function.}
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}. Journal of the American Statistical Association, 113(522): 767-779. \doi{doi:10.1080/01621459.2017.1285776}.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019a. \href{http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf}{nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}. Journal of Statistical Software, 91(8). \doi{doi:10.18637/jss.v091.i08}.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019b. \href{https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_CEopt.pdf}{Coverage Error Optimal Confidence Intervals for Local Polynomial Regression}. Working Paper. 

Fan, J., and Gijbels, I. 1996. Local polynomial modelling and its applications, London: Chapman and Hall.

Wand, M., and Jones, M. 1995. Kernel Smoothing, Florida: Chapman & Hall/CRC.
}
\author{
Sebastian Calonico, Columbia University, New York, NY.  \email{sebastian.calonico@columbia.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{cattaneo@princeton.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
}
\examples{
x   <- rnorm(500)
est <- kdbwselect(x)
summary(est)
}
\seealso{
  \code{\link{kdrobust}}
}
