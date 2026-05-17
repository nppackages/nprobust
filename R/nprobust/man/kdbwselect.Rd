\name{kdbwselect}
\alias{kdbwselect}
\title{Bandwidth Selection Procedures for Kernel Density Estimation and Inference}

\description{
\code{\link{kdbwselect}} implements bandwidth selectors for kernel density point estimators and inference procedures developed in Calonico, Cattaneo and Farrell (2018). See also Calonico, Cattaneo and Farrell (2022) for related optimality results.
It also implements other bandwidth selectors available in the literature.

Companion commands are: \code{\link{kdrobust}} for kernel density point estimation and inference procedures.

A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2019). For more details and related software useful for empirical analysis, visit \url{https://nppackages.github.io/}.
}
\usage{
kdbwselect(x, eval = NULL, neval = NULL, kernel = "epa",
bwselect = "mse-dpi", bwcheck=21, imsegrid=30, subset = NULL, data = NULL)
}
\arguments{
  \item{x}{independent variable.}
  
  \item{eval}{vector of evaluation point(s). By default it uses 30 quantile-spaced points (deciles 0.1 to 0.9 in equal steps) over the support of \code{x}.}

  \item{neval}{number of quantile-spaced evaluation points on the support of \code{x}. Default is \code{neval=30}.}
  
  \item{kernel}{kernel function used to construct the kernel estimators. Options are \code{epa} for the epanechnikov kernel, and \code{uni} for the uniform kernel. Default is \code{kernel = epa}.}
  
  \item{bwselect}{bandwidth selection procedure to be used. Options are:
  
  \code{mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default option.
	
	\code{imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{imse-rot} ROT implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{ce-dpi} second generation DPI implementation of CE-optimal bandwidth.
	
	\code{ce-rot} ROT implementation of CE-optimal bandwidth.
	
	\code{all} reports all available bandwidth selection procedures.
	
	Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb. For details on implementation see Calonico, Cattaneo and Farrell (2019).}
	
  \item{bwcheck}{if a positive integer is provided, then the selected bandwidth is enlarged so that at least \code{bwcheck} effective observations are available at each evaluation point. Default is \code{bwcheck = 21}.}
  
  \item{imsegrid}{number of evaluations points used to compute the IMSE bandwidth selector. Default is \code{imsegrid = 30}.}
  
  \item{subset}{optional rule specifying a subset of observations to be used.}

  \item{data}{an optional data frame. When supplied, \code{x} and \code{subset} may be given as bare variable names referring to columns of \code{data}.}
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
x   <- rnorm(500)
est <- kdbwselect(x)
summary(est)
}
\seealso{
  \code{\link{kdrobust}}
}
