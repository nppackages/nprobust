\name{nprobust-package}
\alias{nprobust-package}
\alias{nprobust}
\title{Nonparametric Robust Estimation and Inference Methods using Local Polynomial Regression and Kernel Density Estimation}
\description{
   This package provides tools for data-driven statistical analysis using local polynomial regression (LPR) and kernel density estimation (KDE) methods as described in Calonico, Cattaneo and Farrell (2018): \code{\link{lprobust}} for local polynomial point estimation and robust bias-corrected inference, \code{\link{lpbwselect}} for local polynomial bandwidth selection, \code{\link{kdrobust}} for kernel density point estimation and robust bias-corrected inference, \code{\link{kdbwselect}} for kernel density bandwidth selection, and \code{\link{nprobust.plot}} for plotting results. The main methodological and numerical features of this  package are described in Calonico, Cattaneo and Farrell (2019).}
\details{
\tabular{ll}{
Package: \tab nprobust\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2020-03-24\cr
License: \tab GPL-2\cr 
}
Function for LPR estimation and inference: \code{\link{lprobust}}\cr
Function for LPR bandwidth selection: \code{\link{lpbwselect}}\cr
Function for KDE estimation and inference: \code{\link{kdrobust}}\cr
Function for KDE bandwidth selection: \code{\link{kdbwselect}}\cr
Function for graphical analysis: \code{\link{nprobust.plot}}\cr
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference}. Journal of the American Statistical Association, 113(522): 767-779. \doi{doi:10.1080/01621459.2017.1285776}.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019. \href{http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf}{nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}. Journal of Statistical Software, 91(8): 1-33. \doi{http://dx.doi.org/10.18637/jss.v091.i08}.
}
\author{
Sebastian Calonico, Columbia University, New York, NY.  \email{sebastian.calonico@columbia.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{cattaneo@princeton.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
}

