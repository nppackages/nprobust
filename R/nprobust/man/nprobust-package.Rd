\name{nprobust-package}
\alias{nprobust-package}
\docType{package}
\alias{nprobust}
\title{Kernel Density and Local Polynomial Regression Methods}
\description{
   This package implements estimation, inference, bandwidth selection, and graphical procedures for kernel density and local polynomial regression methods, including robust bias-corrected confidence intervals as described in Calonico, Cattaneo and Farrell (2018): \code{\link{lprobust}} for local polynomial point estimation and robust bias-corrected inference, \code{\link{lpbwselect}} for local polynomial bandwidth selection, \code{\link{kdrobust}} for kernel density point estimation and robust bias-corrected inference, \code{\link{kdbwselect}} for kernel density bandwidth selection, and \code{\link{nprobust.plot}} for plotting results. The main methodological and numerical features are described in Calonico, Cattaneo and Farrell (2019).}
\details{
\tabular{ll}{
Package: \tab nprobust\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2026-05-17\cr
License: \tab GPL-3\cr
}
Function for LPR estimation and inference: \code{\link{lprobust}}\cr
Function for LPR bandwidth selection: \code{\link{lpbwselect}}\cr
Function for KDE estimation and inference: \code{\link{kdrobust}}\cr
Function for KDE bandwidth selection: \code{\link{kdbwselect}}\cr
Function for graphical analysis: \code{\link{nprobust.plot}}\cr
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
