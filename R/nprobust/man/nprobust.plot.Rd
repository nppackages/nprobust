\name{nprobust.plot}
\alias{nprobust.plot}
\title{Graphical Presentation of Results from \code{nprobust} Package.}
\description{
\code{nprobust.plot} plots estimated density and regression function using the \code{nprobust} package. A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2019).

Companion commands: \code{\link{lprobust}} for local polynomial point estimation and inference procedures, and \code{\link{kdrobust}} for kernel density point estimation and inference procedures.

For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/nppackages/}.
}
\usage{
nprobust.plot(..., alpha = NULL, type = NULL, CItype = NULL,
  title = "", xlabel = "", ylabel = "", lty = NULL, lwd = NULL,
  lcol = NULL, pty = NULL, pwd = NULL, pcol = NULL, CIshade = NULL,
  CIcol = NULL, legendTitle = NULL, legendGroups = NULL)
}
\arguments{
\item{...}{Objects returned by \code{\link{kdrobust}} or \code{\link{lprobust}}.}

\item{alpha}{Numeric scalar between 0 and 1, the significance level for plotting
confidence regions. If more than one is provided, they will be applied to data series
accordingly.}

\item{type}{String, one of \code{"line"} (default), \code{"points"} or \code{"both"}, how
the point estimates are plotted. If more than one is provided, they will be applied to data series
accordingly.}

\item{CItype}{String, one of \code{"region"} (shaded region, default), \code{"line"} (dashed lines),
\code{"ebar"} (error bars), \code{"all"} (all of the previous) or \code{"none"} (no confidence region),
how the confidence region should be plotted. If more than one is provided, they will be applied to data series
accordingly.}

\item{title, xlabel, ylabel}{Strings, title of the plot and labels for x- and y-axis.}

\item{lty}{Line type for point estimates, only effective if \code{type} is \code{"line"} or
\code{"both"}. \code{1} for solid line, \code{2} for dashed line, \code{3} for dotted line.
For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
more than one is provided, they will be applied to data series accordingly.}

\item{lwd}{Line width for point estimates, only effective if \code{type} is \code{"line"} or
\code{"both"}. Should be strictly positive. For other options, see the instructions for
\code{\link{ggplot2}} or \code{\link{par}}. If more than one is provided, they will be applied
to data series accordingly.}

\item{lcol}{Line color for point estimates, only effective if \code{type} is \code{"line"} or
\code{"both"}. \code{1} for black, \code{2} for red, \code{3} for green, \code{4} for blue.
For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
more than one is provided, they will be applied to data series
accordingly.}

\item{pty}{Scatter plot type for point estimates, only effective if \code{type} is \code{"points"} or
\code{"both"}. For options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
more than one is provided, they will be applied to data series
accordingly.}

\item{pwd}{Scatter plot size for point estimates, only effective if \code{type} is \code{"points"} or
\code{"both"}. Should be strictly positive. If more than one is provided, they will be applied to data series
accordingly.}

\item{pcol}{Scatter plot color for point estimates, only effective if \code{type} is \code{"points"} or
\code{"both"}. \code{1} for black, \code{2} for red, \code{3}
for green, \code{4} for blue.
For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
more than one is provided, they will be applied to data series
accordingly.}

\item{CIshade}{Numeric, opaqueness of the confidence region, should be between 0 (transparent) and
1. Default is 0.2. If more than one is provided, they will be applied to data series
accordingly.}

\item{CIcol}{color for confidence region. \code{1} for black, \code{2} for red, \code{3}
for green, \code{4} for blue.
For other options, see the instructions for \code{\link{ggplot2}} or \code{\link{par}}. If
more than one is provided, they will be applied to data series
accordingly.}
\item{legendTitle}{String, title of legend.}
\item{legendGroups}{String Vector, group names used in legend.}
}
\value{
A standard \code{\link{ggplot2}} object is returned, hence can be used for further customization.
}
\details{
Companion command: \code{\link{lprobust}} for local polynomial-based regression functions and derivatives estimation.
}
\examples{
x   <- runif(500) 
y   <- sin(4*x) + rnorm(500)
est <- lprobust(y,x)
nprobust.plot(est)
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2019. \href{http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2019_JSS.pdf}{nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference}. Journal of Statistical Software, 91(8): 1-33. \doi{http://dx.doi.org/10.18637/jss.v091.i08}.
}
\seealso{
\code{\link{lprobust}}, \code{\link{kdrobust}}, \code{\link{ggplot2}}
}
\author{
Sebastian Calonico, Columbia University, New York, NY.  \email{sebastian.calonico@columbia.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{cattaneo@princeton.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
}
