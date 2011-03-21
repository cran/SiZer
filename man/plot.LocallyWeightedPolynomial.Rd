\name{plot.LocallyWeightedPolynomial}
\alias{plot.LocallyWeightedPolynomial}

\title{Plot a LocallyWeightedPolynomial object}
\description{
Creates a plot of an object created by \code{locally.weighted.polynomial}.
}
\usage{
\method{plot}{LocallyWeightedPolynomial}(x, derv = 0, CI.method = 2, 
	alpha = 0.05, use.ess = TRUE, draw.points = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{LocallyWeightedPolynomial object}
  \item{derv}{Derivative to be plotted.  Default is 0 - which plots the smoothed function.}
  \item{CI.method}{What method should be used to calculate the confidence interval about the estimated line.  The
		methods are from Hannig and Marron (2006), where 1 is the point-wise estimate, and 2 is the row-wise estimate.}
  \item{alpha}{The CI has a 1-alpha/2 level of significance.}
  \item{use.ess}{ESS stands for the estimated sample size.  If at any point along the x-axis, the ESS is too small, then
 		we will not plot unless use.ess=FALSE.}
  \item{draw.points}{Should the data points be included in the graph?}
  \item{\dots}{Additional arguments to be passed to the graphing functions.}
}

\references{
Hannig, J., and J. S. Marron. 2006. Advanced distribution theory for SiZer. 
	Journal of the American Statistical Association 101:484-499.

Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect	thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195.
}

\author{Derek Sonderegger}

\seealso{ \code{\link{locally.weighted.polynomial}} }
\examples{
data('Arkansas')
x <- Arkansas$year
y <- Arkansas$sqrt.mayflies

model <- locally.weighted.polynomial(x,y)
plot(model)

model <- locally.weighted.polynomial(x,y,degree=2)
plot(model, derv=1)
plot(model, derv=2)
}

\keyword{hplot}
