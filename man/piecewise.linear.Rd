\name{piecewise.linear}
\alias{piecewise.linear}
\alias{piecewise.linear.simple}
\alias{logLik.PiecewiseLinear}
\alias{predict.PiecewiseLinear}
\alias{plot.PiecewiseLinear}
\alias{print.PiecewiseLinear}
\alias{piecewise.linear.likelihood}

\title{Fit a piecewise linear model}

\description{
Fit a degree 1 spline with 1 knot point where the location of the 
knot point is unknown.
}

\usage{
piecewise.linear(x, y, middle = 1, CI = FALSE, 
	bootstrap.samples = 1000, sig.level = 0.05)
}

\arguments{
  \item{x}{ Vector of data for the x-axis.}
  \item{y}{ Vector of data for the y-axis }
  \item{middle}{A scalar in \eqn{[0,1]}.  This represents the range that the change-point
 	can occur in.  \eqn{0} means the change-point must occur at the middle of the range of 
	x-values.  \eqn{1} means that the change-point can occur anywhere along the range of the 
	x-values.}
  \item{CI}{Whether or not a bootstrap confidence interval should be calculated.}
  \item{bootstrap.samples}{The number of bootstrap samples to take.}
  \item{sig.level}{What significance level to use for the confidence intervals.}
}

\details{
The bootstrap samples are taken by resampling the raw data points.  Often a more appropriate
bootstrap sample would be to calculate the residuals and then add a randomly selected residual
to each y-value.
}

\value{
A list of 5 elements is returned:
	\item{change.point}{The estimate of \eqn{\alpha}.}
	\item{model}{The resulting \code{lm} object once \eqn{\alpha} is known. }
	\item{x}{The x-values used.}
	\item{y}{The y-values used.}
	\item{CI}{Whether or not the confidence interval was calculated.}
	\item{intervals}{If the CIs where calculated, this is a matrix of the
		upper and lower intervals.}
  ...
}
\references{
Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression 
	theory and applications. Journal of the American Statistical Association 
	101:542-553.

Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for 
		identifying ecological thresholds. Ecology 84:2034-2041. 
}	

\author{Derek Sonderegger}

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(Arkansas)
x <- Arkansas$year
y <- Arkansas$sqrt.mayflies

model <- piecewise.linear(x,y, CI=FALSE)
plot(model)
print(model)
predict(model, 2001)
}

\keyword{regression}
