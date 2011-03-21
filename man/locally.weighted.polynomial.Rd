\name{locally.weighted.polynomial}
\alias{locally.weighted.polynomial}
\alias{calc.CI.LocallyWeightedPolynomial}
\alias{find.state}
\alias{find.state.changes}
\alias{find.states}
\alias{kernel.h}
\alias{positive.part}
\alias{x.grid.create}

\title{ Locally-Weighted Polynomial Regression Smoother}

\description{
Smoothes the given bivariate data using kernel regression.
}

\usage{
locally.weighted.polynomial(x, y, h = NA, x.grid = NA, 
	degree = 1, kernel.type = "Normal")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of data for the independent variable}
  \item{y}{Vector of data for the dependent variable}
  \item{h}{The bandwidth for the kernel}
  \item{x.grid}{What x-values should the value of the smoother be calculated at.}
  \item{degree}{The degree of the polynomial to be fit at each x-value. The default is
	to fit a linear regression, ie degree=1.}
  \item{kernel.type}{What kernel to use.  Valid choices are 'Normal',
        'Epanechnikov', 'biweight', and 'triweight'}
}

\details{
The confidence intervals are created using the row-wise method of Hannig and Marron (2006).

Notice that the derivative to be estimated must be less than or equal to the degree of the polynomial
initially fit to the data.

If the bandwidth is not given, the Sheather-Jones bandwidth selection method is used.
}

\value{
Returns a \code{LocallyWeightedPolynomial} object that has the following elements:
\item{data}{A structure of the data used to generate the smoothing curve}
\item{h}{The bandwidth used to generate the smoothing curve.}
\item{x.grid}{The grid of x-values that we have estimated function value
          and derivative(s) for.}
\item{degrees.freedom}{The effective sample size at each grid point}
\item{Beta}{A matrix of estimated beta values.  The number of rows is
          degrees+1, while the number of columns is the same as the length
          of x.grid. Notice that 
		\deqn{ \hat{f}(x_i)   = \beta[1,i] }
        \deqn{ \hat{f'}(x_i)  = \beta[2,i]*1! }
        \deqn{ \hat{f''}(x_i) = \beta[3,i]*2! }  
and so on...}
\item{Beta.var}{Matrix of estimated variances for \code{Beta}.  Same structure as \code{Beta}.}
}

\author{Derek Sonderegger}
\references{
Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures 
	in curves. Journal of the American Statistical Association 94 807-823. 

Hannig, J., and J. S. Marron. 2006. Advanced distribution theory for SiZer. 
	Journal of the American Statistical Association 101 484-499.

Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect	thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195
}

\seealso{\code{\link{SiZer}}, \code{\link{plot.LocallyWeightedPolynomial}},  \code{spm} in package 'SemiPar',
	\code{\link{loess}}, \code{\link{smooth.spline}}, \code{\link{interpSpline}}}

\examples{
data(Arkansas)
x <- Arkansas$year
y <- Arkansas$sqrt.mayflies

layout(cbind(1,2,3))
model <- locally.weighted.polynomial(x,y)
plot(model, main='Smoothed Function', xlab='Year', ylab='Sqrt.Mayflies')

model2 <- locally.weighted.polynomial(x,y,h=.5)
plot(model2, main='Smoothed Function', xlab='Year', ylab='Sqrt.Mayflies')

model3 <- locally.weighted.polynomial(x,y, degree=1)
plot(model3, derv=1, main='First Derivative', xlab='Year', ylab='1st Derivative')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
