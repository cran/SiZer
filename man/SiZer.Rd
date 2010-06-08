\name{SiZer}

\alias{SiZer}
\alias{h.grid.create}


\title{Calculate SiZer Map}

\description{Calculates the SiZer map from a given set of X and Y variables.}

\usage{ SiZer(x, y, h=NA, x.grid=NA, degree=NA, derv=1, grid.length=41) }

\arguments{
	\item{x}{data vector for the independent axis}
	\item{y}{data vector for the dependent axis}
	\item{h}{An integer representing how many bandwidths should be considered, or
			vector of length 2 representing the upper and lower limits h should take,
			or a vector of length greater than two indicating which bandwidths to examine.}
	\item{x.grid}{An integer representing how many bins to use along the x-axis, or
			a vector of length 2 representing the upper and lower limits the x-axis
			should take, or a vector of length greater than two indicating which 
			x-values the derivative should be evaluated at.}
	\item{grid.length}{The default length of the \code{h.grid} or \code{x.grid} if the length of
			either is not given.}
	\item{derv}{The order of derivative for which to make the SiZer map.}
	\item{degree}{The degree of the local weighted polynomial used to smooth the data.
			This must be greater than or equal to \code{derv}.}
}

\details{
SiZer stands for the Significant Zero crossings of the derivative.  There are two dominate
approaches in smoothing bivariate data: locally weighted regression or penalized splines.
Both approaches require the use of a 'bandwidth' parameter that controls how much smoothing
should be done.  Unfortunately there is no uniformly best bandwidth selection procedure.
SiZer (Chaudhuri and Marron, 1999) is a procedure that looks across a range of bandwidths
and classifies the p-th derivative of the smoother into one of three states: significantly
increasing (blue), possibly zero (purple), or significantly negative (red).
}

\value{
Returns an SiZer object which has the following components:
	\item{x.grid}{Vector of x-values at which the derivative was evaluated.}
	\item{h.grid}{Vector of bandwidth values for which a smoothing function was calculated.} 
	\item{slopes}{Matrix of what category a particular x-value and bandwidth falls into 
	(Increasing=1, Possibly Zero=0, Decreasing=-1, Not Enough Data=2).}
}

\references{
Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures 
	in curves. Journal of the American Statistical Association 94:807-823. 

Hannig, J., and J. S. Marron. 2006. Advanced distribution theory for SiZer. 
	Journal of the American Statistical Association 101:484-499.

Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect	thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195.
}

\author{Derek Sonderegger}

\seealso{\code{\link{plot.SiZer}}, \code{\link{locally.weighted.polynomial}}}

\examples{
data('Arkansas')
x <- Arkansas$year
y <- Arkansas$sqrt.mayflies

plot(x,y)

# Calculate the SiZer map for the first derivative
SiZer.1 <- SiZer(x, y, h=c(.5,10), degree=1, derv=1)
plot(SiZer.1)

# Calculate the SiZer map for the second derivative
SiZer.2 <- SiZer(x, y, h=c(.5,10), degree=2, derv=2);
plot(SiZer.2)

# By setting the grid.length larger, we get a more detailed SiZer
# map but it takes longer to compute. 
#
# SiZer.3 <- SiZer(x, y, h=c(.5,10), grid.length=100, degree=1, derv=1)
# plot(SiZer.3)
}

\keyword{hplot}
\keyword{smooth}
