\name{Arkansas}
\alias{Arkansas}
\docType{data}

\title{ Time Series of Macroinvertabrates Abundance in the Arkansas River.}

\description{
A time series of 16 years (5 replicates per year) of mayfly 
(Ephemeroptera:Heptageniidae) abundance in the fall at the 
monitoring station AR1. 
}

\usage{data(Arkansas)}

\format{
  A data frame with 90 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{a numeric vector}
    \item{\code{sqrt.mayflies}}{a numeric vector of abundance values}
  }
}


\source{
Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect	thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195.
}

\examples{
data(Arkansas)
plot(Arkansas)
}
\keyword{datasets}
