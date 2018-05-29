#' Time Series of Macroinvertabrates Abundance in the Arkansas River.
#' 
#' A time series of 16 years (5 replicates per year) of mayfly 
#' (Ephemeroptera:Heptageniidae) abundance in the fall at the monitoring station AR1
#' on the Arkansas River in Colorado, USA. 
#' 
#' @format A data frame with 90 observations on the following 2 variables.
#' \describe{
#'   \item{year}{The year of observation}
#'   \item{sqrt.mayflies}{The Square root of observed abundance.}
#' }
#' 
#' @source Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. 
#' Using SiZer to detect thresholds in ecological data. Frontiers in Ecology and 
#' the Environment 7:190-195.
#' 
#' @keywords datasets
#' 
#' @examples 
#' require(ggplot2)
#' 
#' data(Arkansas)
#' ggplot(Arkansas, aes(x=year, y=sqrt.mayflies)) + 
#'    geom_point()
#'
"Arkansas"


