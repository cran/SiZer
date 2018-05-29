#' Calculate SiZer Map
#' 
#' Calculates the SiZer map from a given set of X and Y variables.
#' 
#' @usage SiZer(x, y, h=NA, x.grid=NA, degree=NA, derv=1, grid.length=41)
#' 
#' @param x data vector for the independent axis
#' @param y data vector for the dependent axis
#' @param h An integer representing how many bandwidths should be considered, or
#'   vector of length 2 representing the upper and lower limits h should take,
#'   or a vector of length greater than two indicating which bandwidths to examine.
#' @param x.grid An integer representing how many bins to use along the x-axis, or
#'   a vector of length 2 representing the upper and lower limits the x-axis
#'   should take, or a vector of length greater than two indicating which 
#'   x-values the derivative should be evaluated at
#' @param grid.length The default length of the \code{h.grid} or \code{x.grid} 
#'   if the length of either is not given.
#' @param derv The order of derivative for which to make the SiZer map.
#' @param degree The degree of the local weighted polynomial used to smooth the data.
#'   This must be greater than or equal to \code{derv}.
#'   
#' @details SiZer stands for the Significant Zero crossings of the derivative.  There are two 
#' dominate approaches in smoothing bivariate data: locally weighted regression or penalized splines.
#' Both approaches require the use of a 'bandwidth' parameter that controls how much smoothing
#' should be done.  Unfortunately there is no uniformly best bandwidth selection procedure.
#' SiZer (Chaudhuri and Marron, 1999) is a procedure that looks across a range of bandwidths
#' and classifies the p-th derivative of the smoother into one of three states: significantly
#' increasing (blue), possibly zero (purple), or significantly negative (red).
#' 
#' @return Returns list object of type SiZer which has the following components:
#'  \describe{
#'    \item{x.grid}{Vector of x-values at which the derivative was evaluated.}
#'    \item{h.grid}{Vector of bandwidth values for which a smoothing function was calculated.} 
#'    \item{slopes}{Matrix of what category a particular x-value and bandwidth falls into 
#'                  (Increasing=1, Possibly Zero=0, Decreasing=-1, Not Enough Data=2).}
#'  }
#' @references 
#'   Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures
#'     in curves. Journal of the American Statistical Association 94:807-823. 
#'     
#'   Hannig, J., and J. S. Marron. 2006. Advanced distribution theory for SiZer. 
#'     Journal of the American Statistical Association 101:484-499.
#'     
#'   Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect
#'     thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195.  
#' 
#' @author Derek Sonderegger
#' @seealso \code{\link{plot.SiZer}}, \code{\link{locally.weighted.polynomial}}
#' @export
#' @examples 
#' data('Arkansas')
#' x <- Arkansas$year
#' y <- Arkansas$sqrt.mayflies
#' 
#' plot(x,y)
#' 
#' # Calculate the SiZer map for the first derivative
#' SiZer.1 <- SiZer(x, y, h=c(.5,10), degree=1, derv=1, grid.length=21)
#' plot(SiZer.1)
#' 
#' # Calculate the SiZer map for the second derivative
#' SiZer.2 <- SiZer(x, y, h=c(.5,10), degree=2, derv=2, grid.length=21);
#' plot(SiZer.2)
#' 
#' # By setting the grid.length larger, we get a more detailed SiZer
#' # map but it takes longer to compute. 
#' #
#' # SiZer.3 <- SiZer(x, y, h=c(.5,10), grid.length=100, degree=1, derv=1)
#' # plot(SiZer.3)  
#'   
SiZer <- function(x, y, h=NA, x.grid=NA, degree=NA, derv=1, grid.length=41){
  
  # calculate x.grid and h.grid from what was passed in.
  x.grid <- x.grid.create(x.grid, x, y, grid.length);
  h.grid <- h.grid.create(h, x, y, grid.length);
  
  # set up degree if it wasn't passed in
  if( is.na(degree) ){
  	 degree = derv+1;
  }
  
	row <- 1;
	slopes <- matrix(nrow=length(h.grid), ncol=length(x.grid));
	for( h in h.grid ){
		print(h)
		model <- locally.weighted.polynomial(x, y, h=h, x.grid=x.grid, degree=degree);     
  		intervals <- 
  	      calc.CI.LocallyWeightedPolynomial(model, derv=derv);
  		slopes[row,] <- find.states(intervals);         
		row <- row + 1;
	}

  out <- NULL;
  out$x.grid <- x.grid;
  out$h.grid <- h.grid;
  out$slopes <- slopes;
  class(out) <- 'SiZer';
  return(out);
}		


#' Plot a SiZer map
#' Plot a \code{SiZer} object that was created using \code{SiZer()}
#' 
#' @param x An object created using \code{SiZer()}
#' @param ylab What the y-axis should be labled.
#' @param colorlist What colors should be used.  This is a vector that 
#'        corresponds to 'decreasing', 'possibley zero', 'increasing', 
#'        and 'insufficient data'.
#' @param \dots Any other parameters to be passed to the function \code{image}.
#' 
#' @details The white lines in the SiZer map give a graphical representation 
#'          of the bandwidth.  The horizontal distance between the lines is \eqn{2h}.
#'  
#' @references 
#'   Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures
#'     in curves. Journal of the American Statistical Association 94:807-823. 
#'     
#'   Hannig, J., and J. S. Marron. 2006. Advanced distribution theory for SiZer. 
#'     Journal of the American Statistical Association 101:484-499.
#'     
#'   Sonderegger, D.L., Wang, H., Clements, W.H., and Noon, B.R. 2009. Using SiZer to detect
#'     thresholds in ecological data. Frontiers in Ecology and the Environment 7:190-195.  
#' 
#' @author Derek Sonderegger
#' @seealso \code{\link{plot.SiZer}}, \code{\link{locally.weighted.polynomial}}
#' @export
#' @examples 
#' data('Arkansas')
#' x <- Arkansas$year
#' y <- Arkansas$sqrt.mayflies
#' 
#' plot(x,y)
#' 
#' # Calculate the SiZer map for the first derivative
#' SiZer.1 <- SiZer(x, y, h=c(.5,10), degree=1, derv=1, grid.length=21)
#' plot(SiZer.1)
#' 
#' # Calculate the SiZer map for the second derivative
#' SiZer.2 <- SiZer(x, y, h=c(.5,10), degree=2, derv=2, grid.length=21);
#' plot(SiZer.2)
#' 
#' # By setting the grid.length larger, we get a more detailed SiZer
#' # map but it takes longer to compute. 
#' #
#' # SiZer.3 <- SiZer(x, y, h=c(.5,10), grid.length=100, degree=1, derv=1)
#' # plot(SiZer.3)  
#'   
plot.SiZer <- function(x, ylab=expression(log[10](h)), 
			colorlist=c('red', 'purple', 'blue', 'grey'), ...){
	temp <- factor(x$slopes);
    final.colorlist <- NULL;
	if( is.element( '-1', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[1]);
	if( is.element( '0', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[2]);
	if( is.element( '1', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[3]);
	if( is.element( '2', levels(temp) ) )
	   final.colorlist <- c(final.colorlist, colorlist[4]);

	# Convert the slopes to a factor list to match up with final.colorlist
	temp <- matrix( as.integer(factor(x$slopes)), nrow=dim(x$slopes)[1] ) 

	graphics::image( x$x.grid, log(x$h.grid,10), t(temp), 
			col=final.colorlist, ylab=ylab, ...)
			
    # draw the bandwidth lines
    x.midpoint <- diff(range(x$x.grid))/2 + min(x$x.grid);
    graphics::lines( x.midpoint + x$h.grid, log(x$h.grid, 10), col='white' );
    graphics::lines( x.midpoint - x$h.grid, log(x$h.grid, 10), col='white' );                     
}	


h.grid.create <- function(h.grid, x, y, grid.length){
	foo <- grid.length;
	h.max <- diff( range(x) ) * 2;
	h.min <- max( diff(sort(x)) );
	if(all(is.na(h.grid))){
		# if we are passed an NA
		out <- 10^seq(log(h.min,10), log(h.max,10), length=grid.length);
	}else if( length(h.grid) == 1 ){
		# if we are passed a single integer
		out <- 10^seq(log(h.min,10), log(h.max,10), length=h.grid);
	}else if( length(h.grid) == 2 ){
		# passed two parameters: assume min and max for h.grid
		out <- 10 ^ seq(log(h.grid[1],10), 
                 		log(h.grid[2],10), length=grid.length);
	}else{
		# otherwise just return what we are sent
		out <- h.grid;
	} 
	return(out);
}


