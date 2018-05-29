#' Creates a piecewise linear model
#' 
#' Fit a degree 1 spline with 1 knot point where the location of the 
#' knot point is unknown.
#' 
#' @param x Vector of data for the x-axis.
#' @param y Vector of data for the y-axis
#' @param middle A scalar in \eqn{[0,1]}. This represents the range that the change-point
#'   can occur in.  \eqn{0} means the change-point must occur at the middle of the range of 
#'   x-values.  \eqn{1} means that the change-point can occur anywhere along the range of the 
#'   x-values. 
#' @param CI Whether or not a bootstrap confidence interval should be calculated. Defaults to
#'   FALSE because the interval takes a non-trivial amount of time to calculate
#' @param bootstrap.samples The number of bootstrap samples to take when calculating the CI.
#' @param sig.level What significance level to use for the confidence intervals.
#' 
#' @details The bootstrap samples are taken by resampling the raw data points.  
#'   Sometimes a more appropriate bootstrap sample would be to calculate the residuals and 
#'   then add a randomly selected residual to each y-value.
#' 
#' @return  A list of 5 elements is returned: \describe{
#'   \item{change.point}{The estimate of \eqn{\alpha}.}
#'   \item{model}{The resulting \code{lm} object once \eqn{\alpha} is known. }
#'   \item{x}{The x-values used.}
#'   \item{y}{The y-values used.}
#'   \item{CI}{Whether or not the confidence interval was calculated.}
#'   \item{intervals}{If the CIs where calculated, this is a matrix of the
#'                    upper and lower intervals.}
#' }
#' 
#' @references 
#' Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression 
#' theory and applications. Journal of the American Statistical Association 
#' 101:542-553.
#' 
#' Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for 
#' identifying ecological thresholds. Ecology 84:2034-2041. 
#' 
#' @seealso The package \code{segmented} has a much more general implementation
#'   of this analysis and users should preferentially use that package.
#'   
#' @examples 
#' data(Arkansas)
#' x <- Arkansas$year
#' y <- Arkansas$sqrt.mayflies
#' 
#' model <- piecewise.linear(x,y, CI=FALSE)
#' plot(model)
#' print(model)
#' predict(model, 2001)
#' @export
piecewise.linear <- function(x, y, middle=1, CI=FALSE, bootstrap.samples=1000, 
                             sig.level=.05){
  
  
  alpha <- piecewise.linear.simple(x,y,middle);
  w <- x-alpha;
  w[w<0] <- 0;
  model <- stats::lm(y~x+w);
  out <- NULL;
  out$change.point <- alpha;
  out$model <- model;
  out$x <- seq(min(x), max(x), length=200);
  w <- out$x-alpha;
  w[w<0] <- 0;
  out$y <- stats::predict(out$model, data.frame(x=out$x, w=w) );
  out$CI <- CI;
  class(out) <- 'PiecewiseLinear';

  # if the user requests confidence intervals for the change point
  # then we'll bootstrap.  Otherwise just return what we have.
  if(CI == FALSE){
    return(out);
  }else{
    data <- data.frame(x=x, y=y);
    # define a function that will return the statistics we wish to bootstrap     
    my.cp <- function(data, index){
      x <- data[index, 1];
      y <- data[index, 2];
      cp <- piecewise.linear.simple(x,y);
      w <- x-cp;
      w[w<0] <- 0;
      model <- stats::lm(y~x+w);
      out <- c(cp, model$coefficients[2],  model$coefficients[3],
                   model$coefficients[2] + model$coefficients[3] );
      return(out);
    }
    boot.result <- boot::boot(data, my.cp, R=bootstrap.samples);
    out$intervals <- apply(boot.result$t, 2, stats::quantile, probs=c(sig.level/2, 1-sig.level/2));
    colnames(out$intervals) <- c('Change.Point', 'Initial.Slope', 
                          'Slope.Change', 'Second.Slope');
    out$CI <- t(out$CI);
    return(out);
  }
}



# This will perform a search for the best changepoint as defined by the 
# maximum likelihood surface.  It uses a search algorithm as opposed 
# to an exhaustive search.  
#
# the middle parameter can be used to narrow the search for a threshold
# by confining the search to the middle data points.  That is, middle=1 and
# all points along the x-axis (between the smallest and largest x values) 
# are possible thresholds.  Middle=.5 confines the search to the middle 50%
# of the range of x values.
piecewise.linear.simple <- function(x, y, middle=1){

  # This is the function that we wish to optimize (over alpha)
  piecewise.linear.likelihood <- function(alpha, x, y){
    N <- length(x);
    w <- (x-alpha);
    w[w<0] <- 0;
    fit <- stats::lm(y ~ x + w);
    Beta <- stats::coefficients(fit);
    Mu <- Beta[1] + Beta[2]*x + Beta[3]*w;    
    SSE <- sum(fit$residuals^2);
    sigma2 <- SSE/N;                    # MLE estimate of sigma^2
    likelihood <- sum( log( stats::dnorm(y, mean=Mu, sd=sqrt(sigma2)) ) );
    return(likelihood);
  }

  r <- range(x);
  offset <- r * (1-middle)/2;
  low <- min(x)  + offset;
  high <- max(x) - offset;
  temp <- stats::optimize(piecewise.linear.likelihood, c(low, high), x=x, y=y, maximum=TRUE);
  return(temp$maximum);
}
  



#' Prints out the model form for a Piecewise linear model
#' 
#' @param x A \code{PiecewiseLinear} object
#' @param \dots Unused at this time.
#' @export
print.PiecewiseLinear <- function(x, ...){
  print(paste('Threshold alpha:', x$change.point));
  print("");
  print("Model coefficients: Beta[0], Beta[1], Beta[2]");
  print(x$model$coefficients);
  if( x$CI == TRUE){
    print(x$intervals);
  }
}

#' Plots a piecewise linear model
#' 
#' @param x A \code{PiecewiseLinear} object
#' @param xlab The label for the x-axis
#' @param ylab The label for the y-axis
#' @param \dots Any further options to be passed to the \code{plot} function
#' @export
plot.PiecewiseLinear <- function(x, xlab='X', ylab='Y', ...){
	graphics::plot(x$model$model$x, x$model$model$y, xlab=xlab, ylab=ylab, ...);
	graphics::lines(x$x, x$y, col='red');
}

#' Calculates predicted values from a piecewise linear object
#' 
#' @param object A \code{PiecewiseLinear} object
#' @param x A vector of x-values in which to calculate the y.
#' @param \dots Unused at this time.
#' @export
predict.PiecewiseLinear <- function( object, x, ...){
	alpha <- object$change.point;
	beta <- object$model$coefficients;
	w <- x - alpha;
	w[w<0] <- 0;
	yhat <- beta[1] + beta[2]*x + beta[3]*w;
	return(yhat);
}

#' Calculates the log-Likelihood value
#' 
#' @param object A \code{PiecewiseLinear} object
#' @param \dots Unused at this time.
#' @export
logLik.PiecewiseLinear <- function(object, ...){
	out <- stats::logLik(object$model)
	attr(out,'df') <- 5  # 3 betas, 1 alpha, 1 sigma
	return(out)
}
