# this implements the bent-cable model
# and it finds the MLE of the threshold and
# bend width by exhaustively searching the
# parameter space.  This is rather time 
# intensive.  
#
# What we should probably do is create a wrapper
# function that first tries a search algorithm
# and if it fails then resorts to the exhaustive
# search.   


#' Fits a bent-cable model to the given data
#' Fits a bent-cable model to the given data by exhaustively searching
#' the 2-dimensional parameter space to find the maximum likelihood
#' estimators for \eqn{\alpha} and \eqn{\gamma}.
#' 
#' @param x The independent variable
#' @param y The dependent variable
#' @param grid.size How many \eqn{\alpha} and \eqn{gamma} values to examine.  
#'                  The total number of parameter combinations examined is 
#'                  \code{grid.size} squared.
#' @details Fit the model which is essentially a piecewise linear model with a 
#'          quadratic curve of length \eqn{2\gamma} connecting the two linear pieces.
#'          
#'          The reason for searching the space exhaustively is because the bent-cable 
#'          model often has a likelihood surface with a very flat ridge instead of 
#'          definite peak.  While the exhaustive search is slow, at least it is possible 
#'          to examine the contour plot of the likelihood surface.
#'          
#'  @return 	A list of 7 elements: \describe{
#'    \item{log.likelihood}{A matrix of log-likelihood values.}
#'    \item{SSE}{A matrix of sum-of-square-error values.}
#'    \item{alphas}{A vector of alpha values examined.}
#'    \item{gammas}{A vector of gamma values examined.}
#'    \item{alpha}{The MLE estimate of alpha.}
#'    \item{gamma}{The MLE estimate of gamma.}
#'    \item{model}{The \code{lm} fit after \eqn{alpha} and \eqn{gamma} are known. }
#'    }
#'    
#' @references 
#' Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression 
#' theory and applications. Journal of the American Statistical Association 
#' 101:542-553.
#' 
#' Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for 
#' identifying ecological thresholds. Ecology 84:2034-2041.   
#' 
#' @author Derek Sonderegger
#' @seealso \code{\link{piecewise.linear}}
#' @export
#' @examples 
#' data(Arkansas)
#' x <- Arkansas$year
#' y <- Arkansas$sqrt.mayflies
#' 
#' # For a more accurate estimate, increase grid.size
#' model <- bent.cable(x,y, grid.size=20)
#' plot(x,y)
#' x.grid <- seq(min(x), max(x), length=200)
#' lines(x.grid, predict(model, x.grid), col='red')
#' 
bent.cable <- function(x,y, grid.size=100){
  r <- range(x);
  alpha.min <- r[1];
  alpha.max <- r[2];
  gamma.max <- (r[2] - r[1])/2;
  gamma.min <- gamma.max / grid.size;

  A <- seq(alpha.min, alpha.max, length.out=grid.size);
  G <- seq(gamma.min, gamma.max, length.out=grid.size);

  L <- matrix(nrow=grid.size, ncol=grid.size);
  SSE <- matrix(nrow=grid.size, ncol=grid.size);

  max.llikelihood <- -Inf;
  hat.alpha <- 0;
  hat.gamma <- 0;

  count.a <- 1;
  for(a in A){
    count.g <- 1;
    for(g in G){
        I.middle <- 1* (abs(x-a) < g);
        I.left   <- 1* ( x >= a+g );
        q <- (x-a+g)^2 / (4*g) * I.middle +
             (x-a)             * I.left;
        fit <- stats::lm(y ~ x + q);
        Beta <- fit$coefficients;
        sigma <- sqrt( sum(fit$residuals^2)/length(x) );
        mu <- Beta[1] + Beta[2]*x + Beta[3]*q;

        llikelihood <- sum( log( stats::dnorm(y, mean=mu, sd=sigma) ) );
        L[count.a, count.g] <- llikelihood;
		SSE[count.a, count.g] <- sum(fit$residuals^2);
        if(llikelihood > max.llikelihood){
  	      max.llikelihood <- llikelihood;
          hat.alpha <- a;
          hat.gamma <- g;
        }
      count.g <- count.g + 1;	 
    }
    count.a <- count.a + 1;
  }
  out <- NULL;
  out$log.likelihood <- L;
  out$SSE <- SSE;
  out$alphas <- A;
  out$gammas <- G;
  out$alpha <- hat.alpha; 
  out$gamma <- hat.gamma; 

  I.middle <- 1* (abs(x-hat.alpha) < hat.gamma);
  I.left   <- 1* ( x >= hat.alpha + hat.gamma );
  q <- (x-hat.alpha+hat.gamma)^2 / (4*hat.gamma) * I.middle +
       (x-hat.alpha)             * I.left;
  fit <- stats::lm(y ~ x + q);
  
  out$model <- fit;
  class(out) <- 'bent_cable';
  return(out);
}





#' Return model predictions for fitted bent-cable model
#' 
#' @param object A bent-cable model
#' @param x The set x-values for which predictions are desired
#' @param \dots A placeholder that is currently ignored.  
#' @export
predict.bent_cable <- function(object, x, ...){
   alpha <- object$alpha;
   gamma <- object$gamma;
   beta <- object$model$coefficients;

  I.middle <- 1* (abs(x-alpha) < gamma);
  I.left   <- 1* ( x >= alpha + gamma );
  q <- (x - alpha + gamma)^2 / (4*gamma) * I.middle +
       (x - alpha) * I.left;

  yhat <- beta[1] + beta[2]*x + beta[3]*q;
  return(yhat);
}

#' Return the log-Likelihood value for a fitted bent-cable model.
#' 
#' @param object A bent-cable model
#' @param \dots Unused at this time.
#' @export
logLik.bent_cable <- function(object, ...){
	out <- stats::logLik(object$model)
	attr(out,'df') <- 6  # 3 betas, 1 gamma, 1 alpha, 1 sigma
	return(out)
}



