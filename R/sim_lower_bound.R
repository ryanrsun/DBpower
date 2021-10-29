#' sim_lower_bound.R
#'
#' Simulate the lower bounds R1 and R2.
#'
#' @param n Number of simulations.
#' @param dim 1 for R1 or 2 for R2
#' @param muVec Mean vector of MVN (under the alternative).
#' @param sigMat Covariance matrix of MVN (under the alternative).
#' @param bounds A d*1 vector of bounds on the magnitudes of the test statistics, where
#' the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
#' If provided, will calculate power from crossing the bounds.
#'
#' @return Lower bound on the probability.
#'
#' @export
#' @examples
#' myVariance <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(myVariance) <- 1
#' myBounds <- set_GBJ_bounds(alpha = 0.01, d=5, sig_vec = myVariance[lower.tri(myVariance)])
#' sim_lower_bound(n=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_lower_bound<- function(n, dim, muVec, sigMat, bounds) {
  
  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=n, mean=muVec, sigma=sigMat)
  J <- length(bounds)
  
  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsLower <- function(x, dim, bounds) {
    topIdx <- kit::topn(abs(x), n = dim, dereasing=TRUE)
    topX <- abs(x)[topIdx]
    topBounds <- bounds[(J - dim + 1):J]
    if (sum(as.numeric(topX > topBounds)) == 0) {
      return(0)
    } else {
      return(1)
    }
  }
  
  crossVec <- apply(zSample, 1, checkBoundsLower)
  crossBoundsProb <- mean(crossVec)
  
  # return
  return( 1 - crossBoundsProb )
}



