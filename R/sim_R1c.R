#' sim_R1c.R
#'
#' Simulate the upper bound R1c.
#'
#' @param n Number of simulations.
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
#' sim_R1c(n=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_R1c<- function(n, muVec, sigMat, bounds) {
  
  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=n, mean=muVec, sigma=sigMat)
  J <- length(bounds)
  
  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsR1c <- function(x, bounds) {
    newBounds <- bounds
    newBounds[1:(J-1)] <- bounds[1]
    sortedX <- sort(abs(x), decreasing = FALSE)
    if (sum(as.numeric(sortedX > newBounds)) == 0) {
      return(0)
    } else {
      return(1)
    }
  }
  
  crossVec <- apply(zSample, 1, checkBoundsR1c)
  boundsPower <- mean(crossVec)
  
  # return
  return( boundsPower = boundsPower )
}



