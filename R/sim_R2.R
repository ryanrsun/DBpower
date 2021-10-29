#' sim_R2.R
#'
#' Simulate the lower bound R2.
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
#' sim_R2(n=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_R2<- function(n, muVec, sigMat, bounds) {
  
  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=n, mean=muVec, sigma=sigMat)
  J <- length(bounds)
  
  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsR2 <- function(x, bounds) {
    top2idx <- kit::topn(abs(x), n=2, dereasing=TRUE)
    top2X <- abs(x)[top2idx]
    top2Bounds <- bounds[(J-1):J]
    if (sum(as.numeric(top2X > top2Bounds)) == 0) {
      return(0)
    } else {
      return(1)
    }
  }
  
  crossVec <- apply(zSample, 1, checkBoundsR2)
  boundsPower <- mean(crossVec)
  
  # return
  return( boundsPower = boundsPower )
}



