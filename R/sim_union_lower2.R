#' sim_union_lower2.R
#'
#' Simulate the lower union bound R2.
#'
#' @param B Number of simulations.
#' @param muVec Mean vector of MVN (under the alternative).
#' @param sigMat Covariance matrix of MVN (under the alternative).
#' @param bounds A d*1 vector of bounds on the magnitudes of the test statistics, where
#' the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
#' If provided, will calculate power from crossing the bounds.
#'
#' @return Simulated lower bound quantity.
#'
#' @export
#' @examples
#' myVariance <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(myVariance) <- 1
#' myBounds <- set_GBJ_bounds(alpha = 0.01, d=5, sig_vec = myVariance[lower.tri(myVariance)])
#' sim_union_lower2(B=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_union_lower2 <- function(B, muVec, sigMat, bounds) {

  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=B, mean=muVec, sigma=sigMat)
  J <- length(bounds)

  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkUnionLower2 <- function(x, bounds) {
    # the union bound, either of first two larger
    top2idx <- kit::topn(abs(x), n=2, decreasing=TRUE)
    top2X <- sort(abs(x[top2idx]), decreasing = FALSE)
    top2Bounds <- bounds[(J-1):J]
    if (sum(as.numeric(top2X > top2Bounds)) > 0) {
      return(1)
    } else {
      return(0)
    }
  }

  # apply and average
  lowerCross <- apply(zSample, 1, checkUnionLower2, bounds = bounds)
  lowerBound <- mean(lowerCross)

  # return
  return(lowerBound = lowerBound)
}



