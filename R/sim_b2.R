#' sim_b2.R
#'
#' Simulate the lower or upper bound R2/S2.
#'
#' @param lower Boolean, if true sim lower bound.
#' @param upper Boolean, if true sim upper bound.
#' @param B Number of simulations.
#' @param muVec Mean vector of MVN (under the alternative).
#' @param sigMat Covariance matrix of MVN (under the alternative).
#' @param bounds A d*1 vector of bounds on the magnitudes of the test statistics, where
#' the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
#' If provided, will calculate power from crossing the bounds.
#'
#' @return A list with the elements:
#' \item{lowerBound}{Lower bound on power.}
#' \item{upperBound}{Upper bound on power.}
#'
#' @export
#' @examples
#' myVariance <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(myVariance) <- 1
#' myBounds <- set_GBJ_bounds(alpha = 0.01, d=5, sig_vec = myVariance[lower.tri(myVariance)])
#' sim_b2(n=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_b2<- function(lower=TRUE, upper=FALSE, B, muVec, sigMat, bounds) {

  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=B, mean=muVec, sigma=sigMat)
  J <- length(bounds)

  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsLower2 <- function(x, bounds) {
    # This is calculating R2c, for lower bounds
    # In other words, the intersection of not crossing bounds
    top2idx <- kit::topn(abs(x), n=2, decreasing=TRUE)
    top2X <- sort(abs(x[top2idx]), decreasing = FALSE)
    top2Bounds <- bounds[(J-1):J]
    if (sum(as.numeric(top2X > top2Bounds)) == 0) {
      return(1)
    } else {
      return(0)
    }
  }

  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsUpper2 <- function(x, bounds) {
    # This is calculating S2c, for upper bounds
    # In other words, the intersection of not crossing bounds
    newBounds <- bounds
    newBounds[1:(J-2)] <- bounds[1]
    sortedX <- sort(abs(x), decreasing = FALSE)
    if (sum(as.numeric(sortedX > newBounds)) == 0) {
      return(1)
    } else {
      return(0)
    }
  }

  # upper or lower
  if (lower) {
    lowerNonCross <- apply(zSample, 1, checkBoundsLower2, bounds = bounds)
    lowerBound <- 1 - mean(lowerNonCross)
  } else {
    lowerBound <- NA
  }
  if (upper) {
    upperNonCross <- apply(zSample, 1, checkBoundsUpper2, bounds = bounds)
    upperBound <- 1 - mean(upperNonCross)
  } else {
    upperBound <- NA
  }

  # return
  return( list(lowerBound = lowerBound, upperBound = upperBound) )
}



