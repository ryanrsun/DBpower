#' sim_R1.R
#'
#' Simulate the lower bound R1 or upper bound S1.
#'
#' @param lower Boolean, if true then simulate lower bound.
#' @param upper Boolean, if true then simulate upper bound.
#' @param n Number of simulations.
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
#' sim_b1(n=5000, muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
sim_b1<- function(n, muVec, sigMat, bounds) {

  # simulate MVN
  zSample <- mvtnorm::rmvnorm(n=n, mean=muVec, sigma=sigMat)

  # function to check if bounds are violated
  # bounds should be smallest to largest
  checkBoundsLower1 <- function(x, bounds) {
    # This is calculating R1c, for lower bounds
    # The probability of not crossing the bounds, 1 minus this is the power
    maxX <- max(abs(x))
    maxBound <- max(bounds)
    return(as.numeric(maxX < maxBound))
  }

  # function to check if bounds are violated
  checkBoundsUpper1 <- function(x, bounds) {
    # This is calculating S1c, for upper bounds
    # The probability of not crossing the bounds, 1 minus this is the power
    newBounds <- bounds
    newBounds[1:(J-1)] <- bounds[1]
    sortedX <- sort(abs(x), decreasing = FALSE)
    if (sum(as.numeric(sortedX > newBounds)) == 0) {
      return(1)
    } else {
      return(0)
    }
  }

  # upper or lower
  if (lower) {
    lowerNonCross <- apply(zSample, 1, checkBoundsLower1)
    lowerBound <- 1 - mean(lowerNonCross)
  } else {
    lowerBound <- NA
  }
  if (upper) {
    upperNonCross <- apply(zSample, 1, checkBoundsUpper1)
    upperBound <- 1 - mean(upperNonCross)
  } else {
    upperBound <- NA
  }

  # return
  return( list(lowerBound = lowerBound, upperBound = upperBound) )
}



