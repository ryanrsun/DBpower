#' calc_union_lower1.R
#'
#' Calculate lower bound R1. Only difference from calc_b1 lower option is change in first elements of integral bounds.
#'
#' @param muVec Mean vector of MVN (under the alternative).
#' @param sigMat Covariance matrix of MVN (under the alternative).
#' @param bounds A d*1 vector of bounds on the magnitudes of the test statistics, where
#' the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
#' If provided, will calculate power from crossing the bounds.
#'
#' @return A list with the elements:
#' \item{allProbsLower}{J*1 vector of all component probabilities.}
#' \item{lowerProb}{Lower bound.}
#'
#' @export
#' @examples
#' myVariance <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(myVariance) <- 1
#' myBounds <- set_GBJ_bounds(alpha = 0.01, d=5, sig_vec = myVariance[lower.tri(myVariance)])
#' calc_union_lower1(muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
calc_union_lower1 <- function(muVec, sigMat, bounds) {

  J <- length(muVec)

  # we need the distribution of Mj*Z, make the Mj matrix here
  createMj <- function(j, size) {
    topHalf <- diag(rep(-1, size))[-j, ]
    topHalf[, j] <- 1
    bottomHalf <- abs(topHalf)

    # put it together
    Mj <- rbind(0, topHalf, bottomHalf)
    Mj[1, j] <- 1
    return(Mj)
  }

  # apply this function to perform the MVN integral each time
  performIntegralUnion1 <- function(j, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2) {
    transformMat <- createMj(j = j, size = length(muVec))
    newMu <- as.numeric(transformMat %*% muVec)
    newSig <- transformMat %*% sigMat %*% t(transformMat)
    intOut1 <- pmvnorm(lower = lBounds1, upper = uBounds1, mean = newMu, sigma = newSig)
    intOut2 <- pmvnorm(lower = lBounds2, upper = uBounds2, mean = newMu, sigma = newSig)
    return(intOut1[1] + intOut2[1])
  }

  # bounds - only difference from calc_b1 lower option is first element of bounds
  lowerBounds1 <- c(bounds[J], rep(0, 2*J - 2))
  upperBounds1 <- c(Inf, rep(Inf, 2*J - 2))
  lowerBounds2 <- c(-Inf, rep(-Inf, 2*J - 2))
  upperBounds2 <- c(-bounds[J], rep(0, 2*J - 2))

  # if calculating lower bound
  allProbsLower <- sapply(X = 1:J, FUN = performIntegralUnion1, muVec = muVec, sigMat = sigMat,
                          lBounds1 = lowerBounds1, uBounds1 = upperBounds1, lBounds2 = lowerBounds2,
                          uBounds2 = upperBounds2)
  lowerProb <- sum(allProbsLower[which(allProbsLower > 0)])

  # return
  return(list(allProbsLower = allProbsLower, lowerProb = lowerProb))
}


