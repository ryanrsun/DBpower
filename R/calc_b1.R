#' calc_b1.R
#'
#' Calculate lower bound 1 - R1 or upper bound 1 - S1.
#'
#' @param lower Boolean, whether to calculate lower bound.
#' @param upper Boolean, whether to calculate upper bound.
#' @param muVec Mean vector of MVN (under the alternative).
#' @param sigMat Covariance matrix of MVN (under the alternative).
#' @param bounds A d*1 vector of bounds on the magnitudes of the test statistics, where
#' the first element is the bound for |Z|_(1) and the last element is the bound for |Z|_(d).
#' If provided, will calculate power from crossing the bounds.
#'
#' @return A list with the elements:
#' \item{allProbsLower}{J*1 vector of all component probabilities.}
#' \item{lowerProb}{Lower bound.}
#' \item{allProbsUpper}{J*1 vector of all component probabilities.}
#' \item{upperProb}{Upper bound.}
#'
#' @export
#' @examples
#' myVariance <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(myVariance) <- 1
#' myBounds <- set_GBJ_bounds(alpha = 0.01, d=5, sig_vec = myVariance[lower.tri(myVariance)])
#' calcb1(muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
calcb1 <- function(lower = TRUE, upper = FALSE, muVec, sigMat, bounds) {
  
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
  performIntegralLower1 <- function(j, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2) {
    transformMat <- createMj(j = j, size = length(muVec))
    newMu <- as.numeric(transformMat %*% muVec)
    newSig <- transformMat %*% sigMat %*% t(transformMat)
    intOut1 <- pmvnorm(lower = lBounds1, upper = uBounds1, mean = newMu, sigma = newSig) 
    intOut2 <- pmvnorm(lower = lBounds2, upper = uBounds2, mean = newMu, sigma = newSig) 
    return(intOut1[1] + intOut2[1])
  }
  # for upper
  performIntegralUpper1 <- function(j, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2) {
    # need to add to the transform mat
    transformMatTop <- createMj(j = j, size = length(muVec)) 
    addBottom <- diag(rep(1, length(muVec)))[-j, ]
    transformMat <- rbind(transformMatTop, addBottom)

    newMu <- as.numeric(transformMat %*% muVec)
    newSig <- transformMat %*% sigMat  %*% t(transformMat)
    intOut1 <- pmvnorm(lower = lBounds1, upper = uBounds1, mean = newMu, sigma = newSig) 
    intOut2 <- pmvnorm(lower = lBounds2, upper = uBounds2, mean = newMu, sigma = newSig) 
    return(intOut1[1] + intOut2[1])
  }
  
  # bounds
  lowerBounds1 <- rep(0, 2*J - 1)
  upperBounds1 <- c(bounds[J], rep(Inf, 2*J - 2))
  lowerBounds2 <- c(-bounds[J], rep(-Inf, 2*J - 2))
  upperBounds2 <- rep(0, 2*J - 1)
  
  # if calculating lower bound
  if (lower) {
    allProbsLower <- sapply(X = 1:J, FUN = performIntegralLower1, muVec = muVec, sigMat = sigMat, 
                       lBounds1 = lowerBounds1, uBounds1 = upperBounds1, lBounds2 = lowerBounds2,
                       uBounds2 = upperBounds2)
    lowerProb <- 1 - sum(allProbsLower[which(allProbsLower > 0)])
  } else {
    allProbsLower <- NA
    lowerProb <- NA
  }
  
  # if calculating upper bound
  if (upper) {
    # add to the bounds the -b1 \leq x \leq b1 part
    lowerExt <- rep(-bounds[1], J - 1)
    upperExt <- rep(bounds[1], J - 1)
    allProbsUpper <- sapply(X = 1:J, FUN = performIntegralUpper1, muVec = muVec, sigMat = sigMat, 
                            lBounds1 = c(lowerBounds1, lowerExt), uBounds1 = c(upperBounds1, upperExt), 
                            lBounds2 = c(lowerBounds2, lowerExt), uBounds2 = c(upperBounds2, upperExt))
    upperProb <- 1 - sum(allProbsUpper[which(allProbsUpper > 0)])
  } else {
    allProbsUpper <- NA
    upperProb <- NA
  }
  
  # return
  return(list(allProbsLower = allProbsLower, lowerProb = lowerProb,
              allProbsUpper = allProbsUpper, upperProb = upperProb))
}


