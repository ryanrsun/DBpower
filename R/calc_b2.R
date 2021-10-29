#' calc_b2.R
#'
#' Calculate lower bound 1 - R2 or upper bound 1 - S2.
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
#' calcb2(muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
calcb2 <- function(lower = TRUE, upper = FALSE, muVec, sigMat, bounds) {
  
  J <- length(muVec)
  muVec <- as.numeric(muVec)
  
  # we need the distribution of Mjk*Z, make the Mjk matrix here
  createMjk <- function(j, k, size) {
    topHalf <- diag(rep(-1, size))[-c(j, k), ]
    topHalf[, k] <- 1
    bottomHalf <- abs(topHalf)
    
    # put it together
    Mjk <- rbind(matrix(data=0, nrow=3, ncol=size), topHalf, bottomHalf)
    Mjk[1, j] <- 1
    Mjk[2, k] <- 1
    Mjk[3, j] <- 1
    Mjk[3, k] <- -1
    return(Mjk)
  }
  
  # apply this function to perform the MVN integral each time for lower bound
  performIntegralLower2 <- function(x, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2,
                                    lBounds3, uBounds3, lBounds4, uBounds4) {
    j <- x[1]
    k <- x[2]
    # two transform matrices, use the Opp one when Zj and Zk have different signs
    transformMat <- createMjk(j = j, k = k, size = length(muVec))
    transformMatOpp <- transformMat
    transformMatOpp[3, k] <- 1
    newMu <- as.numeric(transformMat %*% muVec)
    newSig <- transformMat %*% sigMat %*% t(transformMat)
    newMuOpp <- as.numeric(transformMatOpp %*% muVec)
    newSigOpp <- transformMatOpp %*% sigMat %*% t(transformMatOpp)
    # do the integrals
    intOut1 <- pmvnorm(lower = lBounds1, upper = uBounds1, mean = newMu, sigma = newSig) 
    intOut2 <- pmvnorm(lower = lBounds2, upper = uBounds2, mean = newMuOpp, sigma = newSigOpp) 
    intOut3 <- pmvnorm(lower = lBounds3, upper = uBounds3, mean = newMuOpp, sigma = newSigOpp) 
    intOut4 <- pmvnorm(lower = lBounds4, upper = uBounds4, mean = newMu, sigma = newSig) 
    
    # return
    return(intOut1[1] + intOut2[1] + intOut3[1] + intOut4[1])
  }

  # apply this function to perform the MVN integral each time for upper bound
  performIntegralUpper2 <- function(x, muVec, sigMat, lBounds1, uBounds1, lBounds2, uBounds2,
                                    lBounds3, uBounds3, lBounds4, uBounds4) {
    j <- x[1]
    k <- x[2]
    # need to add to the transform mat
    transformMatTop <- createMjk(j = j, k = k, size = length(muVec))
    transformMatOppTop <- transformMatTop
    transformMatOppTop[3, k] <- 1
    bottomPart <- diag(rep(1, length(muVec)))[-c(j, k), ]
    transformMat <- rbind(transformMatTop, bottomPart)
    transformMatOpp <- rbind(transformMatOppTop, bottomPart)
    # transform
    newMu <- as.numeric(transformMat %*% muVec)
    newSig <- transformMat %*% sigMat %*% t(transformMat)
    newMuOpp <- as.numeric(transformMatOpp %*% muVec)
    newSigOpp <- transformMatOpp %*% sigMat %*% t(transformMatOpp)
    
    # do integral with additional bounds
    lowerExt <- rep(-bounds[1], J - 1)
    upperExt <- rep(bounds[1], J - 1)
    intOut1 <- pmvnorm(lower = c(lBounds1), upper = c(uBounds1), mean = newMu, sigma = newSig)
    intOut2 <- pmvnorm(lower = c(lBounds2), upper = c(uBounds2), mean = newMuOpp, sigma = newSigOpp)
    intOut3 <- pmvnorm(lower = c(lBounds3), upper = c(uBounds3), mean = newMuOpp, sigma = newSigOpp)
    intOut4 <- pmvnorm(lower = c(lBounds4), upper = c(uBounds4), mean = newMu, sigma = newSig)
    return(intOut1[1] + intOut2[1] + intOut3[1] + intOut4[1])
  }
  
  # make bounds
  lowerBounds1 <- rep(0, 2*J - 1)
  upperBounds1 <- c(bounds[J], bounds[J - 1], rep(Inf, 2*J - 3))
  lowerBounds2 <- c(0, -bounds[J - 1], 0, rep(-Inf, 2*J - 4))
  upperBounds2 <- c(bounds[J], 0, Inf, rep(0, 2*J - 4))
  lowerBounds3 <- c(-bounds[J], 0, -Inf, rep(0, 2*J - 4))
  upperBounds3 <- c(0, bounds[J - 1], 0, rep(Inf, 2*J - 4))
  lowerBounds4 <- c(-bounds[J], -bounds[J - 1], rep(-Inf, 2*J - 3))
  upperBounds4 <- rep(0, 2*J - 1)
  
  # if calculating lower bound
  if (lower) {
    allCombs <- expand.grid(1:J, 1:J) %>%
      as.data.frame(.) %>%
      filter(Var1 != Var2) %>%
      as.matrix(.)
    allProbsLower <- apply(allCombs, 1, performIntegralLower2, muVec = muVec, sigMat = sigMat, 
                            lBounds1 = lowerBounds1, uBounds1 = upperBounds1, lBounds2 = lowerBounds2,
                            uBounds2 = upperBounds2, lBounds3 = lowerBounds3, uBounds3 = upperBounds3,
                            lBounds4 = lowerBounds4, uBounds4 = upperBounds4)
    lowerProb <- 1 - sum(allProbsLower[which(allProbsLower > 0)])
  } else {
    allProbsLower <- NA
    lowerProb <- NA
  }
  
  # if calculating upper bound
  if (upper) {
    # add to the bounds the -b1 \leq x \leq b1 part
    lowerExt <- rep(-bounds[1], J - 2)
    upperExt <- rep(bounds[1], J - 2)
    allProbsUpper <- apply(allCombs, 1, FUN = performIntegralUpper2, muVec = muVec, sigMat = sigMat, 
                            lBounds1 = c(lowerBounds1, lowerExt), uBounds1 = c(upperBounds1, upperExt), 
                            lBounds2 = c(lowerBounds2, lowerExt), uBounds2 = c(upperBounds2, upperExt),
                            lBounds3 = c(lowerBounds3, lowerExt), uBounds3 = c(upperBounds3, upperExt),
                            lBounds4 = c(lowerBounds4, lowerExt), uBounds4 = c(upperBounds4, upperExt))
    upperProb <- 1 - sum(allProbsUpper[which(allProbsUpper > 0)])
  } else {
    allProbsUpper <- NA
    upperProb <- NA
  }
  
  # return
  return(list(allProbsLower = allProbsLower, lowerProb = lowerProb,
              allProbsUpper = allProbsUpper, upperProb = upperProb))
}


