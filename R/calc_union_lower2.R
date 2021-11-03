#' calc_union_lower2.R
#'
#' Calculate lower union bound.
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
#' calc_union_lower2(muVec = c(1, 0, 0, 0, 0), sigMat = myVariance, bounds=myBounds)
#'
calc_union_lower2 <- function(muVec, sigMat, bounds) {

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
  performIntegralUnion2 <- function(x, muVec, sigMat, bounds1, bounds2, bounds3, bounds4) {
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
    intOut11 <- pmvnorm(lower = bounds1[1, ], upper = bounds1[2, ], mean = newMu, sigma = newSig)[1]
    intOut12 <- pmvnorm(lower = bounds1[3, ], upper = bounds1[4, ], mean = newMu, sigma = newSig)[1]
    intOut21 <- pmvnorm(lower = bounds2[1, ], upper = bounds2[2, ], mean = newMuOpp, sigma = newSigOpp)[1]
    intOut22 <- pmvnorm(lower = bounds2[3, ], upper = bounds2[4, ], mean = newMuOpp, sigma = newSigOpp)[1]
    intOut31 <- pmvnorm(lower = bounds3[1, ], upper = bounds3[2, ], mean = newMuOpp, sigma = newSigOpp)[1]
    intOut32 <- pmvnorm(lower = bounds3[3, ], upper = bounds3[4, ], mean = newMuOpp, sigma = newSigOpp)[1]
    intOut41 <- pmvnorm(lower = bounds4[1, ], upper = bounds4[2, ], mean = newMu, sigma = newSig)[1]
    intOut42 <- pmvnorm(lower = bounds4[3, ], upper = bounds4[4, ], mean = newMu, sigma = newSig)[1]

    # return
    formerTerms <- intOut11 + intOut21 + intOut31 + intOut41
    latterTerms <- intOut12 + intOut22 + intOut32 + intOut42
    return(c(formerTerms, latterTerms, formerTerms +latterTerms))
  }

  # make bounds
  lowerBounds1A <- c(bounds[J], rep(0, 2*J - 2))
  upperBounds1A <- rep(Inf, 2*J - 1)
  lowerBounds1B <- c(0, bounds[J - 1], rep(0, 2*J - 3))
  upperBounds1B <- c(bounds[J], rep(Inf, 2*J - 2))
  bounds1 <- rbind(lowerBounds1A, upperBounds1A, lowerBounds1B, upperBounds1B)

  lowerBounds2A <- c(bounds[J], -Inf, 0, rep(-Inf, 2*J - 4))
  upperBounds2A <- c(Inf, 0, Inf, rep(0, 2*J - 4))
  lowerBounds2B <- c(0, -Inf, 0, rep(-Inf, 2*J - 4))
  upperBounds2B <- c(bounds[J], -bounds[J - 1], Inf, rep(0, 2*J - 4))
  bounds2 <- rbind(lowerBounds2A, upperBounds2A, lowerBounds2B, upperBounds2B)

  lowerBounds3A <- c(-Inf, 0, -Inf, rep(0, 2*J - 4))
  upperBounds3A <- c(-bounds[J], Inf, 0, rep(Inf, 2*J - 4))
  lowerBounds3B <- c(-bounds[J], bounds[J - 1], -Inf, rep(0, 2*J - 4))
  upperBounds3B <- c(0, Inf, 0, rep(Inf, 2*J - 4))
  bounds3 <- rbind(lowerBounds3A, upperBounds3A, lowerBounds3B, upperBounds3B)

  lowerBounds4A <- rep(-Inf, 2*J - 1)
  upperBounds4A <- c(-bounds[J], rep(0, 2*J - 2))
  lowerBounds4B <- c(-bounds[J], rep(-Inf, 2*J - 2))
  upperBounds4B <- c(0, -bounds[J - 1], rep(0, 2*J - 3))
  bounds4 <- rbind(lowerBounds4A, upperBounds4A, lowerBounds4B, upperBounds4B)

  # do the apply
  allCombs <- expand.grid(1:J, 1:J) %>%
    as.data.frame(.) %>%
    filter(Var1 != Var2) %>%
    as.matrix(.)
  allProbsLower <- apply(allCombs, 1, performIntegralUnion2, muVec = muVec, sigMat = sigMat,
                         bounds1 = bounds1, bounds2 = bounds2, bounds3 = bounds3, bounds4 = bounds4)
  probsDF <- as.data.frame(allCombs) %>% mutate(formerProbs = allProbsLower[1,],
                                                latterProbs = allProbsLower[2,],
                                                totalProbs = allProbsLower[3,])
  lowerProb <- sum(probsDF$totalProbs[which(probsDF$totalProbs > 0)])

  # return
  return(list(probsDF = probsDF, lowerProb = lowerProb))
}


