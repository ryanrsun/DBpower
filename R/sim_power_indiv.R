#' sim_stats_indiv.R
#'
#' Simulate distribution of test statistics from individual data for
#' multiple explanatory factor setting.
#'
#' @param B Number of simulations.
#' @param sigSq Variance of outcome.
#' @param xMat Design matrix of non-genetic covariates, n*p.
#' @param gMat Matrix of genotypes, n*J.
#' @param alphaVec p*1 vector of regression coefficients for xMat.
#' @param betaVec J*1 vector of regression coefficients for gMat.
#' @param checkpoint Report every 50 simulations.
#'
#' @return A list with the elements:
#' \item{tMat}{B*J matrix of T statistics.}
#' \item{zMat}{B*J matrix of Z statistics.}
#' \item{zVecGBJ}{Vector of Z statistics should match first row of zMat.}
#' \item{iMat}{Innovated statistics using correlation matrix under the null.}
#' @export
#' @examples
#' xMat <- cbind(1, rnorm(n = 1000), rbinom(n = 1000, size=1, prob=0.5))
#' gMat <- matrix(data = rbinom(n=10000, size=2, prob=0.3), nrow=1000)
#' alphaVec <- c(1, 1, 1)
#' betaVec <- rep(0, 10)
#' sim_stats_mef(B=10000, sigSq = 1, xMat = xMat, gMat = gMat, alphaVec = alphaVec, betaVec = betaVec)
sim_stats_mef <- function(B, sigSq, xMat, gMat, alphaVec, betaVec, checkpoint = FALSE) {

  # mean vector
  etaVec <- xMat %*% alphaVec + gMat %*% betaVec
  n <- length(etaVec)
  p <- ncol(xMat)
  q <- ncol(gMat)
  # quantities for calculating statistics and correlation matrix
  pMat <- xMat %*% solve(t(xMat) %*% xMat) %*% t(xMat)
  mafVec <- apply(gMat, 2, mean) / 2
  varG <- 2 * mafVec * (1 - mafVec)
  IminusP <- diag(rep(1, n)) - pMat
  GGterm <- t(gMat) %*% IminusP %*% gMat

  # record
  tMat <- matrix(data=NA, nrow=B, ncol=q)
  zMat <- matrix(data=NA, nrow=B, ncol=q)
  for (sim_it in 1:B) {
    # generate outcome
    yVec <- etaVec + rnorm(n=n, mean = 0, sd = sqrt(sigSq))
    # calculate statistics
    residVec <- yVec - pMat %*% yVec
    sigSqHat <- sum(residVec^2) / (n - p)
    tVec <- t(gMat) %*% residVec / sqrt(n)
    zVec <- sqrt(n) * tVec / sqrt(diag(GGterm) * sigSqHat)

    # do GBJ
    if (sim_it == 1) {
      nullMod <- glm(yVec ~ xMat - 1, family=gaussian())
      #packageOutput <- calc_score_stats(null_model = nullMod, factor_matrix = gMat, link_function = "linear")
      packageOutput <- calc_score_stats(null_model = nullMod, factor_matrix = gMat, link_function = "linear", P_mat = sigSqHat * IminusP)
      zVecGBJ <- packageOutput$test_stats
    }

    # record
    tMat[sim_it, ] <- tVec
    zMat[sim_it, ] <- zVec

    # checkpoint
    if (checkpoint) {
      if (sim_it%%500 == 0) {cat(sim_it)}
    }
  }

  # innovate with estimated correlation matrix under the null
  denomVec <- 1 / sqrt(diag(GGterm))
  corMat <- diag(denomVec) %*% GGterm %*% diag(denomVec)
  decomp <- eigen(corMat)
  iMat <- zMat %*% decomp$vectors %*% diag(1 / sqrt(decomp$values))

  return(list(tMat = tMat, zMat = zMat, zVecGBJ = zVecGBJ, iMat = iMat))
}



#' Simulate power of detection boundary tests starting from individual data for
#' multiple outcomes setting.
#'
#' @param B Number of simulations.
#' @param sigSqVec Variance of outcomes.
#' @param xMat Design matrix of non-genetic covariates, n*p.
#' @param gVec n*1 vector of genotypes.
#' @param alphaMat p*K vector of regression coefficients for xMat.
#' @param gammaVec K*1 vector of regression coefficients for each outcome.
#' @param checkpoint Report every 50 simulations.
#'
#' @return A list with the elements:
#' \item{uMat}{Matrix of U statistics.}
#' \item{vMat}{Matrix of V statistics.}
#' \item{vVecGBJ}{Vector of V statistics should match first row of vMat}
#' \item{iMat}{Innovated statistics using correlation matrix under the null.}
#' @export
#' @examples
#' covY <- matrix(data=0.3, nrow=10, ncol=10); diag(covY) <- 1
#' xMat <- cbind(1, rnorm(n = 1000), rbinom(n = 1000, size=1, prob=0.5))
#' gVec <- rbinom(n= 1000, size = 2, prob=0.3)
#' alphaMat <-matrix(data = 1, nrow=3, ncol=10)
#' gammaVec <- rep(0, 10)
#' sim_stats_mo(B=10000, covY = covY, xMat = xMat, gVec = gVec, alphaMat = alphaMat, gammaVec = gammaVec)
sim_stats_mo <- function(B, covY, xMat, gVec, alphaMat, gammaVec, checkpoint = FALSE) {

  # mean matrix n*K
  etaMat <- xMat %*% alphaMat + gVec %*% t(gammaVec)
  n <- nrow(etaMat)
  K <-length(gammaVec)
  p <- ncol(xMat)

  # quantities for calculating statistics and correlation matrix
  pMat <- xMat %*% solve(t(xMat) %*% xMat) %*% t(xMat)
  IminusP <- diag(rep(1, n)) - pMat
  maf <- mean(gVec) / 2
  varG <- t(gVec) %*% IminusP %*% gVec
  varYvec <- diag(covY)
  corY <- sweep(covY, MARGIN=2, STATS=sqrt(varYvec), FUN="/") %>%
    sweep(., MARGIN=1, STATS=sqrt(varYvec), FUN="/")

  # record
  uMat <- matrix(data=NA, nrow=B, ncol=K)
  vMat <- matrix(data=NA, nrow=B, ncol=K)
  for (sim_it in 1:B) {
    # generate outcome
    yMat <- etaMat + mvtnorm::rmvnorm(n=n, mean=rep(0, K), sigma = covY)
    # calculate statistics
    residMat <- yMat - pMat %*% yMat
    sigSqHat <- apply(residMat^2, 2, sum) / (n - p)
    uVec <- t(gVec) %*% residMat / sqrt(n)
    vVec <- sqrt(n) * uVec / sqrt(rep(as.numeric(varG), K) * sigSqHat)

    # do package
    if (sim_it == 1) {
      zVecPackage <- rep(0, K)
      for (k in 1:K) {
        nullMod <- glm(yMat[, k] ~ xMat - 1, family=gaussian())
        #packageOutput <- calc_score_stats(null_model = nullMod, factor_matrix = gMat, link_function = "linear")
        packageOutput <- calc_score_stats(null_model = nullMod, factor_matrix = matrix(data=rep(gVec, 2), ncol=2), link_function = "linear", P_mat = sigSqHat[k] * IminusP)
        zVecPackage[k] <- packageOutput$test_stats[1]
      }
    }

    # record
    uMat[sim_it, ] <- uVec
    vMat[sim_it, ] <- vVec

    # checkpoint
    if (checkpoint) {
      if (sim_it%%500 == 0) {cat(sim_it)}
    }
  }

  # innovate with estimated correlation matrix under the null
  decomp <- eigen(corY)
  iMat <- vMat %*% decomp$vectors %*% diag(1 / sqrt(decomp$values))

  return(list(uMat = uMat, vMat = vMat, zVecPackage = zVecPackage, iMat = iMat))
}





