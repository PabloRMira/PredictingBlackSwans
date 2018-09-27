#' Compute the Average Variance (AVariance) of an estimator
#'
#' Computes the Average Variance (AVariance) of an estimator
#' in a Monte Carlo simulations study.
#' @param modelFit A matrix of fitted probabilities of dimensions (number of observations x
#' number of simulations).
#' @return A numeric scalar value.
#' @keywords PredictingBlackSwans
#' @export
AVar_fun <- function(modelFit) {
  nSim <- ncol(modelFit)
  n <- nrow(modelFit)
  expectedFit <- rowMeans(modelFit)
  expFitMat <- matrix(rep(expectedFit, nSim), n, nSim)
  sqDev <- (modelFit - expFitMat)^2
  indVar <- rowMeans(sqDev)
  aveVar <- mean(indVar)
  return(aveVar)
}
