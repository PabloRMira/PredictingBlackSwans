#' Compute the Average Mean Squared Predicion Error (AMSPE) of an estimator
#'
#' Computes the Average Mean Squared Prediction Error (AMSPE) of an estimator
#' in a Monte Carlo simulation study.
#' @param modelFit A matrix of fitted probabilities of dimensions (number of observations x
#' number of simulations).
#' @param truth The true probabilities of the data generating process at \eqn{x_i}, i.e.,
#' \eqn{\Lambda (x_i^T \beta)}.
#' @return A numeric scalar value.
#' @keywords PredictingBlackSwans
#' @export
AMSPE_fun <- function(modelFit, truth) {
  n <- length(truth)
  nSim <- ncol(modelFit)
  truthMat <- matrix(rep(truth, nSim), n, nSim)
  errorMat <- modelFit - truthMat
  sqErrorMat <- errorMat^2
  meanSqError <- rowMeans(sqErrorMat)
  aMeanSqError <- mean(meanSqError)
  return(aMeanSqError)
}
