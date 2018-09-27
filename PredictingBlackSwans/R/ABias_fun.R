#' Compute the Squared Average Bias (Squared ABias)
#'
#' Computes the squared average bias (Squared ABias) of an estimator
#' in a Monte Carlo simulation study.
#' @param modelFit A matrix of fitted probabilities of dimensions (number of observations x
#' number of simulations).
#' @param truth The true probabilities of the data generating process at \eqn{x_i}, i.e.,
#' \eqn{\Lambda (x_i^T \beta)}.
#' @return A numeric scalar value.
#' @keywords PredictingBlackSwans
#' @export
ABias_fun <- function(modelFit, truth) {
  expectedFit <- rowMeans(modelFit)
  sqBias <- (expectedFit - truth)^2
  aveSqBias <- mean(sqBias)
  return(aveSqBias)
}
