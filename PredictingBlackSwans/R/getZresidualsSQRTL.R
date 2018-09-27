#' Compute the residuals of the square-root Lasso regressions
#'
#' @param i Column index for the response.
#' @param x Matrix of predictors.
#' @param lambda Either 'BCW' for the simulations method
#' proposed by Belloni, Chernozhukov and Wang (2011) or 'VdG' for
#' the proposed method in van de Geer (2014).
#' @return Vector of residuals.
#' @export
getZresidualsSQRTL <- function(j,
                               x,
                               lambda) {
  # Size of data
  nObs <- nrow(x)
  nVars <- ncol(x)
  # Response
  xj <- x[, j]
  # Predictors
  xminusj <- x[, -j]
  if (lambda == "BCW") {
    # Numbre of simulations
    simNum <- 100
    # Recommended parameters by Belloni, Chernozhukov & Wang (2011)
    c <- 1.1
    alpha <- 0.05
    # Lambda parameter
    epsilonMat <- matrix(rnorm(nObs*simNum), nObs, simNum)
    # Normalize by the norms
    epsilonNorms <- sqrt(colSums(epsilonMat^2))
    epsilonMat <- scale(epsilonMat, center=F, scale=epsilonNorms)
    # Absolute gradient
    absGradient <- abs(t(xminusj) %*% epsilonMat)
    # Maximum absolute gradients
    maxGradient <- apply(absGradient, 2, max)
    # Find lambda
    bestlambda = as.vector(sqrt(nObs) * c * quantile(maxGradient, 1 - alpha))
  } else if (lambda == "VdG") {
    # Van de Geer's recommendation: lambda = sqrt(log(p) / n)
    # But in the software of Belloni, Chernozhukov and Wang, it corresponds to
    # lambda_{VdG} = lambda_{BCW} / n and therefore -> bestlambda = sqrt(log(p) * n)
    # because we input the BCW choice for the squre-root-Lasso software
    bestlambda <- sqrt(log(nVars - 1) * nObs)
  } else {
    stop(paste("Error: The selected lambda", lambda, "is not available."))
  }
  sqrtFit <- PredictingBlackSwans::sqrt_lasso(y=xj, X=xminusj, lambda=bestlambda)
  return(sqrtFit$Z)
}
