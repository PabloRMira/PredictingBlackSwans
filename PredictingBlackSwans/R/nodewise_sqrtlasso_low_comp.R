#' Compute the nodewise Lasso using the square-root Lasso for a low dimensional
#' component
#'
#' Compute the nodewise Lasso using the square-root Lasso to
#' calculate the matrix of residuals Z stemming from the nodewise
#' regressions for a low dimensional component.
#' @param x Predictor matrix.
#' @param lowSel Indixes of the low-dimensional selection of variables to desparsify.
#' @param parallel Should parallel computing be used?
#' @param ncores Number of cores for parallel computing.
#' @param lambda Either 'BCW' for the simulations method
#' proposed by Belloni, Chernozhukov and Wang (2011) or 'VdG' for
#' the proposed method in van de Geer (2014).
#' @return The matrix of residuals of the nodewise square-root Lasso regressions, i.e.
#' \eqn{Z = (Z_1, ..., Z_p)} with \eqn{Z_j \in R^n}.
#' @export
nodewise_sqrtlasso_low_comp <- function(x,
                                        lowSel   = 1:ncol(x),
                                        parallel = TRUE,
                                        ncores   = 2L,
                                        lambda   = "BCW") {
  # Get number of regressors
  p <- ncol(x)
  # Compute the matrix of residuals of the nodewise regressions with
  # the square root Lasso
  if (parallel) {
    Z <- parallel::mcmapply(PredictingBlackSwans::getZresidualsSQRTL, j = lowSel, x = list(x = x),
                            lambda = lambda, mc.cores = ncores)
  } else {
    Z <- mapply(PredictingBlackSwans::getZresidualsSQRTL, j = lowSel, x = list(x = x),
                lambda = lambda)
  }
  return(Z)
}
