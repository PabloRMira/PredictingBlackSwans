#' Compute the residuals of the nodewise Lasso regressions
#'
#' @param i Column index of the response.
#' @param x Matrix of predictors.
#' @param lambda Regularization parameter.
#' @return Vector of residuals.
#' @export
getZresiduals <- function(i,
                          x,
                          lambda) {
  glmnetfit  <- glmnet::glmnet(x[, -i],x[, i], lambda=lambda)
  prediction <- glmnet::predict.glmnet(glmnetfit, x[, -i], s=lambda)
  return(x[,i] - prediction)
}
