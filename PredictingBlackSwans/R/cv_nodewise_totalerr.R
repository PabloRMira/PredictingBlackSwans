#' Compute the total cross-validated error for the nodewise regression
#'
#' @param c Column of the response in the nodewise regression.
#' @param K Number of folds for the cross-validation.
#' @param x Matrix of predictors.
#' @param lambdas Vector / Sequence of regularization parameters.
#' @return A (Number of lambdas x Number of folds) matrix of
#' cross-validated errors (error on the discarded fold).
#' @export
cv_nodewise_totalerr <- function(c,
                                 dataselects,
                                 x,
                                 lambdas,
                                 K) {
  # Initialize the output
  totalerr <- matrix(nrow = length(lambdas), ncol = K)
  for(i in 1:K){ # Loop over the test sets
    whichj <- dataselects == i # Test part of the data
    glmnetfit <- glmnet::glmnet(x = x[!whichj, -c, drop = FALSE],
                        y = x[!whichj, c, drop = FALSE],
                      lambda = lambdas)
    predictions  <- glmnet::predict.glmnet(glmnetfit,
                                           newx = x[whichj, -c, drop = FALSE],
                                           s = lambdas)
    totalerr[, i] <- colMeans((x[whichj, c] - predictions)^2)
  }
  return(totalerr)
}
