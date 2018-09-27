#' Compute the nodewise Lasso using K-fold cross-validation
#'
#' Compute the nodewise Lasso using K-fold cross-validation and output
#' the matrix of residuals Z stemming from the nodewise regressions.
#' @param x Predictor matrix.
#' @param parallel Should parallel computing be used?
#' @param ncores Number of cores for parallel computing.
#' @param lambda Either "lambda.1se" or "lambda.min" defined as
#' in the glmnet-Package. See e.g. ?glmnet::cv.glmnet.
#' @param cvfolds Number of folds for the cross-validation.
#' @return The matrix of residuals of the nodewise Lasso regressions, i.e.
#' \eqn{Z = (Z_1, ..., Z_p)} with \eqn{Z_j \in R^n}.
#' @export
nodewise_cv <- function(x,
                        parallel = TRUE,
                        ncores   = 2L,
                        lambda   = "lambda.min",
                        cvfolds  = 5) {

  # Get dimensions of data
  n <- nrow(x)
  p <- ncol(x)

  # Number of folds for cross-validation
  K <- cvfolds

  # Get a sequence of lambda values common for all nodewise regressions
  # For this we use the algorithm given in the package glmnet for finding
  # the sequences for each nodewise regression and take a sequence using
  # 100 values (percentiles) of the combined lambda sequences

  # Number of lambda values as in glmnet
  nlambda <- 100
  # Get ratio between maximum lambda and minimum lambda as in glmnet
  lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
  # Get the lambda.max values for all nodewise regressions
  absProdx <- abs(t(scale(x, scale=apply(x, 2, PredictingBlackSwans::mysd))) %*% x)
  diag(absProdx) <- rep(0, ncol(x))
  lambda.max.vec <- apply(absProdx, 2, max) / n
  # Build the lambda sequences in a matrix
  lambdaMat <- matrix(NA, nlambda, p)
  for (j in 1:ncol(x)) {
    lambdaMat[, j] <- exp(seq(log(lambda.min.ratio * lambda.max.vec[j]),
                              log(lambda.max.vec[j]), length=nlambda))
  }
  lambdas <- quantile(lambdaMat, probs=seq(0, 1, length.out = nlambda))
  # Common lambda sequence for all nodewise regression
  lambdas <- sort(lambdas, decreasing = TRUE)

  # Based on code from cv.glmnet for sampling the data for the folds
  dataselects <- sample(rep(1:K, length = n))

  # Compute the error for each fold for each lambda for each nodewise regression
  if (parallel) {
    totalerr <- parallel::mcmapply(PredictingBlackSwans::cv_nodewise_totalerr,
                         c = 1:p,
                         dataselects = list(dataselects = dataselects),
                         x = list(x = x),
                         lambdas = list(lambdas = lambdas),
                         K = K,
                         mc.cores = ncores,
                         SIMPLIFY = FALSE)
  } else {
    totalerr <- mapply(PredictingBlackSwans::cv_nodewise_totalerr,
                       c = 1:p,
                       K = K,
                       dataselects = list(dataselects = dataselects),
                       x = list(x = x),
                       lambdas = list(lambdas = lambdas),
                       SIMPLIFY = FALSE)
  }

  # Convert into suitable array
  # 1. dim: Lambda; 2. dim: Fold; 3. dim predictor variable
  err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
  # Mean error for each lambda and for each predictor variable
  err.mean   <- apply(err.array, 1, mean)
  # Calculate mean for every lambda x fold combination (= average over p)
  # for every lambda then get the standard errors (over folds)
  err.se     <- apply(apply(err.array, c(1, 2), mean), 1, sd) / sqrt(K)
  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]
  stderr.lambda.min <- err.se[pos.min]
  lambda.1se = max(lambdas[err.mean < (min(err.mean) + stderr.lambda.min)])
  # Take best lambda
  if (lambda == "lambda.min") bestlambda <- lambda.min else
    if (lambda == "lambda.1se") bestlambda <- lambda.1se else
      stop("Error: Lambda option is wrong.")

  # Compute the matrix of residuals of the nodewise regressions with
  # the common lambda computed before by 10-fold cross-validation
  Z <- matrix(NA, n, p)
  if (parallel) {
    Z <- parallel::mcmapply(PredictingBlackSwans::getZresiduals, i = 1:p, x = list(x = x),
                  lambda = bestlambda, mc.cores = ncores)
  } else {
    Z <- mapply(PredictingBlackSwans::getZresiduals, i = 1:p, x = list(x = x),
                lambda = bestlambda)
  }
  return(Z)
}
