#' Compute the Desparsified Lasso Estimator
#'
#' Compute the Desparsified Lasso Estimator for Logistic Regression
#' @param x Matrix of predictors.
#' @param y Response variable.
#' @param cvfolds Number of folds for the cross-validation for both the initial
#' estimator and for the CV-nodewise-Lasso (if this is chosen).
#' @param parallel Should parallel computing be used when possible?
#' @param ncores Number of cores to be used when parallel is TRUE.
#' @export
despLasso <- function(x,
                      y,
                      nodewise = c("cv", "sqrtLasso"),
                      lambda = c("BCW", "VdG"),
                      cvfolds = 5,
                      parallel = TRUE,
                      ncores = getOption("mc.cores", 2L))
{
  # Set default parameters
  standardize <- TRUE # Always standardize the predictors before the computations

  # Multiple correction method: Bonferroni-Holm procedure
  multiplecorr.method <- "holm"

  # Stop if error
  if (!all(nodewise %in% c("cv", "sqrtLasso"))) stop(paste("Error: The nodewise option", nodewise, "is not allowed"))

  ###################################################################################
  ### Preprocessing ###
  ###################################################################################

  # Dimensions of the predictor matrix
  n <- nrow(x)
  p <- ncol(x)

  # Get the standard deviation of the predictors
  # NOTE: This are the standard deviation with factor 1 / (n -1)
  # from the default sd function in base R
  if (standardize) sds <- apply(x, 2, sd) else sds <- rep(1, p)

  # Center and scale the predictor matrix to have all predictors with mean 0 and variance 1
  x <- scale(x, center = TRUE, scale = standardize)

  ###################################################################################
  ### Compute the Lasso Estimator ###
  ###################################################################################
  # Compute the Lasso estimator and the weighted predictor matrix

  # Compute the Lasso estimator by cross-validation
  fitnet <- glmnet::cv.glmnet(x, y, family="binomial", standardize=standardize, nfolds=cvfolds)
  # Get the predicted probabilities
  pihat <- as.vector(glmnet::predict.cv.glmnet(object=fitnet, newx=x,
                                               s=fitnet$lambda.min, type="response"))
  # Extract the Lasso coefficients with the lambda minimizing the cross-validated error curve
  betahat <- as.vector(glmnet::predict.cv.glmnet(fitnet, x, s = fitnet$lambda.min,
                                                 type = "coefficients"))
  diagW <- pihat * (1 - pihat) # Get the diagonal of the weight matrix
  # Weighted predictor matrix (X_\beta in van de Geer (2014))
  xw <- sqrt(diagW) * x
  # Center the columns of the weighted predictor matrix to get
  # rid of the intercept in the nodewise lasso regressions
  x <- scale(xw, center = TRUE, scale = FALSE)

  ###########################################################################################
  ### Nodewise Regression ###
  ###########################################################################################
  # Calculate the matrix of residuals of the nodewise regressions Z = (Z_1, ..., Z_p), Z_j \in \R^n
  if (all(nodewise == "cv")) {
    # NOTE: We use a common lambda sequence for all nodewise regressions
    # performing cross-validation on all nodewise regressions
    Z <- PredictingBlackSwans::nodewise_cv(x        = x,
                                           parallel = parallel,
                                           ncores   = ncores,
                                           lambda   = "lambda.min",
                                           cvfolds  = cvfolds)
    # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
    scaleZ <- diag(crossprod(Z,x))
    Z      <- scale(Z, center = FALSE, scale = scaleZ)
  } else if (all(nodewise == "sqrtLasso")) {
    if (length(lambda) == 1) {
      Z <- PredictingBlackSwans::nodewise_sqrtlasso(x = x,
                                                    parallel=T,
                                                    ncores=ncores,
                                                    lambda=lambda)
      # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
      scaleZ <- diag(crossprod(Z,x))
      Z      <- scale(Z, center = FALSE, scale = scaleZ)
    } else if (length(lambda) == 2) {
      # Lambda recommended by Bernoulli, Chernozhukov & Wang (2011)
      Z <- PredictingBlackSwans::nodewise_sqrtlasso(x = x,
                                                    parallel=T,
                                                    ncores=ncores,
                                                    lambda="BCW")
      # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
      scaleZ <- diag(crossprod(Z,x))
      Zbcw      <- scale(Z, center = FALSE, scale = scaleZ)
      # Lambda recommended by Van de Geer (2014)
      Z <- PredictingBlackSwans::nodewise_sqrtlasso(x = x,
                                                    parallel=T,
                                                    ncores=ncores,
                                                    lambda="VdG")
      # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
      scaleZ <- diag(crossprod(Z,x))
      Zvdg      <- scale(Z, center = FALSE, scale = scaleZ)
    } else {
      stop("Error: You have more than two lambda choices")
    }
  } else if (all(c("cv", "sqrtLasso") %in% nodewise)) {
    # CV Lasso
    Z <- PredictingBlackSwans::nodewise_cv(x        = x,
                                           parallel = parallel,
                                           ncores   = ncores,
                                           lambda   = "lambda.min",
                                           cvfolds  = cvfolds)
    # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
    scaleZ <- diag(crossprod(Z,x))
    Zcv      <- scale(Z, center = FALSE, scale = scaleZ)
    # Lambda recommended by Bernoulli, Chernozhukov & Wang (2011)
    Z <- PredictingBlackSwans::nodewise_sqrtlasso(x = x,
                                                  parallel=T,
                                                  ncores=ncores,
                                                  lambda="BCW")
    # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
    scaleZ <- diag(crossprod(Z,x))
    Zbcw      <- scale(Z, center = FALSE, scale = scaleZ)
    # Lambda recommended by Van de Geer (2014)
    Z <- PredictingBlackSwans::nodewise_sqrtlasso(x = x,
                                                  parallel=T,
                                                  ncores=ncores,
                                                  lambda="VdG")
    # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
    scaleZ <- diag(crossprod(Z,x))
    Zvdg      <- scale(Z, center = FALSE, scale = scaleZ)
  }

  ###########################################################################################
  ### Desparsified Lasso Estimator ###
  ###########################################################################################

  # Adjusted response
  yadj <- (y - pihat) / sqrt(diagW)

  if (length(nodewise) == 1) {
    # Debiasing term
    debias <- crossprod(Z, yadj)

    # Desparsified Lasso
    dLasso <- betahat[-1] + debias

    # Standard errors. Note: As in the hdi software using the "desparsified Hessian"
    se <- sqrt(colSums(Z^2))

    # "Z-Statistic" for the desparsified Lasso
    dLassoRescaled <- dLasso / se

    # Calculate p-values
    pval <- 2 * pnorm(abs(dLassoRescaled), lower.tail = FALSE)

    # Multiple testing correction: Bonferroni-Holm
    pcorr <- p.adjust(pval, method = multiplecorr.method)

    # Return information
    dLassoOutput <- list(pval        = as.vector(pval),
                         pcorr       = pcorr,
                         standardize = standardize,
                         sds         = sds,
                         dLasso      = dLasso / sds,
                         se          = se / sds,
                         Zstatistic  = dLassoRescaled)
  } else {
    # With 'cv'-nodewise matrix
    Z <- Zcv
    # Debiasing term
    debias <- crossprod(Z, yadj)

    # Desparsified Lasso
    dLasso <- betahat[-1] + debias

    # Standard errors. Note: As in the hdi software using the "desparsified Hessian"
    se <- sqrt(colSums(Z^2))

    # "Z-Statistic" for the desparsified Lasso
    dLassoRescaled <- dLasso / se

    # Calculate p-values
    pval <- 2 * pnorm(abs(dLassoRescaled), lower.tail = FALSE)

    # Multiple testing correction: Bonferroni-Holm
    pcorr <- p.adjust(pval, method = multiplecorr.method)

    # Return information
    dLassoCV <- list(pval        = as.vector(pval),
                     pcorr       = pcorr,
                     standardize = standardize,
                     sds         = sds,
                     dLasso      = dLasso / sds,
                     se          = se / sds,
                     Zstatistic  = dLassoRescaled)

    # With 'sqrtLasso'-nodewise matrix and lambda recommended by BCW
    Z <- Zbcw
    # Debiasing term
    debias <- crossprod(Z, yadj)

    # Desparsified Lasso
    dLasso <- betahat[-1] + debias

    # Standard errors. Note: As in the hdi software using the "desparsified Hessian"
    se <- sqrt(colSums(Z^2))

    # "Z-Statistic" for the desparsified Lasso
    dLassoRescaled <- dLasso / se

    # Calculate p-values
    pval <- 2 * pnorm(abs(dLassoRescaled), lower.tail = FALSE)

    # Multiple testing correction: Bonferroni-Holm
    pcorr <- p.adjust(pval, method = multiplecorr.method)

    # Return information
    dLassoSqrtBCW <- list(pval        = as.vector(pval),
                          pcorr       = pcorr,
                          standardize = standardize,
                          sds         = sds,
                          dLasso      = dLasso / sds,
                          se          = se / sds,
                          Zstatistic  = dLassoRescaled)

    # With 'sqrtLasso'-nodewise matrix and lambda recommended by Van de Geer
    Z <- Zvdg
    # Debiasing term
    debias <- crossprod(Z, yadj)

    # Desparsified Lasso
    dLasso <- betahat[-1] + debias

    # Standard errors. Note: As in the hdi software using the "desparsified Hessian"
    se <- sqrt(colSums(Z^2))

    # Calculate p-values
    pval <- 2 * pnorm(abs(dLassoRescaled), lower.tail = FALSE)

    # Multiple testing correction: Bonferroni-Holm
    pcorr <- p.adjust(pval, method = multiplecorr.method)

    # Return information
    dLassoSqrtVdG <- list(pval        = as.vector(pval),
                          pcorr       = pcorr,
                          standardize = standardize,
                          sds         = sds,
                          dLasso      = dLasso / sds,
                          se          = se / sds,
                          Zstatistic  = dLassoRescaled)

    # Return complete information
    dLassoOutput <- list(dLassoCV      = dLassoCV,
                         dLassoSqrtBCW = dLassoSqrtBCW,
                         dLassoSqrtVdG = dLassoSqrtVdG)
  }
  return(dLassoOutput)
}
