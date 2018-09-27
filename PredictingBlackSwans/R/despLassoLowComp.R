#' Compute the Desparsified Lasso Estimator for a Low Dimensional Component
#'
#' Compute the Desparsified Lasso Estimator for Logistic Regression for a low
#' dimensional subset of variables of interest. The square-root Lasso with
#' Belloni, Chernozhukov & Wang (2011) tuning parameter is used.
#' @param x Matrix of predictors.
#' @param y Response variable.
#' @param cvfolds Number of folds for the cross-validation for the initial
#' estimator.
#' @param lowSel Low dimensional selection of variables of interest for which
#' the desparsified Lasso will yield estimates.
#' @param parallel Should parallel computing be used when possible?
#' @param ncores Number of cores to be used when parallel is TRUE.
#' @export
despLassoLowComp <- function(x,
                             y,
                             cvfolds,
                             lowSel=seq(1, ncol(x)),
                             parallel = TRUE,
                             ncores = getOption("mc.cores", 2L))
{
  # Set default parameters
  standardize <- TRUE # Always standardize the predictors before the computations

  # Multiple correction method: Bonferroni-Holm procedure
  multiplecorr.method <- "holm"

  ###################################################################################
  ### Preprocessing ###
  ###################################################################################

  # Dimensions of the predictor matrix
  n <- nrow(x)
  p <- ncol(x)

  # Convert response to numeric if it is factor
  if (is.factor(y)) y <- as.numeric(as.character(y))

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

  # Stratify the folds because of rarity of the event of a financial crisis
  indexDT <- data.table(y=y)
  dataselects0 <- sample(rep(1:cvfolds, length = nrow(indexDT[y==0]))) # Folds for 0's
  dataselects1 <- sample(rep(1:cvfolds, length = nrow(indexDT[y==1]))) # Folds for 1's
  indexDT[y==0, foldidcv:=dataselects0]
  indexDT[y==1, foldidcv:=dataselects1]
  # Compute the Lasso estimator by cross-validation
  fitnet <- glmnet::cv.glmnet(x, y, family="binomial", standardize=standardize, foldid=indexDT$foldidcv)
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
  # Lambda recommended by Bernoulli, Chernozhukov & Wang (2011)
  Z <- nodewise_sqrtlasso_low_comp(x = x,
                                   lowSel = lowSel,
                                   parallel=T,
                                   ncores=ncores,
                                   lambda="BCW")
  # Discard all other variables
  x <- x[, lowSel]
  # Rescale Z such that Z_j^T X_j = 1 (the denominator) for all j to simplify computations
  scaleZ <- diag(crossprod(Z,x))
  Z      <- scale(Z, center = FALSE, scale = scaleZ)

  ###########################################################################################
  ### Desparsified Lasso Estimator ###
  ###########################################################################################

  # Adjusted response
  yadj <- (y - pihat) / sqrt(diagW)

  # Debiasing term
  debias <- crossprod(Z, yadj)

  # Desparsified Lasso
  dLasso <- betahat[-1][lowSel] + debias

  # Standard errors. Note: As in the hdi software using the "desparsified Hessian"
  se <- sqrt(colSums(Z^2))

  # "Z-Statistic" for the desparsified Lasso
  dLassoRescaled <- dLasso / se

  # Calculate p-values
  pval <- 2 * pnorm(abs(dLassoRescaled), lower.tail = FALSE)

  # Multiple testing correction: Bonferroni-Holm
  pcorr <- p.adjust(pval, method = multiplecorr.method)

  # Only the low dimensional component of standard deviations
  sds <- sds[lowSel]

  # Return information
  dLassoOutput <- list(pval        = as.vector(pval),
                       pcorr       = pcorr,
                       standardize = standardize,
                       sds         = sds,
                       dLasso      = dLasso / sds,
                       se          = se / sds,
                       Zstatistic  = dLassoRescaled)
  return(dLassoOutput)
}
