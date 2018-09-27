#' Main function for the inference simulation study
#' using the "desparsified outer product of the gradient"
#'
#' Computes all the results for the simulation study using the "desparsified outer product
#' of the gradient" as in the original paper (van de Geer (2014)).
#' @param path Path to export the results.
#' @param n Number of observations.
#' @param p Number of predictors.
#' @param rho Correlation parameter of the Toeplitz covariance matrix.
#' @param nSim Number of simulations.
#' @param nomSize Nominal size (type I error a.k.a. alpha).
#' @param cvfolds Number of folds for the cross-validations.
#' @param seed Seed for replication purposes.
#' @param parallel Should parallel computing be used? Note: It only works for UNIX systems.
#' @param ncores How many cores should be used for parallel computing?
#' @export
#' @keywords PredictingBlackSwans
inferenceSimMainPaper <- function(path        = getwd(),
                                  n           = 100,
                                  p           = 150,
                                  rho         = 0.5,
                                  nSim        = 100,
                                  nomSize     = 0.05,
                                  cvfolds     = 5,
                                  seed        = 182,
                                  parallel    = T,
                                  ncores      = getOption("mc.cores", 2L)
) {

  # Default parameters
  propY <- 0.5 # Expected proportion of 1's in the data.

  # Message: Begin of the simulations
  startSim <- Sys.time()
  print(paste("Begin of the simulations:", startSim))

  # Get folder structure
  simPath <- path
  simTex <- file.path(simPath, "TeX")
  simHtml <- file.path(simPath, "Html")

  # Set working directory for exporting the graphs
  setwd(simPath)

  # Set the seed for reproducibility
  set.seed(seed)

  #####################################################################################
  ### Define the data generating process ###
  #####################################################################################

  # Population covariance matrix
  covMat <- toeplitz(rho^seq(0, (p-1)))

  # Design matrix: Centered Regressors
  X <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=covMat)

  # Standardized the predictor matrix
  # Mean zero and Variance 1
  X <- base::scale(X, center=TRUE, scale=TRUE)

  # Define the non-zero coefficients
  posNz <- c(1, 2, 3, 20, 50, 80)
  nzCoef <- c(1, -0.5, 2, 1, -0.5, 2)
  # Cardinality of the active set
  cardActive <- length(nzCoef)
  # Define the parameter vector
  betaVec <- numeric(p)
  betaVec[posNz] <- nzCoef

  # Define the positions of the subset the infeasible LR knows:
  # 12 variables: 6 with true non-zero coefficient + 6 noise variables
  posInfeas <- c(seq(1, 6), seq(20, 21), seq(50, 51), seq(80, 81))
  posInfeasSim <- which(posInfeas %in% posNz) # Positions for the estimated model in the simulations later
  nVarInf <- length(posInfeas)

  # Choose the intercept such that the proportion
  # of 1's approximates propY
  # Note: Numerical solution
  target <- function(alpha) {
    linFun <- alpha + X[, 1:length(nzCoef)] %*% nzCoef
    trueProbs <- 1 / (1 + exp(-linFun))
    targetValue <- (mean(trueProbs) - propY)^2
    return(targetValue)
  }
  alpha <- nleqslv::nleqslv(x=0, fn=target)$x

  # True Linear function
  # for the mapping to the true probabilities
  linFun <- alpha + X %*% betaVec

  # True probabilities for the Data Generating Process
  trueProbs <- 1 / (1 + exp(-linFun))

  # Initialize containers
  # Number of methods to compare: OracleLR, InfeasibleLR, DLassoCV, DLassoSqrtBCW, DLassoSqrtVdG
  nModels <- 5
  # Array for the coefficients of the estimators
  coefArray <- array(0, c(p, nSim, nModels))
  # Coverage: Check if the true value is in the confidence interval
  coverageCheck <- array(NA, c(p, nSim, nModels))
  # Power of the Oracle Logistic Regression
  PowerLogistic <- matrix(NA, nSim, length(nzCoef))
  # Power of the infeasible logistic regression estimator
  PowerInfeasible <- matrix(NA, nSim, length(nzCoef))
  # Power of the desparsified Lasso
  PowerDLasso <- array(NA, c(nSim, length(nzCoef), 3)) # 3 Desparsified Lasso's
  # Family-wise Error Rate of the Oracle logistic regression
  FWERLogistic <- matrix(NA, nSim, length(nzCoef))
  # Family-wise Error Rate of the infeasible logistic regression estimator
  FWERInfeasible <- matrix(NA, nSim, nVarInf - length(nzCoef))
  # Family-wise Error Rate of the Desparsified Lasso
  FWERDlasso <- array(NA, c(nSim, length(betaVec) - length(nzCoef), 3)) # 3 Desparsified Lasso's
  # Family-wise Error Rate of the Desparsified Lasso for the low-dimensional test
  FWERDlassoLow <- array(NA, c(nSim, nVarInf - length(nzCoef), 3))
  # Power of the desparsified Lasso for the low-dimensional test
  PowerDLassoLow <- array(NA, c(nSim, length(nzCoef), 3))
  # Save also the standard errors
  seArray <- array(0, c(p, nSim, nModels))
  # Data
  df <- data.frame(X=X)
  yPos <- ncol(df) + 1
  # Loop over simulations
  i <- 1
  while (i <= nSim) {
    print(paste("Simulation Number:", i, "of", nSim))

    # Generate the response from true probabilities
    df$y <- rbinom(n=n, size=1, prob=trueProbs)

    #########################################################################################
    ### Infeasible Logistic Regression Estimator ###
    #########################################################################################
    # Infeasible Logistic Regression: Knows a low-dimensional
    # subset including the variables with non-zero coefficient

    # Note: First this one, because if it does not converge because of
    # complete separation, then the simulation is repeated

    # Estimate logistic regression model
    model <- glm(y ~ ., family="binomial", data=df[c(posInfeas, yPos)])

    # Check if there is complete separation and if yes, then repeat the iteration
    if (model$converged == FALSE |
        (round(mean(model$fitted.values[which(df$y==1)]), 9)==1 &
         round(mean(model$fitted.values[which(df$y==0)]), 9)==0)) {
      print(paste("Repeating iteration because of complete separation!"))
      next
    }

    # Get coefficients
    coeffi <- as.vector(model$coefficients)[2:(nVarInf+1)]
    coefArray[posInfeas, i, 2] <- coeffi

    # Get standard errors
    stdErrors <- as.vector(summary(model)$coefficients[, 2])[2:(nVarInf+1)]
    seArray[posInfeas, i, 2]

    # Confidence interval for all variables in the true model (true active set)
    confIntLow <- coeffi - qnorm(1 - (nomSize / 2)) * stdErrors
    confIntUp <- coeffi + qnorm(1 - (nomSize / 2)) * stdErrors

    # Check whether the true value is in the confidence interval
    # for each variable
    coverageCheck[posInfeas, i, 2] <- ifelse(confIntLow <= betaVec[posInfeas] &
                                               confIntUp >= betaVec[posInfeas], 1, 0)

    # p-values for the null hypothesis that the coefficient is equal to zero
    pvalues <- 2 * pnorm(abs(coeffi / stdErrors), lower.tail=F)

    # Adjusted p-values with the Bonferroni-Holm procedure
    pvalcorr <- p.adjust(pvalues, method="holm")

    # Check whether the null hypothesis is rejected (H_0: \beta_j = 0)
    rejectCheck <- ifelse(pvalcorr <= nomSize, 1, 0)

    # Power of the Infeasible logistic regression
    PowerInfeasible[i, ] <- rejectCheck[posInfeasSim]

    # Family-wise Error Rate for the Infeasible LR
    FWERInfeasible[i, ] <- rejectCheck[setdiff(seq(1, nVarInf), posInfeasSim)]

    ##########################################################################################
    ### Oracle estimator ###
    ##########################################################################################

    # Estimate logistic regression model
    model <- glm(y ~ ., family="binomial", data=df[c(posNz, yPos)])

    # Get coefficients without intercept
    coeffi <- as.vector(model$coefficients)[2:(length(nzCoef)+1)]
    coefArray[posNz, i, 1] <- coeffi # Save them

    # Get standard errors
    stdErrors <- as.vector(summary(model)$coefficients[, 2])[2:(length(nzCoef)+1)]
    seArray[posInfeas, i, 1]

    # Confidence interval for all variables in the true model (true active set)
    confIntLow <- coeffi - qnorm(1 - (nomSize / 2)) * stdErrors
    confIntUp <- coeffi + qnorm(1 - (nomSize / 2)) * stdErrors

    # Check whether the true value is in the confidence interval for each variable
    coverageCheck[posNz, i, 1] <- ifelse(confIntLow <= nzCoef & confIntUp >= nzCoef, 1, 0)

    # p-values for the null hypothesis that the coefficient is equal to zero
    pvalues <- 2 * pnorm(abs(coeffi / stdErrors), lower.tail=F)

    # Adjusted p-values with the Bonferroni-Holm procedure
    pvalcorr <- p.adjust(pvalues, method="holm")

    # Power of the Oracle
    PowerLogistic[i, ] <- ifelse(pvalcorr <= nomSize, 1, 0)

    # p-values for the hypothesis that the estimated coefficient is equal to the true value
    pvalues <- 2 * pnorm(abs((coeffi - nzCoef) / stdErrors), lower.tail=F)

    # Adjusted p-values with Bonferroni-Holm
    pvalcorr <- p.adjust(pvalues, method="holm")

    # Family-wise Error Rate for the Oracle
    FWERLogistic[i, ] <- ifelse(pvalcorr <= nomSize, 1, 0)

    #########################################################################################
    ### Desparsified Lasso Estimator ###
    #########################################################################################

    dLassoList <- despLassoPaper(x           = X,
                                 y           = df$y,
                                 nodewise    = c("cv", "sqrtLasso"),
                                 lambda      = c("BCW", "VdG"),
                                 cvfolds     = cvfolds,
                                 parallel    = parallel,
                                 ncores      = ncores) # Seed for the cross-validation

    ### Desparsified Lasso with Cross-validation ###
    dLasso <- dLassoList$dLassoCV

    # Get Desparsified Lasso Coefficients
    dCoef <- dLasso$dLasso
    coefArray[, i, 3] <- dCoef # Save them

    # Get standard errors
    stdErrors <- dLasso$se
    seArray[, i, 3] <- stdErrors

    # Confidence intervals for all variables
    confIntLow <- dCoef - qnorm(1 - (nomSize / 2)) * stdErrors
    confIntUp <- dCoef + qnorm(1 - (nomSize / 2)) * stdErrors

    # Check whether the true value is in the confidence interval
    coverageCheck[, i, 3] <- ifelse(confIntLow <= betaVec & confIntUp >= betaVec, 1, 0)

    # Extract corrected p-values
    pcorr <- dLasso$pcorr
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test: Check whether the test (correctly) rejects
    # the null hypothesis that the estimated value is
    # equal to zero for all non-zero coefficients
    PowerDLasso[i, ,1] <- rejectCheck[posNz]

    # Familywise Error Rate: Check how many times the test (erroneously) rejects the null hypothesis
    # that the estimated coefficients are equal to zero for all zero-coefficients
    FWERDlasso[i, ,1] <- rejectCheck[-posNz]

    ### Low-dimensional test ###
    # Extract p-values for the same coefficients as in the infeasible LR
    pval <- dLasso$pval[posInfeas]

    # Correct them with the Bonferroni-Holm procedure
    pcorr <- p.adjust(pval, method="holm")
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test
    PowerDLassoLow[i,,1] <- rejectCheck[posInfeasSim]

    # Family-wise Error Rate
    FWERDlassoLow[i,,1] <- rejectCheck[setdiff(seq(1, nVarInf), posInfeasSim)]

    ### Desparsified Lasso with Square-Root Lasso ###
    ### and lambda given by BCW ###
    dLasso <- dLassoList$dLassoSqrtBCW

    # Get Desparsified Lasso Coefficients
    dCoef <- dLasso$dLasso
    coefArray[, i, 4] <- dCoef # Save them

    # Get standard errors
    stdErrors <- dLasso$se
    seArray[, i, 4] <- stdErrors

    # Confidence intervals for all variables
    confIntLow <- dCoef - qnorm(1 - (nomSize / 2)) * stdErrors
    confIntUp <- dCoef + qnorm(1 - (nomSize / 2)) * stdErrors

    # Check whether the true value is in the confidence interval
    coverageCheck[, i, 4] <- ifelse(confIntLow <= betaVec & confIntUp >= betaVec, 1, 0)

    # Extract corrected p-values
    pcorr <- dLasso$pcorr
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test: Check whether the test (correctly) rejects
    # the null hypothesis that the estimated value is
    # equal to zero for all non-zero coefficients
    PowerDLasso[i, ,2] <- rejectCheck[posNz]

    # Familywise Error Rate: Check how many times the test (erroneously) rejects the null hypothesis
    # that the estimated coefficients are equal to zero for all zero-coefficients
    FWERDlasso[i, ,2] <- rejectCheck[-posNz]

    ### Low-dimensional test ###
    # Extract p-values for the same coefficients as in the infeasible LR
    pval <- dLasso$pval[posInfeas]

    # Correct them with the Bonferroni-Holm procedure
    pcorr <- p.adjust(pval, method="holm")
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test
    PowerDLassoLow[i,,2] <- rejectCheck[posInfeasSim]

    # Family-wise Error Rate
    FWERDlassoLow[i,,2] <- rejectCheck[setdiff(seq(1, nVarInf), posInfeasSim)]

    ### Desparsified Lasso with Square-Root Lasso ###
    ### and lambda given by van de Geer ###
    dLasso <- dLassoList$dLassoSqrtVdG

    # Get Desparsified Lasso Coefficients
    dCoef <- dLasso$dLasso
    coefArray[, i, 5] <- dCoef # Save them

    # Get standard errors
    stdErrors <- dLasso$se
    seArray[, i, 5] <- stdErrors

    # Confidence intervals for all variables
    confIntLow <- dCoef - qnorm(1 - (nomSize / 2)) * stdErrors
    confIntUp <- dCoef + qnorm(1 - (nomSize / 2)) * stdErrors

    # Check whether the true value is in the confidence interval
    coverageCheck[, i, 5] <- ifelse(confIntLow <= betaVec & confIntUp >= betaVec, 1, 0)

    # Extract corrected p-values
    pcorr <- dLasso$pcorr
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test: Check whether the test (correctly) rejects
    # the null hypothesis that the estimated value is
    # equal to zero for all non-zero coefficients
    PowerDLasso[i, ,3] <- rejectCheck[posNz]

    # Familywise Error Rate: Check how many times the test (erroneously) rejects the null hypothesis
    # that the estimated coefficients are equal to zero for all zero-coefficients
    FWERDlasso[i, ,3] <- rejectCheck[-posNz]

    ### Low-dimensional test ###
    # Extract p-values for the same coefficients as in the infeasible LR
    pval <- dLasso$pval[posInfeas]

    # Correct them with the Bonferroni-Holm procedure
    pcorr <- p.adjust(pval, method="holm")
    rejectCheck <- ifelse(pcorr <= nomSize, 1, 0)

    # Power of the test
    PowerDLassoLow[i,,3] <- rejectCheck[posInfeasSim]

    # Family-wise Error Rate
    FWERDlassoLow[i,,3] <- rejectCheck[setdiff(seq(1, nVarInf), posInfeasSim)]
    i <- i + 1 # Update iterator
  }
  # Save also raw tables without formatting
  rawTables <- list()

  # Mean Coverage tables: Average coverage of the active set
  # and average coverage of the non-active set
  mCoverage <- data.frame(OracleLR=c(mean(coverageCheck[posNz,,1]),
                                     0), # Will change this after rounding
                          InfeasibleLR=c(mean(coverageCheck[posNz,,2]),
                                         mean(coverageCheck[setdiff(posInfeas, posNz),,2])),
                          DLassoCV=c(mean(coverageCheck[posNz,,3]),
                                     mean(coverageCheck[setdiff(seq(1, p), posNz),,3])),
                          DLassoSqrtBCW=c(mean(coverageCheck[posNz,,4]),
                                          mean(coverageCheck[setdiff(seq(1, p), posNz),,4])),
                          DLassoSqrtVdG=c(mean(coverageCheck[posNz,,5]),
                                          mean(coverageCheck[setdiff(seq(1, p), posNz),,5])))
  row.names(mCoverage) <- c("Active Set", "Non-active Set")
  rawTables[["mCoverage"]] <- mCoverage
  mCoverage <- round(mCoverage, 2)
  mCoverage[2, 1] <- "" # Since the Oracle does not estimate any coefficient in the non-active set

  # Test results, general
  tResults <- data.frame(OracleLR=c(mean(apply(FWERLogistic, 1, max)),
                                    mean(PowerLogistic)),
                         InfeasibleLR=c(mean(apply(FWERInfeasible, 1, max)),
                                        mean(PowerInfeasible)),
                         DLassoCV=c(mean(apply(FWERDlasso[,,1], 1, max)),
                                    mean(PowerDLasso[,,1])),
                         DLassoSqrtBCW=c(mean(apply(FWERDlasso[,,2], 1, max)),
                                         mean(PowerDLasso[,,2])),
                         DLassoSqrtVdG=c(mean(apply(FWERDlasso[,,3], 1, max)),
                                         mean(PowerDLasso[,,3])))
  row.names(tResults) <- c("FWER", "Power")
  rawTables[["tResults"]] <- tResults
  tResults <- round(tResults, 2)

  # Test results, low-dimensional
  tResultsLow <- data.frame(OracleLR=c(mean(apply(FWERLogistic, 1, max)),
                                       mean(PowerLogistic)),
                            InfeasibleLR=c(mean(apply(FWERInfeasible, 1, max)),
                                           mean(PowerInfeasible)),
                            DLassoCV=c(mean(apply(FWERDlassoLow[,,1], 1, max)),
                                       mean(PowerDLassoLow[,,1])),
                            DLassoSqrtBCW=c(mean(apply(FWERDlassoLow[,,2], 1, max)),
                                            mean(PowerDLassoLow[,,2])),
                            DLassoSqrtVdG=c(mean(apply(FWERDlassoLow[,,3], 1, max)),
                                            mean(PowerDLasso[,,3])))
  row.names(tResultsLow) <- c("FWER", "Power")
  rawTables[["tResultsLow"]] <- tResultsLow
  tResultsLow <- round(tResultsLow, 2)

  # Print results with stargazer in html
  # HTML
  setwd(simHtml)
  stargazer::stargazer(mCoverage,
                       type="html",
                       out=paste0("Coverage_Results_Paper_n_", n, "_p_", p, "_rho_", rho, ".html"),
                       title="Mean Coverage Results",
                       summary=FALSE,
                       notes=c(paste("Number of observations:", n),
                               paste("Number of variables:", p),
                               paste("Correlation strength:", rho),
                               paste("Number of non-zero coefficients:", length(nzCoef)),
                               paste("Nominal size of the test:", nomSize),
                               paste("Average Proportion of 1's in the data:", propY),
                               paste("Number of simulations", nSim)),
                       digits=2)

  stargazer::stargazer(tResults,
                       type="html",
                       out=paste0("Hypothesis_Testing_Paper_n_", n, "_p_", p, "_rho_", rho, ".html"),
                       title="Hypothesis Testing",
                       summary=FALSE,
                       notes=c(paste("Number of observations:", n),
                               paste("Number of variables:", p),
                               paste("Correlation strength:", rho),
                               paste("Number of non-zero coefficients:", length(nzCoef)),
                               paste("Nominal size of the test:", nomSize),
                               paste("Average Proportion of 1's in the data:", propY),
                               paste("Number of simulations", nSim)),
                       digits=2)

  stargazer::stargazer(tResultsLow,
                       type="html",
                       out=paste0("Hypothesis_Testing_LowDim_Paper_n_", n, "_p_", p, "_rho_", rho, ".html"),
                       title="Hypothesis Testing (Low dimensional testing)",
                       summary=FALSE,
                       notes=c(paste("Number of observations:", n),
                               paste("Number of variables:", p),
                               paste("Correlation strength:", rho),
                               paste("Number of non-zero coefficients:", length(nzCoef)),
                               paste("Nominal size of the test:", nomSize),
                               paste("Average Proportion of 1's in the data:", propY),
                               paste("Number of simulations", nSim)),
                       digits=2)

  # Plot sampling distributions
  plotList <- list() # Save the plots
  setwd(simPath)
  methodNames <- c("Oracle LR", "Infeasible LR", "DLasso (CV)", "DLasso (BCW)", "DLasso (VdG)")
  DTauxList <- list()
  for (i in 1:length(methodNames)) {
    DTaux <- as.data.table(coefArray[posNz,,i])
    DTaux[, CoefNumber:=paste0("Coef Number: ", posNz)]
    mDTaux <- melt(DTaux, id.vars="CoefNumber")
    mDTaux[, Model:=methodNames[i]]
    DTauxList[[i]] <- mDTaux
  }
  coefDT <- rbindlist(DTauxList)
  coefDT[, CoefNumber:=factor(CoefNumber, levels = paste0("Coef Number: ", posNz))]
  nzCoefDT <- data.table(CoefNumber=factor(paste0("Coef Number: ", posNz),
                                           levels = paste0("Coef Number: ", posNz)),
                         value=nzCoef)

  pt <- ggplot(data=coefDT, aes(x=value, fill=Model, color=Model)) +
    geom_density(alpha=.35, bw="SJ") +
    geom_vline(data=nzCoefDT, aes(xintercept=value)) +
    facet_wrap(~ CoefNumber) +
    theme_bw() +
    labs(y="Density", x="Coefficients")
  plotList[["AllMethods_NzCoefs"]] <- pt
  ggsave(paste0("InferSim_AllMethods_NonzeroCoefs_Paper_n_", n, "_p_", p, "_rho_", rho, ".png"), plot=pt, units="cm",
         height=10, width=15)

  pt <- ggplot(data=coefDT[Model != "Infeasible LR"], aes(x=value, fill=Model, color=Model)) +
    geom_density(alpha=.35, bw="SJ") +
    geom_vline(data=nzCoefDT, aes(xintercept=value)) +
    facet_wrap(~ CoefNumber) +
    theme_bw() +
    labs(y="Density", x="Coefficients")
  plotList[["WithoutInfeasible_NzCoefs"]] <- pt
  ggsave(paste0("InferSim_WithoutInfeasible_Paper_n_", n, "_p_", p, "_rho_", rho, ".png"), plot=pt,
         units="cm", height=10, width=15)

  # Plot the sampling distribution of the zero-coefficients
  posInfeasCoef <- setdiff(posInfeas, posNz)
  DTauxList <- list()
  for (i in 1:length(methodNames)) {
    DTaux <- as.data.table(coefArray[posInfeasCoef,,i])
    DTaux[, CoefNumber:=paste0("Coef Number: ", posInfeasCoef)]
    mDTaux <- melt(DTaux, id.vars="CoefNumber")
    mDTaux[, Model:=methodNames[i]]
    DTauxList[[i]] <- mDTaux
  }
  coefDT <- rbindlist(DTauxList)
  coefDT[, CoefNumber:=factor(CoefNumber, levels = paste0("Coef Number: ", posInfeasCoef))]
  zCoefDT <- data.table(CoefNumber=factor(paste0("Coef Number: ", posInfeasCoef),
                                          levels = paste0("Coef Number: ", posInfeasCoef)),
                        value=betaVec[posInfeasCoef])
  pt <- ggplot(data=coefDT[Model != "Oracle LR"], aes(x=value, fill=Model, color=Model)) +
    geom_density(alpha=.35, bw="SJ") +
    geom_vline(data=zCoefDT, aes(xintercept=value)) +
    facet_wrap(~ CoefNumber) +
    theme_bw() +
    labs(y="Density", x="Coefficients")
  plotList[["AllMethods_zCoefs"]] <- pt
  ggsave(paste0("InferSim_Zero_Coefficients_Paper_n_", n, "_p_", p, "_rho_", rho, ".png"), plot=pt,
         units="cm", height=10, width=15)

  # End of the simulations
  endSim <- Sys.time()
  # Time information
  print(paste("Begin of the simulations:", startSim))
  print(paste("End of the simulations:", endSim))
  # Time difference
  print(endSim - startSim)

  # Save the data
  timeInfo <- data.frame(begin      = startSim,
                         end        = endSim,
                         timeNeeded = endSim - startSim)
  parInfo <- data.frame(path    = path,
                        n       = n,
                        p       = p,
                        rho     = rho,
                        nSim    = nSim,
                        propY   = propY,
                        nomSize = nomSize,
                        seed    = seed)

  # Save results in a list
  inferenceList <- list()
  listName <- paste0("n_", n, "_p_", p, "_rho_", rho)
  inferenceList[[listName]] <- list(mCoverage   = mCoverage,
                                    tResults    = tResults,
                                    tResultsLow = tResultsLow,
                                    parInfo     = parInfo,
                                    timeInfo    = timeInfo,
                                    plotList    = plotList,
                                    rawTables   = rawTables,
                                    coefArray   = coefArray,
                                    seArray     = seArray,
                                    posNz       = posNz,
                                    betaVec     = betaVec,
                                    posInfeas   = posInfeas)

  saveRDS(inferenceList, file=paste0("Inference_Results_Paper_", listName, ".rds"))
}
