#' Prediction Simulation Study
#'
#' Replicate the results of our simulations part regarding prediction accuracy.
#' @param path Path to export the results.
#' @param n Number of observations.
#' @param pList Vector of number of covariates for each scenario.
#' @param coefConfig Either 'big' or 'small' or c('big', 'small') for
#' the scenarios corresponding to big or small coefficients.
#' @param rho Correlation strength of the Toeplitz covariance matrix.
#' @param nSim Number of simulations.
#' @param cvfolds Number of folds for the cross-validation.
#' @param parallel Should parallel computing be used when possible?
#' @param ncores Number of cores for parallel computing.
#' @param seed Seed for replication purposes.
#' @export
#' @keywords PredictingBlackSwans
predictionSim <- function(path         = getwd(),
                          n            = 100,
                          pList        = c(80, 150, 500),
                          coefConfig   = c("big", "small"),
                          rho          = 0.9,
                          nSim         = 100,
                          cvfolds      = 5,
                          parallel     = TRUE,
                          ncores       = 2L,
                          seed         = 182
) {

  # Default arguments
  propY <- 0.5 # Proportion of 1's
  nVarInf <- 20 # Number of variables for the infeasible LR model

  # Message: Begin of the simulations
  startSim <- Sys.time()
  print(paste("Begin of the simulations:", startSim))

  # Save current directory
  curDir <- getwd()

  # Create folder structure if it does not exist
  simPath <- file.path(path, "Prediction_Simulations")
  simTex <- file.path(simPath, "TeX")
  simHtml <- file.path(simPath, "Html")
  simPng <- file.path(simPath, "Png")
  simPdf <- file.path(simPath, "Pdf")

  dir.create(simPath, showWarnings = FALSE)
  dir.create(simTex, showWarnings = FALSE)
  dir.create(simHtml, showWarnings = FALSE)
  dir.create(simPng, showWarnings = FALSE)
  dir.create(simPdf, showWarnings = FALSE)

  # Set working directory for exporting the graphs
  setwd(simPath)

  # Register kernels for parallel computing if parallel = TRUE
  if (parallel == TRUE) doMC::registerDoMC(cores=ncores)

  # Initialize iterators and containers
  predictionList <- list()
  # Loop over p
  for (kk in 1:length(pList)) {
    p <- pList[kk]
    print(paste("Variable combination:", kk, "of", length(pList)))

    #####################################################################################
    ### Define the data generating process ###
    #####################################################################################
    # Set the seed for reproducibility
    set.seed(seed)

    # Correlation matrix
    covMat <- toeplitz(rho^seq(0, (p-1)))

    # Design matrix: Centered Regressors
    X <- MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=covMat)

    # Standardized the predictor matrix
    # Mean zero and Variance 1
    X <- base::scale(X, center=TRUE, scale=TRUE)

    # Loop over coefficient configurations
    for (mm in 1:length(coefConfig)) {
      coefConfigSet <- coefConfig[mm]
      print(paste("Coefficient configuration combination:", mm, "of", length(coefConfig)))

      # Define the non-zero coefficients
      if (coefConfigSet == "big") {
        nzCoef <- c(1, -1, 2, -2, 3, -3)
      } else if (coefConfigSet == "small") {
        nzCoef <- c(0.2, -0.2, 0.5, -0.5, 0.8, -0.8)
      } else {
        stop("Error: Only two 'small' or 'big' are allowed as coefficient configuration.")
      }

      # Define the parameter vector
      betaVec <- numeric(p)
      posNz <- c(1, 10, 20, 30, 50, 70) # Position of the non-zero coefficients
      betaVec[posNz] <- nzCoef

      posInfeas <- c(posNz, seq(2, 9), 11)

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

      # Dataset
      df <- data.frame(X=X)
      yPos <- ncol(df) + 1

      # Initialize containers
      nModels <- 5
      probsArray <- array(NA, c(n, nSim, nModels))
      resultsMat <- matrix(NA, nModels, 3)
      # Loop over simulations
      i <- 1
      while (i <= nSim) {
        if (i %% 10 == 0) print(paste("Simulation Number:", i, "of", nSim))

        # Generate the response from true probabilities
        df$y <- as.factor(rbinom(n=n, size=1, prob=trueProbs))

        # First this one, because if it does not converge because of
        # complete separation, then the simulation is repeated
        # Infeasible Logistic Regression: Know a low-dimensional
        # subset including the variables with non-zero coefficient
        infeasModel <- glm(y ~ ., family="binomial", data=df[c(posInfeas, yPos)])
        # As symptoms for complete separation either the algorithm did not converge
        # or if yes, then the estimated probabilities will numerically either all 1
        # or all 0.
        if (infeasModel$converged == FALSE |
            (round(mean(infeasModel$fitted.values[which(df$y==1)]), 9)==1 &
             round(mean(infeasModel$fitted.values[which(df$y==0)]), 9)==0)) {
          print(paste("Repeating iteration because of complete separation!"))
          next
        }
        probsArray[, i, 2] <- infeasModel$fitted.values

        # Oracle: Knows the variables with non-zero coefficient
        probsArray[, i, 1] <- glm(y ~ ., family="binomial", data=df[c(posNz, yPos)])$fitted.values

        # Lasso model
        lassoModel <- glmnet::cv.glmnet(x        = X,
                                        y        = df$y,
                                        family   = "binomial",
                                        parallel = T,
                                        nfolds   = cvfolds,
                                        alpha    = 1)
        probsArray[, i, 3] <- glmnet::predict.cv.glmnet(lassoModel, X, s="lambda.min", type="response")

        # Ridge model
        ridgeModel <- glmnet::cv.glmnet(x        = X,
                                        y        = df$y,
                                        family   = "binomial",
                                        parallel = T,
                                        nfolds   = cvfolds,
                                        alpha    = 0)
        probsArray[, i, 4] <- glmnet::predict.cv.glmnet(ridgeModel, X, s="lambda.min", type="response")

        # Random forest
        probsArray[, i, 5] <- randomForest::randomForest(x=X, y=df$y)$votes[, 2]
        # probsArray[, i, 5] <- rep(mean(as.numeric(as.character(df$y))), n)
        i <- i + 1 # Update iterator
      }

      # Save the results
      for (s in 1:nModels) {
        # Average Bias
        resultsMat[s, 1] <- PredictingBlackSwans::ABias_fun(probsArray[,,s], trueProbs)
        # Average Variance
        resultsMat[s, 2] <- PredictingBlackSwans::AVar_fun(probsArray[,,s])
        # Average Mean Squared Prediction Error
        resultsMat[s, 3] <- PredictingBlackSwans::AMSPE_fun(probsArray[,,s], trueProbs)
      }

      # Formatted Table
      resultsTable <- as.data.frame(resultsMat)
      names(resultsTable) <- c("Squared ABias", "AVariance", "AMSPE")
      row.names(resultsTable) <- c("Oracle LR",
                                   "Infeasible LR",
                                   "Lasso",
                                   "Ridge",
                                   "Random Forest")
                                   # "Constant")

      # Save results
      fileName <- paste0("PredContest_", coefConfigSet, "_coefficients_p_", p)
      printNotes <- c(paste0("Number of observations: ", n),
                   paste0("Number of variables: ", p),
                   paste0("Coefficient configuration: ", coefConfigSet),
                   paste0("Population matrix: Toeplitz with rho = ", rho),
                   paste0("Number of non-zero coefficients: ", length(nzCoef)),
                   paste0("Number of simulations: ", nSim))
      setwd(simHtml)
      stargazer::stargazer(resultsTable,
                           type="html",
                           out=paste0(fileName, ".html"),
                           title="Prediction Contest",
                           digits=5,
                           summary=FALSE,
                           notes=printNotes)
      setwd(simTex)
      stargazer::stargazer(resultsTable,
                           type="latex",
                           out=paste0(fileName, ".tex"),
                           title="Prediction Contest",
                           digits=5,
                           summary=FALSE,
                           notes=printNotes)

      # Plot results
      df <- resultsTable[-3]
      df$model <- row.names(df)
      mdf <- reshape::melt.data.frame(df, id.vars="model")
      mdf$modelOrd <- factor(mdf$model, levels=c("Oracle LR",
                                                 "Infeasible LR",
                                                 "Lasso",
                                                 "Ridge",
                                                 "Random Forest"))
                                                 # "Constant"))

      if (coefConfigSet == "big") subTitle <- paste0("Number of variables: ", p, "; Big coefficients")
      else subTitle <- paste0("Number of variables: ", p, "; Small coefficients")
      pt <- ggplot(data = mdf, aes(x=modelOrd, y=value, fill=variable)) +
        geom_col() +
        labs(title="AMSPE Decomposition", x="", y="",
             subtitle=subTitle) +
        theme_classic() +
        theme(plot.title=element_text(hjust=0.5),
              plot.subtitle=element_text(hjust=0.5),
              legend.position="bottom",
              legend.title.align=0.5,
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        guides(fill=guide_legend(title="Components", title.position="top")) +
        scale_fill_manual(values=c("grey20", "grey60"))

      setwd(simPng)
      ggsave(filename=paste0(fileName, ".png"),
             plot=pt,
             height=9,
             width=17,
             units="cm")

      setwd(simPdf)
      ggsave(filename=paste0(fileName, ".pdf"),
             plot=pt,
             height=9,
             width=17,
             units="cm")

      # Name for the list
      listName <- paste0("CoefConfig_", coefConfigSet, "_p_", p)

      # Save results in a list
      predictionList[[listName]] <- list(resultsTable = resultsTable,
                                     probsArray   = probsArray,
                                     p            = p,
                                     propY          = propY,
                                     rho      = rho,
                                     plot         = pt)
    }
  }
  # End of the simulations
  endSim <- Sys.time()
  # Time information
  print(paste("Begin of the simulations:", startSim))
  print(paste("End of the simulations:", endSim))
  # Time difference
  print(endSim - startSim)

  # Save data
  timeInfo <- data.frame(begin      = startSim,
                         end        = endSim,
                         timeNeeded = endSim - startSim)
  parInfo <- data.frame(path        = path,
                        n           = n,
                        pList       = pList,
                        nSim        = nSim,
                        propY         = propY,
                        rho = rho,
                        nVarInf     = nVarInf,
                        seed        = seed)
  predSimulationsResults <- list(predictionList = predictionList,
                                 parInfo        = parInfo,
                                 timeInfo       = timeInfo)
  save(predSimulationsResults, file="predSimulationsResults.RData")
  # Return back to the original working directory
  setwd(curDir)
}
