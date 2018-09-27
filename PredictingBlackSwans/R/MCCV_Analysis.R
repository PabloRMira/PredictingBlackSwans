#' Monte Carlo Cross-Validation Analysis for the performance comparison
#'
#' Main function to reproduce the results of our paper regarding the
#' performance comparison by means of Monte Carlo Cross-Validation (MCCV).
#' @param mccvnumber Number of Monte Carlo Cross-Validation runs.
#' @param cvfolds Number of folds for the cross-validation with the Lasso.
#' @param seed Seed for reproducibility.
#' @param parallel Should paralle computing be used? (Note: Only for UNIX computers).
#' @param ncores Number of cores.
#' @export
#' @keywords PredictingBlackSwans
MCCV_Analysis <- function(mccvnumber=100,
                          cvfolds=5,
                          seed=813,
                          parallel=T,
                          ncores=2L) {
  # Seed
  set.seed(seed)

  # Prepare output
  outputTable <- data.frame(Model=c("LR", "Lasso", "Random Forest"), # Model
                            AUC=character(3), # AUC
                            ConfInt=character(3), # 95% Confidence interval for AUC
                            N=numeric(3), # Number of observations
                            Ncrises=numeric(3), # Number of crisis observations
                            Nvars=numeric(3), # Number of variables
                            stringsAsFactors = F)
  names(outputTable) <- c("Model", "AUC", "95% CI", "No. of observations", "No. of crises",
                          "No. of variables")

  ##################################################################################
  ### Data Preparation ###
  ##################################################################################

  # Load the data
  dat <- as.data.table(PredictingBlackSwans::JST_rawData)
  # Remove redundant identifiers
  dat$iso <- NULL
  dat$ifs <- NULL

  # Create variable for the original crisis indicator
  dat[, oricri:=crisisJST]

  # Change codification of crises:
  # Old: 0 <- No crisis; 1 <- Crisis
  # New: 0 <- No crisis; 1 <- Pre-Crisis; 2 <- Crisis
  dat[crisisJST==1, crisisJST:=2]

  # Set the Pre-Crisis period before the crisis period
  dat[, crisislead:=shift(crisisJST, type="lead"), by=country]
  dat[crisisJST==0 & crisislead==2,
      crisisJST:=1L, by=country]
  dat[, crisislead:=NULL]

  prePeriod <- 2 # Period before a crisis to forecast
  for (i in 1:(prePeriod-1)) {
    dat[, crisislead:=shift(crisisJST, type="lead"), by=country]
    dat[crisisJST==0 & crisislead==1, crisisJST:=1L, by=country]
    dat[, crisislead:=NULL]
  }

  # Extend the crisis period since the crisis dates may be very rough
  # and the crisis dynamics would then be mixed with the non-crisis dynamics
  m <- 3 - 1 # Minimum crisis length: 3-years
  for (i in 1:m) {
    dat[, auxvar:=shift(oricri, n=i, type="lag"), by=country]
    dat[is.na(auxvar), auxvar:=0]
    dat[auxvar==1, auxvar:=2]
    dat[, crisisJST:=crisisJST + auxvar]
    dat[, auxvar:=NULL]
  }
  dat[, oricri:=NULL]

  # Exclude observations between the two world wars because of extreme outliers in most countries
  dat <- dat[year < 1914 | (year > 1920 & year < 1938) | year > 1946]

  ##############################################################################################################
  ### LR Analysis ###
  ##############################################################################################################
  # Similar specification as in Ward
  fdat <- copy(dat)

  ### Generate new variables ###
  fdat[, tloans_gdp:=tloans / gdp] # Total loans to GDP
  fdat[, narrowm_gdp:=narrowm / gdp] # Narrow money to GDP

  # Log-differences (~ growth rates for small values)
  varNames <- c("tloans_gdp", "narrowm_gdp", "debtgdp", "gdp", "cpi")
  fdat <- PredictingBlackSwans::logDiffFun(fdat, varNames)

  # Interaction terms (similar as in Ward (2017) but without HP-filter)
  fdat[, pubburden:=debtgdp_lDiff * ltrate] # Public burden
  fdat[, privburden:=tloans_gdp_lDiff * ltrate] # Private burden
  fdat[, jointburden:=tloans_gdp_lDiff * ltrate * debtgdp_lDiff] # Joint burden
  fdat[, tloans_gdp_i_gdp_lDiff:=tloans_gdp * gdp_lDiff]
  fdat[, debtgdp_i_gdp_lDiff:=debtgdp * gdp_lDiff]

  # Variable selection
  varSel <- c("tloans_gdp_lDiff", "debtgdp_lDiff", "narrowm_gdp_lDiff", "ltrate", "gdp_lDiff", "cpi_lDiff",
              "pubburden", "privburden", "jointburden", "tloans_gdp_i_gdp_lDiff", "debtgdp_i_gdp_lDiff")
  fdat <- fdat[, c(varSel, "country", "year", "crisisJST"), with=F]

  # Remove NAs
  fdat <- na.omit(fdat)
  # Remove Infinities
  datAux <- fdat[, -c("crisisJST", "country", "year")]
  fdat <- fdat[is.finite(rowSums(datAux)), ]

  # MCCV
  ci <- c(0.99, 0.95, 0.9) # Confidence intervals
  n.ci <- length(ci)

  # Keep only indicator when it is either 0 or 1 (No crisis - Pre-crisis)
  fdat <- fdat[crisisJST %in% c(0, 1)]
  fdat[, crisisJST:=as.factor(crisisJST)] # Convert the response to factor
  fdat[, rowNumbers:=1:nrow(fdat)] # Get the row numbers for the MCCV
  fmla <- as.formula(paste("crisisJST~", paste(varSel, collapse="+"))) # Formula for the LR model
  aucVec <- numeric(mccvnumber) # Container for the AUC
  ciArray <- array(NA, c(mccvnumber, 2, n.ci)) # Dims: 1: MCCV-iteration; 2: Lower / Upper CI; 3: 99% / 95% / 90%
  rocDTlist <- list() # Container for the ROC-Curves
  print("LR MCCV Analysis starts")
  for (j in 1:mccvnumber) {
    if (j %% 10 == 0) print(paste0("Iteration: ", j, " of ", mccvnumber))
    set.seed(j*seed) # For reproducibility
    # Stratified sampling without replacement
    indexes0 <- sample(fdat[crisisJST==0, rowNumbers], size=0.632*nrow(fdat[crisisJST==0]), replace=F) # Indexes for 0's
    indexes1 <- sample(fdat[crisisJST==1, rowNumbers], size=0.632*nrow(fdat[crisisJST==1]), replace=F) # Indexes for 1's
    train <- fdat[rowNumbers %in% c(indexes0, indexes1)] # Train sample
    test <- fdat[!(rowNumbers %in% c(indexes0, indexes1))] # Test sample
    logitMod <- glm(fmla, data=train, family="binomial") # Logistic Regression
    preds <- predict(logitMod, newdata=test, type="response") # Predictions
    rocAux <- pROC::roc(response=test[, crisisJST], predictor=preds, ci=T, direction="<") # ROC / AUC
    rocDTaux <- data.table(x=(1 - rocAux$specificities), y=rocAux$sensitivities) # Get ROC
    rocDTaux[, MCCViter:=j] # Index for the MCCV iteration
    rocDTlist[[j]] <- rocDTaux # Save the ROC curve
    aucVec[j] <- as.numeric(rocAux$auc) # AUC
    ciArray[j, , 1] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[1]))[c(1, 3)]  # 99% CI for AUC
    ciArray[j, , 2] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[2]))[c(1, 3)]  # 95% CI for AUC
    ciArray[j, , 3] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[3]))[c(1, 3)]  # 90% CI for AUC
  }
  # Average confidence intervals
  ciMat <- matrix(NA, 3, 2)
  ciMat[1, ] <- colMeans(ciArray[,,1])
  ciMat[2, ] <- colMeans(ciArray[,,2])
  ciMat[3, ] <- colMeans(ciArray[,,3])
  aucStar <- mean(aucVec)
  # Formatting for Latex
  for (j in 1:3) {
    if (0.5 < ciMat[j, 1] & j < 3) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
      break
    }
    if (j == 3 & 0.5 < ciMat[j, 1]) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
    } else {
      aucStar <- paste0(round(aucStar, 2))
    }
  }
  outputTable[1, 2] <- aucStar
  outputTable[1, 3] <- paste0("[", round(ciMat[2, 1], 2), ", ", round(ciMat[2, 2], 2), "]")
  outputTable[1, 4] <- nrow(fdat)
  outputTable[1, 5] <- nrow(fdat[crisisJST==1])
  outputTable[1, 6] <- length(varSel)

  # ROC data for all Monte Carlo Cross Validations
  rocDT <- rbindlist(rocDTlist)
  rocDT <- rocDT[order(x, y)]
  # Calculate the mean roc curve and standard deviation across Cross-Validations
  xroc <- seq(0, 1, length=100)
  yrocmat <- matrix(NA, 100, mccvnumber)
  for (i in 1:mccvnumber) {
    rocInterp <- approxfun(x=rocDT[MCCViter==i, x],
                           y=rocDT[MCCViter==i, y], method="linear")
    yrocmat[, i] <- rocInterp(xroc)
  }
  meanRocLogit <- data.table(x=xroc, y=rowMeans(yrocmat))
  meanRocLogit[, Model:="LR"]

  # Find the model with nearest AUC to the mean AUC as a reference model
  # in the MCCV loop to test for AUC equality against the Lasso
  posRef <- which.min((aucVec - mean(aucVec))^2)
  set.seed(posRef*seed) # Seed to find the model with nearest AUC to the mean AUC
  # Stratified sampling without replacement
  indexes0 <- sample(fdat[crisisJST==0, rowNumbers], size=0.632*nrow(fdat[crisisJST==0]), replace=F) # Indexes for 0's
  indexes1 <- sample(fdat[crisisJST==1, rowNumbers], size=0.632*nrow(fdat[crisisJST==1]), replace=F) # Indexes for 1's
  train <- fdat[rowNumbers %in% c(indexes0, indexes1)] # Train sample
  test <- fdat[!(rowNumbers %in% c(indexes0, indexes1))] # Test sample
  logitMod <- glm(fmla, data=train, family="binomial") # Logistic Regression
  preds <- predict(logitMod, newdata=test, type="response") # Predictions
  rocRefLogit <- pROC::roc(response=test[, crisisJST], predictor=preds, ci=T, direction="<") # Reference ROC

  # Housekeeping
  rm(ciMat, ciArray, rocAux, rocDTaux, rocDTlist)

  ##############################################################################################################
  ### Random Forest ###
  ##############################################################################################################
  fdat <- copy(dat)

  # First selection of raw variables
  varsFirst <- c("rconpc", "rgdpmad", "pop", "gdp", "imports", "exports", "ltrate", "stir", "tmort",
                 "debtgdp", "revenue", "expenditure", "tloans", "hpnom", "narrowm", "cpi",
                 "thh", "tbus", "iy", "stocks")
  fdat <- fdat[, c("crisisJST", "country", "year", varsFirst), with=F]

  # Take log-differences
  logDiffVars <- setdiff(varsFirst, c("ltrate", "stir", "ca", "revenue", "exports", "imports",
                                      "expenditure", "narrowm"))
  fdat <- PredictingBlackSwans::logDiffFun(fdat, logDiffVars)

  # Create new variables / features
  fdat[, tloans_gdp:=tloans / gdp] # Total loans to the private sector to GDP
  fdat[, tmort_gdp:=tmort / gdp] # Total mortgages to GDP
  fdat[, thh_gdp:=thh / gdp] # Total credit to households to GDP
  fdat[, tbus_gdp:=tbus / gdp] # Total credit to businesses to GDP
  fdat[, stocks_gdp:=stocks / gdp] # Stock prices to GDP

  fdat[, tloans_cpi:=tloans / cpi] # Real total credit to the private sector
  fdat[, tmort_cpi:=tmort / cpi] # Real total mortgages
  fdat[, thh_cpi:=thh / cpi] # Real total credit to households
  fdat[, tbus_cpi:=tbus / cpi] # Real total credit to businesses
  fdat[, stocks_cpi:=stocks / cpi] # Real stock prices

  fdat[, netexport:=exports - imports] # Net export
  fdat[, govbalance:=revenue - expenditure] # Government balance
  fdat[, ratediff:=ltrate - stir] # Long-term - Short-term interest rate

  fdat[, netexport_gdp:=netexport / gdp] # Net export to GDP
  fdat[, govbalance_gdp:=govbalance / gdp] # Government balance to GDP

  ### Detrend new features ###
  finVars <- c("tloans", "tmort", "thh", "tbus", "stocks") # Financial variables
  ratioVars <- c("gdp", "cpi") # Variables for the ratios (e.g. "to GDP", real variables)
  # Detrend by taking log-differences
  logDiffVars <- paste0(rep(finVars, each=2), "_", rep(ratioVars, length(finVars)))
  fdat <- PredictingBlackSwans::logDiffFun(fdat, logDiffVars)

  # Detrend by taking simple differences
  newFeatures <- c("netexport", "govbalance")
  diffVars <- c("ltrate", "stir", paste0(newFeatures, "_", "gdp"))
  fdat <- PredictingBlackSwans::diffFun(fdat, diffVars)

  # Remove no longer needed variables because of being rather non-stationary
  # Note: Keep investment to GDP and the interest rates because of these variables
  # being naturally / economically bounded and therefore stationary
  exVars <- c(newFeatures, paste0(newFeatures, "_", "gdp"),
              logDiffVars, setdiff(varsFirst, c("iy", "ltrate", "stir")))
  fdat <- fdat[, -c(exVars), with=F]

  # Take lags of all variables in the selection
  lagVars <- setdiff(names(fdat), c("crisisJST", "country", "year"))
  fdat <- PredictingBlackSwans::lagFun(fdat, inputVars=lagVars, maxLag=2)

  # Selection
  varSel <- setdiff(names(fdat), c("crisisJST", "country", "year"))

  # Remove NAs
  fdat <- na.omit(fdat)
  # Remove Infinities
  datAux <- fdat[, -c("crisisJST", "country", "year")]
  fdat <- fdat[is.finite(rowSums(datAux)), ]

  fdat <- fdat[crisisJST %in% c(0, 1)] # Remove crisis periods: only non-crises and precrises
  fdat[, crisisJST:=as.factor(crisisJST)] # Response as factor

  # Random forest
  print("Random Forest OOB Analysis starts")
  strataSampling <- as.vector(ceiling(table(fdat$crisisJST) * 0.632)) # Stratified Sampling
  rfMod <- randomForest::randomForest(x=as.matrix(fdat[, varSel, with=F]),
                                      y=fdat$crisisJST,
                                      strata=fdat$crisisJST,
                                      sampsize=strataSampling)
  preds <- predict(rfMod, type="prob")[,2] # Out-of-bag predictions
  rocRf <- pROC::roc(fdat$crisisJST, preds, ci=T, direction="<") # Out-of-bag ROC analysis

  ciMat <- matrix(NA, 3, 2)
  ciMat[1, ] <- as.numeric(pROC::ci.auc(rocRf, conf.level=ci[1]))[c(1, 3)] # 99%
  ciMat[2, ] <- as.numeric(pROC::ci.auc(rocRf, conf.level=ci[2]))[c(1, 3)] # 95%
  ciMat[3, ] <- as.numeric(pROC::ci.auc(rocRf, conf.level=ci[3]))[c(1, 3)] # 90%
  aucStar <- as.numeric(pROC::auc(rocRf))
  for (j in 1:3) {
    if (0.5 < ciMat[j, 1] & j < 3) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
      break
    }
    if (j == 3 & 0.5 < ciMat[j, 1]) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
    } else {
      aucStar <- paste0(round(aucStar, 2))
    }
  }
  outputTable[3, 2] <- aucStar
  outputTable[3, 3] <- paste0("[", round(ciMat[2, 1], 2), ", ", round(ciMat[2, 2], 2), "]")
  outputTable[3, 4] <- nrow(fdat)
  outputTable[3, 5] <- nrow(fdat[crisisJST==1])
  outputTable[3, 6] <- length(varSel)

  # Mean ROC
  meanRocRf <- data.table(x=(1 - rocRf$specificities), y=rocRf$sensitivities)
  meanRocRf[, Model:="Random Forest"]

  #########################################################################################
  ### Lasso High-Dimensional ###
  #########################################################################################
  # High-dimensional quadratic function
  fdat <- copy(dat)

  # First selection of raw variables
  varsFirst <- c("rconpc", "rgdpmad", "pop", "gdp", "imports", "exports", "ltrate", "stir", "tmort",
                 "debtgdp", "revenue", "expenditure", "tloans", "hpnom", "narrowm", "cpi",
                 "thh", "tbus", "iy", "stocks")
  fdat <- fdat[, c("crisisJST", "country", "year", varsFirst), with=F]

  # Take log-differences
  logDiffVars <- setdiff(varsFirst, c("ltrate", "stir", "ca", "revenue", "exports", "imports",
                                      "expenditure", "narrowm"))
  fdat <- PredictingBlackSwans::logDiffFun(fdat, logDiffVars)

  # Create new variables / features
  fdat[, tloans_gdp:=tloans / gdp] # Total loans to the private sector to GDP
  fdat[, tmort_gdp:=tmort / gdp] # Total mortgages to GDP
  fdat[, thh_gdp:=thh / gdp] # Total credit to households to GDP
  fdat[, tbus_gdp:=tbus / gdp] # Total credit to businesses to GDP
  fdat[, stocks_gdp:=stocks / gdp] # Stock prices to GDP

  fdat[, tloans_cpi:=tloans / cpi] # Real total credit to the private sector
  fdat[, tmort_cpi:=tmort / cpi] # Real total mortgages
  fdat[, thh_cpi:=thh / cpi] # Real total credit to households
  fdat[, tbus_cpi:=tbus / cpi] # Real total credit to businesses
  fdat[, stocks_cpi:=stocks / cpi] # Real stock prices

  fdat[, netexport:=exports - imports] # Net export
  fdat[, govbalance:=revenue - expenditure] # Government balance
  fdat[, ratediff:=ltrate - stir] # Long-term - Short-term interest rate

  fdat[, netexport_gdp:=netexport / gdp] # Net export to GDP
  fdat[, govbalance_gdp:=govbalance / gdp] # Government balance to GDP

  ### Detrend new features ###
  finVars <- c("tloans", "tmort", "thh", "tbus", "stocks") # Financial variables
  ratioVars <- c("gdp", "cpi") # Variables for the ratios (e.g. "to GDP", real variables)
  # Detrend by taking log-differences
  logDiffVars <- paste0(rep(finVars, each=2), "_", rep(ratioVars, length(finVars)))
  fdat <- PredictingBlackSwans::logDiffFun(fdat, logDiffVars)

  # Detrend by taking simple differences
  newFeatures <- c("netexport", "govbalance")
  diffVars <- c("ltrate", "stir", paste0(newFeatures, "_", "gdp"))
  fdat <- PredictingBlackSwans::diffFun(fdat, diffVars)

  # Remove no longer needed variables because of being rather non-stationary
  # Note: Keep interest rates and the investment to GDP because of these variables
  # being naturally bounded and therefore stationary
  exVars <- c(newFeatures, paste0(newFeatures, "_", "gdp"),
              logDiffVars, setdiff(varsFirst, c("iy", "ltrate", "stir")))
  fdat <- fdat[, -c(exVars), with=F]

  # Interactions. Keeping it simple to not overload the model because as shown in the
  # prediction simulation study, the Lasso get worse with too much variables
  # Therefore, we only build interactions of the credit variables and the interest rates
  # with each other
  interVars <- setdiff(names(fdat), c("crisisJST", "country", "year", "iy", "rconpc_lDiff",
                                      "rgdpmad_lDiff", "pop_lDiff", "iy_lDiff",
                                      "netexport_gdp_diff", "govbalance_gdp_diff"))
  fdat <- interactionsFun(inputData=fdat, inputVars=interVars)

  # Take lags of all variables in the selection
  lagVars <- setdiff(names(fdat), c("crisisJST", "country", "year"))
  fdat <- PredictingBlackSwans::lagFun(inputData=fdat, inputVars=lagVars, maxLag=2)

  # Selection
  varSel <- setdiff(names(fdat), c("crisisJST", "country", "year"))

  # Remove NAs
  fdat <- na.omit(fdat)
  # Remove Infinities
  datAux <- fdat[, -c("crisisJST", "country", "year")]
  fdat <- fdat[is.finite(rowSums(datAux)), ]

  linEffectsDT <- copy(fdat) # Data.table for the linear effects produced later
  fdat <- fdat[crisisJST %in% c(0, 1)] # Remove crisis periods
  fdat[, crisisJST:=as.factor(crisisJST)] # Response as factor
  fdat[, rowNumbers:=1:nrow(fdat)] # Get the row numbers for the MCCV
  aucVec <- numeric(mccvnumber) # Container for the AUC
  ciArray <- array(NA, c(mccvnumber, 2, n.ci)) # Dims: 1: MCCV-iteration; 2: Lower / Upper CI; 3: 99% / 95% / 90%
  rocDTlist <- list() # Container for the ROC-Curves
  if (parallel == T) doMC::registerDoMC(cores=ncores)
  print("Lasso MCCV Analysis")
  for (j in 1:mccvnumber) {
    print(paste0("MCCV Iteration: ", j, " of ", mccvnumber))
    set.seed(j*seed) # For reproducibility
    # Stratified sampling without replacement
    indexes0 <- sample(fdat[crisisJST==0, rowNumbers], size=0.632*nrow(fdat[crisisJST==0]), replace=F) # Indexes for 0's
    indexes1 <- sample(fdat[crisisJST==1, rowNumbers], size=0.632*nrow(fdat[crisisJST==1]), replace=F) # Indexes for 1's
    train <- fdat[rowNumbers %in% c(indexes0, indexes1)] # Train sample
    test <- fdat[!(rowNumbers %in% c(indexes0, indexes1))] # Test sample
    # Stratify the folds for the cross-validation
    dataselects0 <- sample(rep(1:cvfolds, length = nrow(train[crisisJST==0]))) # Folds for 0's
    dataselects1 <- sample(rep(1:cvfolds, length = nrow(train[crisisJST==1]))) # Folds for 1's
    train[crisisJST==0, foldidcv:=dataselects0]
    train[crisisJST==1, foldidcv:=dataselects1]
    lassoMod <- glmnet::cv.glmnet(x=as.matrix(train[, varSel, with=F]),
                                  y=train[, crisisJST],
                                  family="binomial",
                                  foldid=train$foldidcv,
                                  parallel=parallel,
                                  alpha=1)
    preds <- glmnet::predict.cv.glmnet(lassoMod, newx=as.matrix(test[, varSel, with=F]),
                                       s="lambda.1se", type="response") # Predictions
    rocAux <- pROC::roc(response=test[, crisisJST], predictor=as.vector(preds), ci=T, direction="<")
    rocDTaux <- data.table(x=(1 - rocAux$specificities), y=rocAux$sensitivities)
    rocDTaux[, MCCViter:=j]
    rocDTlist[[j]] <- rocDTaux
    aucVec[j] <- as.numeric(rocAux$auc) # AUC
    ciArray[j, , 1] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[1]))[c(1, 3)]  # 99% CI
    ciArray[j, , 2] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[2]))[c(1, 3)]  # 95% CI
    ciArray[j, , 3] <- as.numeric(pROC::ci.auc(rocAux, conf.level=ci[3]))[c(1, 3)]  # 90% CI
  }
  ciMat <- matrix(NA, 3, 2)
  ciMat[1, ] <- colMeans(ciArray[,,1])
  ciMat[2, ] <- colMeans(ciArray[,,2])
  ciMat[3, ] <- colMeans(ciArray[,,3])
  aucStar <- mean(aucVec)
  for (j in 1:3) {
    if (0.5 < ciMat[j, 1] & j < 3) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
      break
    }
    if (j == 3 & 0.5 < ciMat[j, 1]) {
      aucStar <- paste0(round(aucStar, 2), "$^{", paste(rep("\\star", (4 - j)), collapse=""),
                        "}$", collapse = "")
    } else {
      aucStar <- paste0(round(aucStar, 2))
    }
  }
  outputTable[2, 2] <- aucStar
  outputTable[2, 3] <- paste0("[", round(ciMat[2, 1], 2), ", ", round(ciMat[2, 2], 2), "]")
  outputTable[2, 4] <- nrow(fdat)
  outputTable[2, 5] <- nrow(fdat[crisisJST==1])
  outputTable[2, 6] <- length(varSel)

  # ROC data for all Monte Carlo Cross Validations
  rocDT <- rbindlist(rocDTlist)
  rocDT <- rocDT[order(x, y)]
  # Calculate the mean roc curve and standard deviation across Cross-Validations
  # using an interpolation
  xroc <- seq(0, 1, length=100)
  yrocmat <- matrix(NA, 100, mccvnumber)
  for (i in 1:mccvnumber) {
    rocInterp <- approxfun(x=rocDT[MCCViter==i, x],
                           y=rocDT[MCCViter==i, y], method="linear")
    yrocmat[, i] <- rocInterp(xroc)
  }
  meanRocLasso <- data.table(x=xroc,
                             y=rowMeans(yrocmat))
  meanRocLasso[, Model:="Lasso"]

  # Find the model with nearest AUC to the mean AUC in the MCCV loop as a reference model
  # to test for AUC equality against the Random Forest
  posRef <- which.min((aucVec - mean(aucVec))^2)
  set.seed(posRef*seed)
  # Sample without replacement but stratified
  indexes0 <- sample(fdat[crisisJST==0, rowNumbers], size=0.632*nrow(fdat[crisisJST==0]), replace=F) # Indexes for 0's
  indexes1 <- sample(fdat[crisisJST==1, rowNumbers], size=0.632*nrow(fdat[crisisJST==1]), replace=F) # Indexes for 1's
  train <- fdat[rowNumbers %in% c(indexes0, indexes1)] # Train sample
  test <- fdat[!(rowNumbers %in% c(indexes0, indexes1))] # Test sample
  # Stratify the folds for the cross-validation
  dataselects0 <- sample(rep(1:cvfolds, length = nrow(train[crisisJST==0]))) # Folds for 0's
  dataselects1 <- sample(rep(1:cvfolds, length = nrow(train[crisisJST==1]))) # Folds for 1's
  train[crisisJST==0, foldidcv:=dataselects0]
  train[crisisJST==1, foldidcv:=dataselects1]
  lassoMod <- glmnet::cv.glmnet(x=as.matrix(train[, varSel, with=F]),
                                y=train[, crisisJST],
                                family="binomial",
                                foldid=train$foldidcv,
                                parallel=parallel,
                                alpha=1)
  preds <- glmnet::predict.cv.glmnet(lassoMod, newx=as.matrix(test[, varSel, with=F]),
                                     s="lambda.1se", type="response")
  rocRefLasso <- pROC::roc(response=test[, crisisJST], predictor=as.vector(preds), ci=T, direction="<")

  # Analysis on whole dataset for interpretability analyses
  # Stratify the folds for the cross-validation
  dataselects0 <- sample(rep(1:cvfolds, length = nrow(fdat[crisisJST==0]))) # Folds for 0's
  dataselects1 <- sample(rep(1:cvfolds, length = nrow(fdat[crisisJST==1]))) # Folds for 1's
  fdat[crisisJST==0, foldidcv:=dataselects0]
  fdat[crisisJST==1, foldidcv:=dataselects1]
  lassoMod <- glmnet::cv.glmnet(x=as.matrix(fdat[, varSel, with=F]),
                                y=fdat[, crisisJST],
                                family="binomial",
                                foldid=fdat$foldidcv,
                                parallel=parallel,
                                alpha=1)

  # Variable importance
  lassoCoef <- glmnet::coef.cv.glmnet(lassoMod, s="lambda.1se")
  coefDT <- data.table(coefName=lassoCoef@Dimnames[[1]][lassoCoef@i+1][-1],
                       coefValue=as.vector(lassoCoef)[lassoCoef@i+1][-1])
  # Standardized coefficients
  sds <- sqrt(colMeans(as.matrix(train[, varSel, with=F])^2) - (colMeans(as.matrix(train[, varSel, with=F])))^2)
  coefDT[, sds:=sds[coefName]]
  coefDT[, standardized_coef:=coefValue * sds]
  coefDT[, st_coef_abs:=abs(standardized_coef)]
  coefDT[, Sign:=ifelse(standardized_coef >= 0, "POS", "NEG")]
  pVarImp <- ggplot(data=coefDT, aes(x=reorder(coefName, +st_coef_abs), y=st_coef_abs, fill=Sign, color=Sign)) +
    geom_col(alpha=.75) +
    theme_classic() +
    labs(x="Variable", y="Standardized Coefficients") +
    coord_flip()

  # (Pre-)crisis probabilities for all countries
  testCri <- linEffectsDT[year >= 2000]
  testCri$probs <- predict.cv.glmnet(lassoMod, newx=as.matrix(testCri[, varSel, with=F]),
                                     s="lambda.1se", type="response")
  crisDT <- testCri[, .(country, year, probs)]
  countryVec <- unique(as.character(testCri$country))
  crisProbsList <- list()
  for (i in 1:length(countryVec)) {
    p <- PredictingBlackSwans::plotCrisProb(crisDT, testCri, countryVec[i])
    crisProbsList[[countryVec[i]]] <- p$plot
  }

  # Linear Effects
  linEffectsList <- list()
  testCri <- linEffectsDT[year >= 2000]
  testCri <- testCri[, c("country", "year", "crisisJST", varSel), with=F]
  driverAux <- as.matrix(testCri[, varSel, with=F])
  tmpCoefs <- as.vector(glmnet::coef.glmnet(lassoMod, s="lambda.1se"))
  intercept <- tmpCoefs[1]
  tmpCoefs <- tmpCoefs[-1]
  driverTab <- t(t(driverAux) * tmpCoefs)
  driverTab <- as.data.table(driverTab)
  driverTab[, country:=testCri[, country]]
  driverTab[, year:=testCri[, year]]
  tmpCoefs <- glmnet::coef.glmnet(lassoMod, s="lambda.1se")
  nzCoefs <- tmpCoefs@Dimnames[[1]][tmpCoefs@i[-1] + 1]
  driverTab <- driverTab[, c("country", "year", nzCoefs), with=F]

  countryVec <- unique(as.character(testCri$country))
  for (i in 1:length(countryVec)) {
    p <- PredictingBlackSwans::plotLinEffects(driverTable=driverTab,
                                              fullDT=testCri,
                                              countryName=countryVec[i])
    linEffectsList[[countryVec[i]]] <- p
  }

  # Housekeeping
  rm(ciMat, ciArray, rocAux, rocDTaux, rocDTlist)

  ############################################################################################
  ### Test of equality of AUCs for the models as in Ward (2017) ###
  ############################################################################################

  # Prepare the output
  outputLatex <- outputTable

  # Test of equality of AUCs: LR vs Lasso
  # H_0: AUC_Lasso = AUC_Logit      H_1: AUC_Lasso > AUC_Logit
  testLassoLogit <- pROC::roc.test(rocRefLasso, rocRefLogit, method="delong", alternative="greater")
  pval <- testLassoLogit$p.value

  if (pval < 0.05) {
    outputLatex[2,1] <- paste0(outputLatex[2,1], "$^\\S$")
  }

  # Test of equality of AUCs: Lasso vs Random Forest
  # H_0 AUC_RF = AUC_Lasso         H_1: AUC_RF > AUC_Lasso
  testRfLasso <- pROC::roc.test(rocRf, rocRefLasso, method="delong", alternative="greater")
  pval <- testRfLasso$p.value

  if (pval < 0.05) {
    outputLatex[3,1] <- paste0(outputLatex[3,1], "$^\\S$")
  }

  ############################################################################################
  ### ROC comparison ###
  ############################################################################################
  # Order columns
  meanRocLogit <- meanRocLogit[order(x, y)]
  meanRocRf <- meanRocRf[order(x, y)]
  meanRocLasso <- meanRocLasso[order(x, y)]

  rocComp <- rbindlist(list(meanRocLogit, meanRocRf, meanRocLasso))
  rocComp[, ModOrd:=factor(Model, levels=c("LR", "Lasso", "Random Forest"))]

  # Plot
  pRoc <- ggplot(data=rocComp, aes(x=x, y=y, color=ModOrd)) +
    geom_step() +
    geom_segment(x=0, y=0, xend=1, yend=1, linetype=2, color="black") +
    labs(x="False Positive Rate (FPR)", y="True Positive Rate (TPR)",
         title="") +
    theme_bw() +
    theme(plot.title=element_text(hjust=.5)) +
    guides(color=guide_legend(title = "Model")) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent)

  names(outputLatex)[3] <- "95\\% CI"
  outputList <- list(outputLatex=outputLatex,
                     ROCPlot=pRoc,
                     pVarImp=pVarImp,
                     pCrisProbs=crisProbsList,
                     pLinEffects=linEffectsList)
  return(outputList)
}
