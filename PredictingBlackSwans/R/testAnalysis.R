#' Test Analysis
#'
#' Main function to reproduce the results of our paper regarding the
#' test analysis in the prediction part.
#' @param mccvnumber Number of Monte Carlo Cross-Validation runs.
#' @param cvfolds Number of folds for the cross-validation with the Lasso.
#' @param seed Seed for reproducibility.
#' @param parallel Should parallel computing be used? (Note: Only for UNIX computers)
#' @param ncores Number of cores for parallel computing
#' @export
#' @keywords PredictingBlackSwans
testAnalysis <- function(mccvnumber=100,
                         cvfolds=5,
                         seed=813,
                         parallel=T,
                         ncores=2L) {
  # Seed
  set.seed(seed)

  # Prepare output
  outputTableTest <- data.frame(Model=c("LR", "Lasso", "Random Forest"), # Model
                                AUC=character(3), # AUC
                                N=numeric(3), # Number of observations
                                Ncrises=numeric(3), # Number of crisis observations
                                Nvars=numeric(3), # Number of variables
                                stringsAsFactors = F)
  names(outputTableTest) <- c("Model", "AUC", "No. of observations", "No. of crises", "No. of variables")

  # Begin of test set
  beginTest <- 1998 # To facilitate comparison with Ward (2017)

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

  # Train and test
  train <- fdat[year < beginTest]
  test <- fdat[beginTest <= year]

  # Discard country, year for train sample
  train <- train[, setdiff(names(train), c("country", "year")), with=F]

  # Keep only the pre-crisis indicator when it is either 0 or 1 (No crisis - Pre-crisis)
  train <- train[crisisJST %in% c(0, 1)]
  train[, crisisJST:=as.factor(crisisJST)] # Convert the response to factor
  train[, rowNumbers:=1:nrow(train)] # Get the row numbers for the MCCV
  fmla <- as.formula(paste("crisisJST~", paste(varSel, collapse="+"))) # Formula for the LR model
  ### Model performance on test data ###
  # Train with whole train dataset, test on test dataset (> 1998)
  logitMod <- glm(fmla, data=train, family="binomial") # Logistic Regression
  # Keep only the pre-crisis indicator when it is either 0 or 1 (No crisis - Pre-crisis)
  twt <- test[crisisJST %in% c(0, 1)]
  twt$preds <- predict(logitMod, newdata=twt, type="response")
  rocTest <- pROC::roc(response=twt$crisisJST, predictor=twt$preds, direction="<")
  rocDTtestLogit <- data.table(x=(1 - rocTest$specificities), y=rocTest$sensitivities)
  rocDTtestLogit[, Model:="LR"]
  outputTableTest[1, 2] <- round(as.numeric(rocTest$auc), 2)
  outputTableTest[1, 3] <- nrow(train)
  outputTableTest[1, 4] <- nrow(train[crisisJST==1])
  outputTableTest[1, 5] <- length(varSel)

  probsLogit <- copy(test)
  probsLogit$probs <- predict(logitMod, newdata=probsLogit, type="response")
  probsLogit <- probsLogit[, .(country, year, crisisJST, probs)]
  probsLogit[!(crisisJST %in% c(0, 1)), probs:=NA]

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
  # Note: Keep interest rates and the investment to GDP because of these variables
  # being naturally bounded and therefore stationary
  exVars <- c(newFeatures, paste0(newFeatures, "_", "gdp"),
              logDiffVars, setdiff(varsFirst, c("iy", "ltrate", "stir")))
  fdat <- fdat[, -c(exVars), with=F]

  # Take lags of all variables in the selection
  lagVars <- setdiff(names(fdat), c("crisisJST", "country", "year"))
  fdat <- fdat <- PredictingBlackSwans::lagFun(fdat, inputVars=lagVars, maxLag=2)

  # Selection
  varSel <- setdiff(names(fdat), c("crisisJST", "country", "year"))

  # Remove NAs
  fdat <- na.omit(fdat)
  # Remove Infinities
  datAux <- fdat[, -c("crisisJST", "country", "year")]
  fdat <- fdat[is.finite(rowSums(datAux)), ]

  # Train and test
  train <- fdat[year < beginTest]
  test <- fdat[beginTest <= year]

  # Discard country, year for train sample
  train <- train[, setdiff(names(train), c("country", "year")), with=F]

  # Keep only indicator when it is either 0 or 1 (No crisis - Pre-crisis)
  train <- train[crisisJST %in% c(0, 1)]
  train[, crisisJST:=as.factor(crisisJST)] # Convert the response to factor

  # Random forest
  strataSampling <- as.vector(ceiling(table(train$crisisJST) * 0.632)) # Stratified Sampling
  rfMod <- randomForest::randomForest(x=as.matrix(train[, varSel, with=F]),
                                      y=train$crisisJST,
                                      strata=train$crisisJST,
                                      sampsize=strataSampling)
  ### Model performance on test data ###
  twt <- test[crisisJST %in% c(0, 1)] # Test sample only with non-crisis and pre-crisis periods
  twt$preds <- predict(rfMod, newdata=twt, type="prob")[, 2] # Predictions
  rocTest <- pROC::roc(response=twt$crisisJST, predictor=twt$preds, direction="<") # ROC
  rocDTtestRf <- data.table(x=(1 - rocTest$specificities), y=rocTest$sensitivities)
  rocDTtestRf[, Model:="Random Forest"]
  outputTableTest[3, 2] <- round(as.numeric(rocTest$auc), 2)
  outputTableTest[3, 3] <- nrow(train)
  outputTableTest[3, 4] <- nrow(train[crisisJST==1])
  outputTableTest[3, 5] <- length(varSel)

  probsRf <- copy(test)
  probsRf$probs <- predict(rfMod, newdata=probsRf, type="prob")[, 2]
  probsRf <- probsRf[, .(country, year, crisisJST, probs)]
  probsRf[!(crisisJST %in% c(0, 1)), probs:=NA]

  #########################################################################################
  ### Lasso High-Dimensional ###
  #########################################################################################
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
  fdat <- PredictingBlackSwans::interactionsFun(inputData=fdat, inputVars=interVars)

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

  # Train and test
  train <- fdat[year < beginTest]
  test <- fdat[beginTest <= year]

  # Discard country, year for train sample
  train <- train[, setdiff(names(train), c("country", "year")), with=F]

  # Keep only indicator when it is either 0 or 1 (No crisis - Pre-crisis)
  train <- train[crisisJST %in% c(0, 1)]
  train[, crisisJST:=as.factor(crisisJST)] # Convert the response to factor
  train[, rowNumbers:=1:nrow(train)] # Get the row numbers for the MCCV
  if (parallel == T) doMC::registerDoMC(cores=ncores)
  ### Model performance on test data ###
  # Stratified cross-validation
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
  twt <- test[crisisJST %in% c(0, 1)] # Test sample only with non-crisis / pre-crisis observations
  twt$preds <- glmnet::predict.cv.glmnet(lassoMod, newx=as.matrix(twt[, varSel, with=F]),
                                         s="lambda.1se",
                                         type="response")
  rocTest <- pROC::roc(response=twt$crisisJST, predictor=twt$preds, direction="<") # ROC
  rocDTtestLasso <- data.table(x=(1 - rocTest$specificities), y=rocTest$sensitivities)
  rocDTtestLasso[, Model:="Lasso"]
  outputTableTest[2, 2] <- round(as.numeric(rocTest$auc), 2)
  outputTableTest[2, 3] <- nrow(train)
  outputTableTest[2, 4] <- nrow(train[crisisJST==1])
  outputTableTest[2, 5] <- length(varSel)

  probsLasso <- copy(test)
  probsLasso$probs <- as.vector(glmnet::predict.cv.glmnet(lassoMod,
                                                          newx=as.matrix(probsLasso[, varSel, with=F]),
                                                          s="lambda.1se", type="response"))
  probsLasso <- probsLasso[, .(country, year, crisisJST, probs)]
  probsLasso[!(crisisJST %in% c(0, 1)), probs:=NA]

  ############################################################################################
  ### ROC comparison ###
  ############################################################################################

  # Order the columns
  rocDTtestLogit <- rocDTtestLogit[order(x, y)]
  rocDTtestRf <- rocDTtestRf[order(x, y)]
  rocDTtestLasso <- rocDTtestLasso[order(x, y)]
  rocCompTest <- rbindlist(list(rocDTtestLogit, rocDTtestRf, rocDTtestLasso))
  rocCompTest[, ModOrd:=factor(Model, levels=c("LR", "Lasso", "Random Forest"))]

  pRocTest <- ggplot(data=rocCompTest, aes(x=x, y=y, color=ModOrd)) +
    geom_step() +
    geom_segment(x=0, y=0, xend=1, yend=1, linetype=2, color="black") +
    labs(x="False Positive Rate (FPR)", y="True Positive Rate (TPR)") +
    theme_bw() +
    theme(plot.title=element_text(hjust=.5)) +
    scale_x_continuous(labels=scales::percent) +
    scale_y_continuous(labels=scales::percent) +
    guides(color=guide_legend(title = "Model"))

  ################################################################################################
  ### Test Predictions ###
  ################################################################################################

  # LR
  p <- PredictingBlackSwans::plotCrisProb(crisDT=probsLogit[, .(country, year, probs)],
                                          fullDT=test,
                                          countryName="All")
  pLogitProbs <- p$plot + theme_bw() + theme(panel.grid=element_blank())

  # Lasso
  p <- PredictingBlackSwans::plotCrisProb(crisDT=probsLasso[, .(country, year, probs)],
                                          fullDT=test,
                                          countryName="All")
  pLassoProbs <- p$plot + theme_bw() + theme(panel.grid=element_blank())

  # Random Forest
  p <- PredictingBlackSwans::plotCrisProb(crisDT=probsRf[, .(country, year, probs)],
                                          fullDT=test,
                                          countryName="All")
  pRfProbs <- p$plot + theme_bw() + theme(panel.grid=element_blank())

  outputList <- list(outputTableTest=outputTableTest,
                     pRocTest=pRocTest,
                     pLogitProbs=pLogitProbs,
                     pLassoProbs=pLassoProbs,
                     pRfProbs=pRfProbs)
  return(outputList)
}
