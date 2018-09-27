#' Analyzing Black Swans
#'
#' Main function to reproduce the results of our paper regarding the
#' inference analysis.
#' @param path Path to export the results
#' @param cvfolds Number of folds for the cross-validation in the desparsified Lasso
#' for the initial estimator (Lasso).
#' @param seed Seed for reproducibility.
#' @param parallel Should parallel computing be used? (Only for UNIX computers)
#' @param ncores Number of cores for parallel-computing.
#' @details A new directory "Inference_Analysis" will be created
#' where the results will be placed in.
#' @export
#' @keywords PredictingBlackSwans
analyzingBS <- function(path=getwd(),
                        cvfolds=5,
                        seed=813,
                        parallel=T,
                        ncores=2L) {
  # Get current directory
  curDir <- getwd()

  # Seed for reproducibility
  set.seed(seed)

  # Create folder structure
  inferPath <- file.path(path, "Inference_Analysis")
  misPath <- file.path(inferPath, "Missing_Values_Analysis")
  misCountryPath <- file.path(misPath, "Country_Analysis")

  dir.create(inferPath, showWarnings = FALSE)
  dir.create(misPath, showWarnings = F)
  dir.create(misCountryPath, showWarnings = F)

  # Initialize regressin output
  regOutput <- data.frame(matrix(NA, 2, 8), stringsAsFactors = F)
  names(regOutput) <- c("Model Specification", "Estimate", "p-value", "Corrected p-value",
                        "", "Estimate", "p-value", "Corrected p-value")
  regOutput[, 1] <- c("L1. log(Real total credit)", "L2. log(Real total credit)")
  regOutput[, 5] <- c("", "") # Empty column

  # Set working directory for exporting the graphs
  setwd(inferPath)

  ##################################################################################
  ### Missing Values Analysis ###
  ##################################################################################

  PredictingBlackSwans::missingValuesAnalysis(misPath)

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

  # Take observations only from 1950 onwards because of being only interested
  # in the second era of financial capitalism (Schularick and Taylor (2012))
  # and because of broken series
  dat <- dat[year >= 1950]

  ##################################################################################
  ### Generate new features ###
  ##################################################################################
  set.seed(seed)
  fdat <- copy(dat)
  # First selection of raw variables
  varsFirst <- c("rconpc", "rgdpmad", "pop", "gdp", "imports", "exports", "ltrate", "stir", "ca",
                 "debtgdp", "revenue", "expenditure", "tloans", "hpnom", "narrowm", "cpi",
                 "iy", "stocks")
  fdat <- fdat[, c("crisisJST", "country", "year", varsFirst), with=F]

  # Impute the only missing values for Denmark for the variable "debtgdp" in 1997 using
  # a linear interpolation between the two 1996 and 1998 for this variable
  fdat[, aux1:=shift(debtgdp), by=country]
  fdat[, aux2:=shift(debtgdp, type="lead"), by=country]
  fdat[country=="Denmark" & year == 1997, debtgdp:=(aux1 + aux2) / 2]
  fdat[, (c("aux1", "aux2")):=NULL]

  # Generate new features
  fdat[, tloans_cpi:=tloans / cpi] # <- Variable of interest (with its lags)
  fdat[, netexport_gdp:=(exports - imports) / gdp] # Net export to GDP
  fdat[, govbalance_gdp:=(revenue - expenditure) / gdp] # Government balance to GDP
  fdat[, ca_gdp:=ca / gdp] # Current account to GDP

  ### Detrend variables using the Hamilton filter (2017) after logarihmizing ###
  # Note: We logarithmized to take into account exponential trends
  logVars <- setdiff(names(fdat), c("country", "year", "crisisJST", "ltrate", "stir", "iy",
                                    "netexport_gdp", "govbalance_gdp", "ca", "ca_gdp", "debtgdp",
                                    "tloans"))
  for (i in 1:length(logVars)) {
    fdat[, paste0("l", logVars[i]):=log(get(logVars[i]))]
    fdat <- PredictingBlackSwans::hamiltonFilter(fdat, paste0("l", logVars[i]), h=3)
  }

  # Hamilton filter without logarithmized variables
  myvars <- c("ltrate", "stir", "netexport_gdp", "govbalance_gdp", "ca_gdp", "iy", "debtgdp")
  for (i in 1:length(myvars)) {
    fdat <- PredictingBlackSwans::hamiltonFilter(fdat, myvars[i], h=3)
  }

  # Take only filtered variables
  varSel <- c("country", "year", "crisisJST", names(fdat)[grep(names(fdat), pattern="_dt")])
  fdat <- fdat[, varSel, with=F]

  # Take only complete time series and not broken ones (through missing values)
  countryVec <- unique(as.character(fdat$country))
  for (i in 1:length(countryVec)) {
    mycountry <- countryVec[i]
    auxDT <- fdat[country==mycountry]
    auxDT[, rowNumber:=1:nrow(auxDT)]
    auxDT2 <- auxDT[, -c("country", "year", "crisisJST"), with=F]
    auxDT2[, yearNA:=rowSums(auxDT2)]
    maxYear <- auxDT[rowNumber==auxDT2[is.na(yearNA), max(rowNumber)], year]
    fdat[country==mycountry & year <= maxYear, (names(fdat)):=NA]
  }
  fdat <- na.omit(fdat)

  # Dickey-Fuller and KPSS tests for stationarity for each country for each variable
  # Count only non-stationarity for the countries that both do not reject the
  # Dickey-Fuller test for unit root and reject the KPSS test for stationarity
  varSel <- setdiff(names(fdat), c("country", "year", "crisisJST"))
  nomSize <- 0.05
  nonStatioTable <- data.frame(Variable=varSel,
                            PercNoStationary=NA,
                            Countries="a", stringsAsFactors = F)
  for (j in 1:length(varSel)) {
    myvar <- varSel[j]
    noStationaryCheck <- numeric(length(countryVec))
    countryNonStatio <- as.character()
    for(i in 1:length(countryVec)) {
      mycountry <- countryVec[i]
      var <- fdat[country==mycountry, myvar, with=F][[1]]
      kpsstest <- tseries::kpss.test(var, null="Level")
      adftest <- tseries::adf.test(var, alternative="stationary")
      noStationaryCheck[i] <- ifelse(adftest$p.value > nomSize & kpsstest$p.value < nomSize, 1, 0)
      if (noStationaryCheck[i] == 1) countryNonStatio <- c(countryNonStatio, as.character(mycountry))
    }
    nonStatioTable[j, 2] <- mean(noStationaryCheck)
    nonStatioTable[j, 3] <- as.character(paste(countryNonStatio, collapse=", "))
  }
  names(nonStatioTable) <- c("Variable", "\\% Non-Stationary", "Countries for which non-stationarity arises")
  varNames <- c("log(Real consumption)", "log(Real GDP per capita)", "log(Population)", "log(GDP)", "log(Imports)",
                "log(Exports)", "log(Revenue)", "log(Expenditure)", "log(House Price Index)", "log(Narrow money)",
                "log(Consumer Price Index)", "log(Stock Price Index)", "log(Real total credit)", "LT interest rate",
                "ST interest rate", "Net export", "Gov. balance", "Current account / GDP", "Investment / GDP", "Gov. debt / GDP")
  nonStatioTable$Variable <- varNames

  setwd(inferPath)
  stargazer::stargazer(nonStatioTable, type="latex", out="Non_Stationarity_Table.tex", summary=F, rownames = F)

  # Take interactions
  interVars <- setdiff(names(fdat), c("crisisJST", "year", "country", "ltloans_cpi_dt"))
  fdat <- PredictingBlackSwans::interactionsFun(fdat, interVars)

  # Take 4 lags for all variables except our variable of interest, ltloans_cpi_dt
  varSel <- setdiff(names(fdat), c("country", "year", "crisisJST", "ltloans_cpi_dt"))
  fdat <- PredictingBlackSwans::lagFun(fdat, varSel, 4)

  # For our variable of interest take only one lag
  # Note: In the paper we claim that we run a regression of the crisis indicator
  # a lagged polynomial of regressors. Our equivalent implementation runs a regression of the
  # 1-year pre-crisis indicator on contemporaneous regression to it. This is why we only take
  # an additional lag of the ltloans_cpi_dt variable since this variable is already a lag
  # in the same sense as in the paper.
  varSel <- "ltloans_cpi_dt"
  fdat <- PredictingBlackSwans::lagFun(fdat, varSel, 1)

  # Remove NAs
  fdat <- na.omit(fdat)

  # Remove the crisis periods
  fdat <- fdat[crisisJST!=2]

  # Convert crisis indicator to factor
  fdat[, crisisJST:=as.factor(crisisJST)]

  logitMod <- glm(crisisJST ~ ltloans_cpi_dt + ltloans_cpi_dt_1L + country,
                  data=fdat,
                  family="binomial")
  summaryLogit <- summary(logitMod)
  regOutput[, 2] <- as.vector(round(summaryLogit$coefficients[,1][2:3], 4))
  pval <- as.vector(round(summaryLogit$coefficients[, 4][2:3], 4))
  pvalueStar <- ifelse(pval > 0.1, paste0(pval), ifelse(pval > 0.05, paste0(pval, "$^\\star$"),
                                                        ifelse(pval > 0.01, paste0(pval, "$^{\\star\\star}$"),
                                                               paste0(pval, "$^{\\star\\star\\star}$"))))
  regOutput[, 3] <- pvalueStar
  pvalcorr <- as.vector(p.adjust(summaryLogit$coefficients[, 4][2:3], method="holm"))
  pval <- round(pvalcorr, 4)
  pvaluecorrStar <- ifelse(pval > 0.1, paste0(pval), ifelse(pval > 0.05, paste0(pval, "$^\\star$"),
                                                            ifelse(pval > 0.01, paste0(pval, "$^{\\star\\star}$"),
                                                                   paste0(pval, "$^{\\star\\star\\star}$"))))
  regOutput[, 4] <- pvaluecorrStar

  # Inference: Desparsified Lasso CV
  x <- as.matrix(fdat[, -c("country", "year", "crisisJST"), with=F])
  y <- fdat$crisisJST
  # Initialize parallel computing
  if (parallel == T) doMC::registerDoMC(cores=ncores)

  # Indexes of the low-dimensional component of interest
  lowSel <- grep(pattern = "ltloans_cpi_dt", x=attributes(x)$dimnames[[2]])
  # Names of the low-dimensional component of interest
  lowNames <- attributes(x)$dimnames[[2]][lowSel]

  # Desparsified Lasso for a low dimensional component
  dLasso <- PredictingBlackSwans::despLassoLowComp(x=x, y=y, lowSel=lowSel, cvfolds=5, parallel=T, ncores=2L)
  regOutput[, 6] <- round(as.vector(dLasso$dLasso), 4)
  pval <- round(as.vector(dLasso$pval), 4)
  pvalueStar <- ifelse(pval > 0.1, paste0(pval), ifelse(pval > 0.05, paste0(pval, "$^\\star$"),
                                                        ifelse(pval > 0.01, paste0(pval, "$^{\\star\\star}$"),
                                                               paste0(pval, "$^{\\star\\star\\star}$"))))
  regOutput[, 7] <- pvalueStar
  pval <- round(as.vector(dLasso$pcorr), 4)
  pvaluecorrStar <- ifelse(pval > 0.1, paste0(pval), ifelse(pval > 0.05, paste0(pval, "$^\\star$"),
                                                               ifelse(pval > 0.01, paste0(pval, "$^{\\star\\star}$"),
                                                                      paste0(pval, "$^{\\star\\star\\star}$"))))
  regOutput[, 8] <- pvaluecorrStar

  # Table
  mynote <- paste("Dependent variable: crisis indicator. Treatment effect: Credit exuberance indicator",
                  "represented by two lags of the deviation from trend of log real total credit.",
                  "Corrected p-values using the Bonferroni-Holm procedure.",
                  "The low-dimensional logistic regression model uses country fixed effects.",
                  "The desparsified Lasso procedure uses the square-root Lasso for the nodewise regression.",
                  "The high-dimensional analysis includes 737 observations and 955 variables.")

  outTexTable <- xtable::xtable(regOutput, type="latex",
                                caption="Analysis of the imbalances prior to a financial crisis",
                                align="llcllp{0.1cm}cll")
  xtable::print.xtable(outTexTable, type="latex", file="DLasso_Application.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(regOutput), nrow(regOutput)),
                                                                     command=c(" \\toprule",
                                                                               " \\multicolumn{1}{l}{} &",
                                                                               " \\multicolumn{3}{c}{LR Model} &",
                                                                               " \\multicolumn{1}{c}{} &",
                                                                               " \\multicolumn{3}{c}{Desparsified Lasso} \\\\",
                                                                               " \\cmidrule(l r){2-4} \\cmidrule(l r){6-8}",
                                                                               " \\bottomrule \\\\[-.3cm]",
                                                                               paste0(" \\multicolumn{8}{l}{\\parbox{20cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))
  # Return back to original working directory
  setwd(curDir)
}
