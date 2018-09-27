#' Implementation of the Hamilton (2017) filter
#'
#' Compute the deviations from trend for a time series using the Hamilton filter for each
#' country separately.
#' @param inputData A data.table.
#' @param inputVar The variable for detrending.
#' @param h Horizon for which we build a prediction.
#' @return A new data.table with an additional variable. This variable has the name
#' of the input variable with a '_dt'-ending.
#' @export
#' @keywords PredictingBlackSwans
hamiltonFilter <- function(inputData,
                           inputVar,
                           h=3) {
  inputAux <- copy(inputData)
  outVar <- paste0(inputVar, "_dt")
  inputAux[, paste0("auxvar", seq(1,4)):=shift(.SD, (seq(0, 3)+h)), .SDcols=inputVar,
           by=country]
  dataAux <- inputAux[, c("country", "year", inputVar, paste0("auxvar", seq(1,4))), with=F]
  dataAux <- na.omit(dataAux)
  countryVec <- as.character(unique(inputAux$country))
  dataList <- list()
  fmla <- as.formula(paste0(inputVar, " ~ auxvar1 + auxvar2 + auxvar3 + auxvar4"))
  for (j in 1:length(countryVec)) {
    countryData <- dataAux[country==countryVec[j]]
    countryData[, (outVar):=lm(fmla, data=countryData)$residuals]
    dataList[[j]] <- countryData[, c("country", "year", outVar), with=F]
  }
  dataComb <- rbindlist(dataList)
  outputData <- merge(inputData, dataComb, by=c("country", "year"), all.x=T)
  return(outputData)
}
