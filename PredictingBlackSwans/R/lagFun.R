#' Compute lags of variables
#'
#' @param inputData A data.table.
#' @param inputVars Variables for which the lags will be computed.
#' @param maxLag Number up to which lags will be computed. E.g. if maxLag = 3
#' then the first, second and third lag will be calculated.
#' @return Additional variables in the data.table. These will have the names
#' of the input variables with an "_jL"-ending, with j = 1, ..., maxAve.
#' @export
#' @keywords PredictingBlackSwans
#' @examples lagFun(inputData=dat, inputVars=c("gdp", "revenue"), maxLag=3)
lagFun <- function(inputData, inputVars, maxLag) {
  for (i in 1:length(inputVars)) {
    inputData[, paste0(inputVars[i], "_", 1:maxLag, "L"):=shift(.SD, 1:maxLag), .SDcols=inputVars[i],
              by=country]
  }
  return(inputData)
}
