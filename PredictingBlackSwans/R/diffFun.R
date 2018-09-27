#' Compute time differences
#'
#' Compute time differences y_t - t_{t-1}.
#' @param inputData A data.table.
#' @param inputVars Variables for which the time differences will be calculated.
#' @return Additional variables in the dataset. These have the names
#' of the input variables with a "_diff"-ending.
#' @export
#' @keywords PredictingBlackSwans
#' @examples diffFun(data=dat, inputVars=c("gdp", "revenue"))
diffFun <- function(inputData, inputVars) {
  data <- copy(inputData)
  for (i in 1:length(inputVars)) {
    data[, auxvar:=shift(get(inputVars[i])), by=country]
    data[, paste0(inputVars[i], "_diff"):=get(inputVars[i]) - auxvar,
         by=country]
  }
  data[, auxvar:=NULL]
  return(data)
}
