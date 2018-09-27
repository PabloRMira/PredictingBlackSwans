#' Compute log-differences
#'
#' Compute log-differences \eqn{log(y_t) - log(y_{t-1})}
#' @param inputData A data.table.
#' @param inputVars Variables for which the log-differences will be calculated.
#' @return Additional variables in the data.table. These have the names
#' of the input variables with a "_lDiff"-ending.
#' @export
#' @keywords PredictingBlackSwans
#' @examples logDiff(data=dat, inputVars=c("gdp", "revenue"))
logDiffFun <- function(inputData, inputVars) {
  data <- copy(inputData)
  for (i in 1:length(inputVars)) {
    data[, auxvar:=shift(get(inputVars[i])), by=country]
    data[, paste0(inputVars[i], "_lDiff"):=(log(get(inputVars[i])) - log(auxvar)),
         by=country]
  }
  data[, auxvar:=NULL]
  return(data)
}
