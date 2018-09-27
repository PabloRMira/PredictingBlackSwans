#' Compute growth rates
#'
#' Compute growth rates as (y_t - y_{t-1}) / y_{t-1}.
#' @param data A data.table.
#' @param inputVars Variables for which the growth rates will be calculated.
#' @return Additional variables in the data.table. These have the names
#' of the input variables with a "_gr"-ending.
#' @export
#' @keywords PredictingBlackSwans
#' @examples growthFun(data=dat, inputVars=c("gdp", "revenue"))
growthFun <- function(data, inputVars) {
  for (i in 1:length(inputVars)) {
    data[, auxvar:=shift(get(inputVars[i])), by=country]
    data[, paste0(inputVars[i], "_gr"):=(get(inputVars[i]) - auxvar) / auxvar,
         by=country]
  }
  data[, auxvar:=NULL]
}
