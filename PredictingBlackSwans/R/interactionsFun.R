#' Compute interactions terms
#'
#' @param inputData A data.table.
#' @param inputVars Variables for which the interactions will be computed.
#' @return Additional variables in the dataset. These have the names
#' of the pairs input variables with an "_i_" in between.
#' @export
#' @keywords PredictingBlackSwans
#' @examples interactionsFun(data=dat, inputVars=c("gdp", "revenue"))
interactionsFun <- function(inputData, inputVars) {
  data <- copy(inputData)
  for (i in 1:(length(inputVars)-1)) {
    for (j in (i+1):length(inputVars)) {
      data[, paste(inputVars[i], inputVars[j], sep="_i_"):=get(inputVars[i]) * get(inputVars[j])]
    }
  }
  return(data)
}
