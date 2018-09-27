#' Auxiliary function to compute the standard deviation with factor (1 / n)
#'
#' @param y A vector.
#' @return The standard deviation of the vector y with factor (1 / n) instead
#' of the default in base R (1 / (n- 1)).
#' @export
mysd <- function(y) {
  sqrt(sum((y - mean(y))^2) / length(y))
}
