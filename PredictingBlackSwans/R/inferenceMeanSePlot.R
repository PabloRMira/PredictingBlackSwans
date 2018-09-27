#' Plot the standard errors of the estimate against the estimates
#'
#' Plot the Monte-Carlo standard errors of the estimate (y-axis) against the
#' Monte-Carlo estimates (x-axis) to assess the drivers of the inference results
#' regarding the worse coverage with simultaneouly better power properties of the tests.
#' @param path Path to look for the input-files.
#' @export
#' @keywords PredictingBlackSwans
inferenceMeanSePlot <- function(path) {

  # Load the data
  setwd(path)

  # n = 100, p = 150, rho = 0.5
  list1 <- readRDS(file="Inference_Results_n_100_p_150_rho_0.5.rds")

  # n = 200, p = 250, rho = 0.5
  list2 <- readRDS(file="Inference_Results_n_200_p_250_rho_0.5.rds")

  # Extract the estimates for the non-zero coefficients for the DLassoCV
  posNz <- list1$n_100_p_150_rho_0.5$posNz
  df1 <- list1[["n_100_p_150_rho_0.5"]]$coefArray[posNz,,3]
  df2 <- list2[["n_200_p_250_rho_0.5"]]$coefArray[posNz,,3]

  # Put the frames in shape
  df1 <- as.data.table(t(df1))
  names(df1) <- paste0("Coef Number: ", posNz)
  mdf1 <- melt(df1, measure.vars = names(df1))
  mdf1[, Scenario:="$n = 100, p = 150$"]

  df2 <- as.data.table(t(df2))
  names(df2) <- paste0("Coef Number: ", posNz)
  mdf2 <- melt(df2, measure.vars = names(df2))
  mdf2[, Scenario:="$n = 200, p = 250$"]

  mdf <- rbindlist(list(mdf1, mdf2))
  names(mdf) <- c("coefNumber", "estimate", "scenario")

  # Now the standard errors
  df1 <- list1[["n_100_p_150_rho_0.5"]]$seArray[posNz,,3]
  df2 <- list2[["n_200_p_250_rho_0.5"]]$seArray[posNz,,3]

  # Put the frames in shape
  df1 <- as.data.table(t(df1))
  names(df1) <- paste0("Coef Number: ", posNz)
  mdf1 <- melt(df1, measure.vars = names(df1))
  mdf1[, Scenario:="$n = 100, p = 150$"]

  df2 <- as.data.table(t(df2))
  names(df2) <- paste0("Coef Number: ", posNz)
  mdf2 <- melt(df2, measure.vars = names(df2))
  mdf2[, Scenario:="$n = 200, p = 250$"]

  mdfse <- rbindlist(list(mdf1, mdf2))
  mdfse[, variable:=NULL]
  mdfse[, Scenario:=NULL]
  names(mdfse) <- "se"

  mdfComb <- cbind(mdf, mdfse)
  mdfComb[, coefNumber:=factor(coefNumber, levels=paste0("Coef Number: ", posNz),
                               labels=paste0("Coef Number: ", posNz))]

  # Add boundary lines for the tests

  # H_0 : \beta_j = 0

  nzCoef <- c(1, -0.5, 2, 1, -0.5, 2)

  null0 <- data.table(coefNumber=paste0("Coef Number: ", posNz),
                      slope=1/qnorm(0.975))
  null0two <- data.table(coefNumber=paste0("Coef Number: ", posNz),
                      slope=-1/qnorm(0.975))
  null0 <- rbindlist(list(null0, null0two))
  null0[, intercept:=0]
  null0[, coefNumber:=factor(coefNumber, levels=paste0("Coef Number: ", posNz),
                               labels=paste0("Coef Number: ", posNz))]
  null0[, boundary:="$H_0: \\beta^0 = 0$"]

  null1 <- data.table(coefNumber=paste0("Coef Number: ", posNz),
                      slope=1/qnorm(0.975),
                      intercept=-nzCoef/qnorm(0.975))
  null1two <- data.table(coefNumber=paste0("Coef Number: ", posNz),
                         slope=-1/qnorm(0.975),
                         intercept=nzCoef/qnorm(0.975))
  null1 <- rbindlist(list(null1, null1two))
  null1[, coefNumber:=factor(coefNumber, levels=paste0("Coef Number: ", posNz),
                             labels=paste0("Coef Number: ", posNz))]
  null1[, boundary:="$H_0: \\beta^0 = \\beta^0$"]

  nullComb <- rbindlist(list(null0, null1))

  tikzDevice::tikz(file="Inference_Plot_Explanation.tex", width=8, height=4.5)
  p <- ggplot(data=mdfComb, aes(x=estimate, y=se, color=scenario, shape=scenario)) +
    geom_point(alpha=.75) +
    geom_abline(data=nullComb, aes(intercept=intercept, slope=slope, linetype=boundary)) +
    facet_wrap(~ coefNumber, scales="free_x") +
    theme_bw() +
    labs(x="$\\hat{\\beta}$", y="se$(\\hat{\\beta})$") +
    guides(color=guide_legend(title="Scenario", title.position="top", title.hjust=.5,
                              order=1),
           shape=guide_legend(title="Scenario", title.position="top", title.hjust=.5,
                              order=1),
           linetype=guide_legend(title="Boundary", title.position = "top",
                                 title.hjust = .5)) +
    theme(legend.position = "bottom", panel.grid = element_blank()) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, NA))
  print(p)
  dev.off()
}
