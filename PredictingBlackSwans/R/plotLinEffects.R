#' Plot Linear Effects driving Crisis Probabilities
#'
#' Plot of the linear effects driving the crisis probabilities.
#' @param crisDT A data table with the countries, years
#' and predicted probability.
#' @param fullDT The full data set.
#' @param countryName Name of the country to plot. It can also be "All".
#' @return A ggplot.
#' @export
#' @keywords PredictingBlackSwans
plotLinEffects <- function(driverTable,
                           fullDT,
                           countryName) {

  # Rename full data set
  dat <- fullDT

  # Preprocess data sets
  countryVarDat <- dat[country==countryName,
                       c("country", "year", "crisisJST"),
                       with=FALSE]
  countryVarDat <- na.omit(countryVarDat)
  countryVarDat <- merge(countryVarDat,
                         driverTable,
                         sort=FALSE,
                         by=c("country", "year"),
                         all.x=TRUE)
  minYear <- min(driverTable[country==countryName, year])
  maxYear <- max(driverTable[country==countryName, year])

  # Dataset
  df <- countryVarDat[year >= minYear & year <= maxYear]
  names(df)[1:3] <- c("country", "year", "crisis")

  # Get crisis and precrisis dates
  dfDates <- df[, c("year", "crisis")]

  # Precrisis dates
  dfPrecrisis <- dfDates[crisis==1]
  if (nrow(dfPrecrisis) > 0) {
    dfPrecrisis[, auxvar1:=shift(year, n=1, type="lag")]
    dfPrecrisis[, auxvar2:=shift(year, n=1, type="lead")]
    dfPrecrisis[, ind1:=ifelse(year - auxvar1 > 1, 1, 0)]
    dfPrecrisis[1, ind1:=1]
    preDates <- dfPrecrisis[ind1==1, .(year)]
    names(preDates) <- "start"

    dfPrecrisis[, ind2:=ifelse(auxvar2 - year > 1, 1, 0)]
    dfPrecrisis[.N, ind2:=1]
    preDates[, end:=dfPrecrisis[ind2==1, year]]
    preDates[, end:=end+1]
    preDates[end > maxYear, end:=maxYear]
    preDates <- preDates - 0.5
  }

  # Crisis dates
  dfCrisis <- dfDates[crisis==2]
  if (nrow(dfCrisis) > 0) {
    dfCrisis[, auxvar1:=shift(year, n=1, type="lag")]
    dfCrisis[, auxvar2:=shift(year, n=1, type="lead")]
    dfCrisis[, ind1:=ifelse(year - auxvar1 > 1, 1, 0)]
    dfCrisis[1, ind1:=1]
    criDates <- dfCrisis[ind1==1, .(year)]
    names(criDates) <- "start"

    dfCrisis[, ind2:=ifelse(auxvar2 - year > 1, 1, 0)]
    dfCrisis[.N, ind2:=1]
    criDates[, end:=dfCrisis[ind2==1, year]]
    criDates <- criDates - 0.5
  }
  dfSel <- df[, -c("country", "crisis"), with=FALSE]
  mVars <- names(dfSel)[-1]
  dfSel[, netEffect:=rowSums(.SD), .SDcols=mVars]
  dfNE <- dfSel[, .(year, netEffect)]
  dfSel[, netEffect:=NULL]
  # Separate positive and negative effects for the area plot to work correctly
  varSel <- names(dfSel)[-1]
  meltDf <- melt(dfSel, id.vars="year", measure.vars=mVars)
  meltDf$positive <- ifelse(meltDf$value >= 0, meltDf$value, 0)
  meltDf$negative <- ifelse(meltDf$value < 0, meltDf$value, -1e-36) # <- This is needed to avoid
  # holes in the graph (see also https://stackoverflow.com/questions/51656490/holes-with-geom-area-ggplot2-for-negative-y-range)

  # Plot with country and variable
  p <- ggplot(data=meltDf) +
    geom_area(aes(x=year, y=positive, fill=variable), position="stack") +
    geom_area(aes(x=year, y=negative, fill=variable), position="stack") +
    geom_line(data=dfNE, aes(x=year, y=netEffect, color="Net Effect")) +
    geom_point(data=dfNE, aes(x=year, y=netEffect, color="Net Effect")) +
    guides(fill=guide_legend(title="Linear Effects", title.position = "top"),
           color=guide_legend(title="")) +
    scale_color_manual(values="black")
  if (nrow(dfPrecrisis) > 0) {
    p <- p + geom_rect(data=preDates,
                       aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
                       inherit.aes=FALSE, alpha=0.3)
  }
  if (nrow(dfCrisis) > 0) {
    p <- p + geom_rect(data=criDates,
                       aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
                       inherit.aes=FALSE, alpha=0.8)
  }
  p <- p + labs(title=countryName,
                x="Time",
                y="Linear Effects") +
    theme_classic() +
    theme(plot.title=element_text(hjust=.5),
          legend.position = "bottom", legend.title = element_text(hjust=.5))
  return(p)
}
