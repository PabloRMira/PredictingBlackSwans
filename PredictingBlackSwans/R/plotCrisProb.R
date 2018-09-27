#' Plot Crisis Probabilities for each Country
#'
#' Plot crisis probabilities for each country in the sample.
#' @param crisDT A data table with the countries, years
#' and predicted probability.
#' @param fullDT The full data set.
#' @param countryName Name of the country to plot. It can also be "All".
#' @return A ggplot.
#' @export
#' @keywords PredictingBlackSwans
plotCrisProb <- function(crisDT,
                         fullDT,
                         countryName) {

  # Rename the full data set
  dat <- fullDT

  # If countryName = All then do a facet wrap
  if (countryName == "All") {

    # Preprocess data set
    varDat <- dat[, c("country" , "year", "crisisJST"), with=FALSE]
    varDat <- na.omit(varDat)
    varDat <- merge(varDat,
                    crisDT,
                    sort=FALSE,
                    by=c("country", "year"),
                    all.x=T)
    minYear <- min(crisDT[, year])
    maxYear <- max(crisDT[, year])

    # Data set
    df <- varDat[year>=minYear & year<=maxYear]
    names(df) <- c("country", "year", "crisis", "probs")

    # Get crisis and precrisis dates
    dfDates <- df[, c("year", "country", "crisis")]

    # Precrisis dates
    dfPrecrisis <- dfDates[crisis==1]
    if (nrow(dfPrecrisis) > 0) {
      dfPrecrisis[, auxvar1:=shift(year, n=1, type="lag"), by="country"]
      dfPrecrisis[, auxvar2:=shift(year, n=1, type="lead"), by="country"]
      dfPrecrisis[, ind1:=ifelse(year - auxvar1 > 1, 1, 0)]
      dfPrecrisis[is.na(ind1), ind1:=1]
      preDates <- dfPrecrisis[ind1==1, c("year", "country")]
      names(preDates)[1] <- "start"

      dfPrecrisis[, ind2:=ifelse(auxvar2 - year > 1, 1, 0)]
      dfPrecrisis[is.na(ind2), ind2:=1]
      preDates[, end:=dfPrecrisis[ind2==1, year]]
      preDates[, end:=end+1]
      preDates[end > maxYear, end:=maxYear]
      # Modify the end for the
      # shadowed area in the plot to look nicer
      preDates$end <- preDates$end - .5
    }

    # Crisis dates
    dfCrisis <- dfDates[crisis==2]
    if (nrow(dfCrisis) > 0) {
      dfCrisis[, auxvar1:=shift(year, n=1, type="lag"), by="country"]
      dfCrisis[, auxvar2:=shift(year, n=1, type="lead"), by="country"]
      dfCrisis[, ind1:=ifelse(year - auxvar1 > 1, 1, 0)]
      dfCrisis[is.na(ind1), ind1:=1]
      criDates <- dfCrisis[ind1==1, c("year", "country")]
      names(criDates)[1] <- "start"

      dfCrisis[, ind2:=ifelse(auxvar2 - year > 1, 1, 0)]
      dfCrisis[is.na(ind2), ind2:=1]
      criDates[, end:=dfCrisis[ind2==1, year]]
      criDates[, end:=end+1]
      criDates[end > maxYear, end:=maxYear]
      # Modify the end and start for the
      # shadowed area in the plot to look nicer
      criDates$start <- criDates$start - .5
      criDates$end <- criDates$end - .5
    }

    # Plot with all countries and corresponding variable
    p <- ggplot() +
      geom_line(data=df, aes(x=year, y=probs, group=country)) +
      geom_point(data=df, aes(x=year, y=probs, group=country))
    if (nrow(dfPrecrisis) > 0) {
      p <- p + geom_rect(data=preDates,
                         aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, group=country),
                         inherit.aes=FALSE, alpha=0.3)
    }
    if (nrow(dfCrisis) > 0) {
      p <- p + geom_rect(data=criDates,
                         aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, group=country),
                         inherit.aes=FALSE, alpha=0.8)
    }
    p <- p + labs(x="Time",
                  y="Crisis Probabilities") +
      theme(plot.title=element_text(hjust=.5)) +
      facet_wrap(~ country)
  } else {
    # Preprocess data set
    countryVarDat <- dat[country==countryName,
                         c("country", "year", "crisisJST"),
                         with=FALSE]
    countryVarDat <- na.omit(countryVarDat)
    countryVarDat <- merge(countryVarDat,
                           crisDT,
                           sort=FALSE,
                           by=c("country", "year"),
                           all.x=T)
    minYear <- min(crisDT[country==countryName, year])
    maxYear <- max(crisDT[country==countryName, year])

    # Dataset
    df <- countryVarDat[year >= minYear & year <= maxYear]
    names(df) <- c("country", "year", "crisis", "probs")

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
      # Modify the end for the
      # shadowed area in the plot to look nicer
      preDates$end <- preDates$end - .5
      preDates$start <- preDates$start - .5
    } else {
      preDates <- NULL
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
      # Modify the end and start for the
      # shadowed area in the plot to look nicer
      criDates$start <- criDates$start - .5
      criDates$end <- criDates$end - .5
    } else {
      criDates <- NULL
    }

    # Plot with country and variable
    p <- ggplot() +
      geom_line(data=df, aes(x=year, y=probs)) +
      geom_point(data=df, aes(x=year, y=probs))
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
                  y="Crisis Probabilities") +
      theme_classic() +
      theme(plot.title=element_text(hjust=.5))
  }
  outputList <- list(plot=p,
                     probsDf=df,
                     preDates=preDates,
                     criDates=criDates)
  return(outputList)
}
