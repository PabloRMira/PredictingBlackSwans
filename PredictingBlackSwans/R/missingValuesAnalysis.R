#' Missing Values Analysis
#'
#' Output missing values information.
#' @param path Path to export the results.
#' @export
#' @keywords PredictingBlackSwans
missingValuesAnalysis <- function(path) {
  setwd(path)
  ### All data ###
  dat <- as.data.table(PredictingBlackSwans::JST_rawData)
  datMis <- dat
  # Code factors as character
  factorVars <- names(sapply(datMis, class)[sapply(datMis, class) == "factor"])
  datMis[, (factorVars):=lapply(.SD, as.character), .SDcols=factorVars]
  datMis[!is.na(datMis)] <- 1
  datMis[is.na(datMis)] <- 0
  datMis[, (factorVars):=lapply(.SD, as.numeric), .SDcols=factorVars]
  datMis[, RowNumber:=seq(1, nrow(datMis))]
  mdt <- melt(datMis, id.vars="RowNumber")
  mdt$valueOrd <- factor(mdt$value, levels=c(0, 1), labels=c("Not Available", "Available"))

  p <- ggplot(data=mdt, aes(x=variable, y=RowNumber, fill=valueOrd)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.title=element_blank()) +
    scale_fill_manual(values=c("darkblue", "white"))
  ggsave("Missing_Values_All_Data_Overview.png", plot=p, units="cm", height=19, width=29)

  misDat <- colSums(datMis[, setdiff(names(datMis), "RowNumber"), with=F])
  mdt <- melt(misDat)
  names(mdt) <- "NotMissing"
  mdt$Variable <- row.names(mdt)
  mdt$Total <- nrow(datMis)
  mdt <- as.data.table(mdt)
  mdt[, Missing:=Total - NotMissing]
  mdt[, PercMissing:=Missing / Total]
  mdt[, PercNotMissing:=NotMissing / Total]
  mdt <- mdt[, .(Variable, Total, NotMissing, Missing, PercNotMissing, PercMissing)]

  # Output as html
  stargazer::stargazer(as.data.frame(mdt),type = "html",title = "Missing Values Overview",summary = F,
            out="Missing_Values_All_Data_Overview.html",digits = 3)

  ### Differentiate for 1s and 0s in the response ###
  # Load the raw data
  dat <- as.data.table(PredictingBlackSwans::JST_rawData)
  datMis <- dat
  # Code factors as character
  factorVars <- names(sapply(datMis, class)[sapply(datMis, class) == "factor"])
  datMis[, (factorVars):=lapply(.SD, as.character), .SDcols=factorVars]
  crisisJST <- datMis[, crisisJST]
  datMis[!is.na(datMis)] <- 1
  datMis[is.na(datMis)] <- 0
  datMis[, (factorVars):=lapply(.SD, as.numeric), .SDcols=factorVars]
  datMis$crisisJST <- crisisJST
  datMis0 <- datMis[crisisJST == 0]
  misDat0 <- colSums(datMis0[, setdiff(names(datMis0), c("RowNumber", "crisisJST")), with=F])
  mdt0 <- melt(misDat0)
  names(mdt0) <- "NotMissing"
  mdt0$Variable <- row.names(mdt0)
  mdt0$Total <- nrow(datMis0)
  mdt0 <- as.data.table(mdt0)
  mdt0[, Missing:=Total - NotMissing]
  mdt0[, PercMissing:=Missing / Total]
  mdt0[, PercNotMissing:=NotMissing / Total]
  mdt0 <- mdt0[, .(Variable, Total, NotMissing, Missing, PercNotMissing, PercMissing)]

  # Output
  stargazer::stargazer(as.data.frame(mdt0),type = "html",title = "Missing Values Overview (Non Crisis observations)",summary = F,
            out="Missing_Values_NonCrisis_Data_Overview.html",digits = 3)

  # 1's
  datMis1 <- datMis[crisisJST == 1]
  misDat1 <- colSums(datMis1[, setdiff(names(datMis1), c("RowNumber", "crisisJST")), with=F])
  mdt1 <- melt(misDat1)
  names(mdt1) <- "NotMissing"
  mdt1$Variable <- row.names(mdt1)
  mdt1$Total <- nrow(datMis1)
  mdt1 <- as.data.table(mdt1)
  mdt1[, Missing:=Total - NotMissing]
  mdt1[, PercMissing:=Missing / Total]
  mdt1[, PercNotMissing:=NotMissing / Total]
  mdt1 <- mdt1[, .(Variable, Total, NotMissing, Missing, PercNotMissing, PercMissing)]

  # Output
  stargazer::stargazer(as.data.frame(mdt1),type = "html",title = "Missing Values Overview (Crisis Observations)",summary = F,
            out="Missing_Values_Crisis_Data_Overview.html",digits = 3)

  # Now per country in time for all variables
  setwd("./Country_Analysis")
  dat <- as.data.table(PredictingBlackSwans::JST_rawData)
  countryVec <- unique(as.character(dat$country))
  datMis <- dat
  # Code factors as character
  factorVars <- names(sapply(datMis, class)[sapply(datMis, class) == "factor"])
  datMis[, (factorVars):=lapply(.SD, as.character), .SDcols=factorVars]
  crisisJST <- datMis[, crisisJST]
  year <- datMis[, year]
  country <- datMis[, country]
  datMis[!is.na(datMis)] <- 1
  datMis[is.na(datMis)] <- 0
  datMis[, (factorVars):=lapply(.SD, as.numeric), .SDcols=factorVars]
  datMis$crisisJST <- crisisJST
  datMis$year <- year
  datMis$country <- country
  # Remove unimportant variables
  datMis$iso <- NULL
  datMis$ifs <- NULL
  # Change the indicators: 0: No crisis / Not available; 1: No Crisis / Available; 2: Crisis / Not available; 3: Crisis / Available
  varVec <- setdiff(names(datMis), c("year", "country", "crisisJST"))
  for (col in varVec) {
    set(datMis, i=which(datMis[[col]]==0 & datMis[["crisisJST"]]==0), j=col, value=0)
    set(datMis, i=which(datMis[[col]]==1 & datMis[["crisisJST"]]==0), j=col, value=1)
    set(datMis, i=which(datMis[[col]]==0 & datMis[["crisisJST"]]==1), j=col, value=2)
    set(datMis, i=which(datMis[[col]]==1 & datMis[["crisisJST"]]==1), j=col, value=3)
  }
  datMis[, (varVec):=lapply(.SD, FUN=function(x) factor(x, levels=c(0,1,2,3),
                                                        labels=c("Non crisis & NA",
                                                                 "Non crisis & Available",
                                                                 "Crisis & NA",
                                                                 "Crisis & Available"))), .SDcols=varVec]
  datMis$crisisJST <- NULL # Remove response
  # For each country
  for (i in 1:length(countryVec)) {
    aCountry <- countryVec[i]
    auxDT <- datMis[country==aCountry, setdiff(names(datMis), "country"), with=F]
    auxMDT <- melt(auxDT, id.vars="year")
    p <- ggplot(data=auxMDT, aes(x=year, y=variable, fill=value)) +
      geom_tile(color="black") +
      labs(x="Year", y="Variable",
           title=paste("Missing Values - Country Overview:", aCountry)) +
      theme_bw() +
      theme(plot.title=element_text(hjust=.5),
            legend.title=element_blank())
    ggsave(paste0("Missing_Values_Country_", aCountry, ".png"), plot=p, units="cm", height=19, width=29)
  }
}
