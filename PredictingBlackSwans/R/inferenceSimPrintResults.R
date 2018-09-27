#' Print the inference simulation results in Latex format
#'
#' Print the results of our simulation study for inference in Latex format as in the
#' master thesis.
#' @param inPath Path to look for the input-files.
#' @param outPath Path to export the results.
#' @export
#' @keywords PredictingBlackSwans
inferenceSimPrintResults <- function(inPath,
                                     outPath
) {

  # Preliminaries
  emptyRow <- data.frame(OracleLR="", InfeasibleLR="", DLassoCV="", DLassoSqrtBCW="",
                         DLassoSqrtVdG="")

  #################################################################################################
  #### Print results ###
  #################################################################################################

  # 1. n = 100, p = 150; Coverage

  # Load the data
  setwd(inPath)

  # 1.(a) n = 100, p = 150, rho = 0.5
  list1 <- readRDS(file="Inference_Results_n_100_p_150_rho_0.5.rds")

  # 1.(b) n = 100, p = 150, rho = 0.9
  list2 <- readRDS(file="Inference_Results_n_100_p_150_rho_0.9.rds")

  # Extract data frames for the coverage
  df1 <- list1[["n_100_p_150_rho_0.5"]]$mCoverage
  df2 <- list2[["n_100_p_150_rho_0.9"]]$mCoverage

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$",
                     "", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$")

  # Write the note for the coverage
  mynote <- paste("Average coverage for the active and non-active set,",
                  "for nominal coverage equal to 0.95.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Average coverage. n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="Coverage_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                         command=c(" \\toprule",
                                                                   " \\multicolumn{1}{l}{} &",
                                                                   " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                   " \\multicolumn{1}{c}{} &",
                                                                   " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                   " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                   " \\bottomrule \\\\[-.3cm]",
                                                                   paste0(" \\multicolumn{6}{l}{\\parbox{12cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 2. n = 100, p = 150; Tests, high-dimensional

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_150_rho_0.5"]]$tResults
  df2 <- list2[["n_100_p_150_rho_0.9"]]$tResults

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the high-dimensional coefficient vector. \\\\ n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestHD_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                         command=c(" \\toprule",
                                                                   " \\multicolumn{1}{l}{} &",
                                                                   " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                   " \\multicolumn{1}{c}{} &",
                                                                   " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                   " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                   " \\bottomrule \\\\[-.3cm]",
                                                                   paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 3. n = 100, p = 150; Tests, low-dimensional component

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_150_rho_0.5"]]$tResultsLow
  df2 <- list2[["n_100_p_150_rho_0.9"]]$tResultsLow

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the low-dimensional component. \\\\ n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestLD_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  #################################################################################################################################

  # 1. n = 200, p = 250; Coverage

  # Load the data
  setwd(inPath)

  # 1.(a) n = 200, p = 250, rho = 0.5
  list1 <- readRDS(file="Inference_Results_n_200_p_250_rho_0.5.rds")

  # 1.(b) n = 100, p = 150, rho = 0.9
  list2 <- readRDS(file="Inference_Results_n_200_p_250_rho_0.9.rds")

  # Extract data frames for the coverage
  df1 <- list1[["n_200_p_250_rho_0.5"]]$mCoverage
  df2 <- list2[["n_200_p_250_rho_0.9"]]$mCoverage

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$",
                     "", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$")

  # Write the note for the coverage
  mynote <- paste("Average coverage for the active and non-active set,",
                  "for nominal coverage equal to 0.95.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Average coverage. n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="Coverage_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                                     command=c(" \\toprule",
                                                                               " \\multicolumn{1}{l}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                               " \\multicolumn{1}{c}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                               " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                               " \\bottomrule \\\\[-.3cm]",
                                                                               paste0(" \\multicolumn{6}{l}{\\parbox{12cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 2. n = 200, p = 250; Tests, high-dimensional

  # Extract data frames for the tests
  df1 <- list1[["n_200_p_250_rho_0.5"]]$tResults
  df2 <- list2[["n_200_p_250_rho_0.9"]]$tResults

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the high-dimensional coefficient vector. \\\\ n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestHD_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 3. n = 200, p = 250; Tests, low-dimensional component

  # Extract data frames for the tests
  df1 <- list1[["n_200_p_250_rho_0.5"]]$tResultsLow
  df2 <- list2[["n_200_p_250_rho_0.9"]]$tResultsLow

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the low-dimensional component. \\\\ n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestLD_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

#############################################################################################

  # 1. n = 100, p = 500; Coverage

  # Load the data
  setwd(inPath)

  # 1.(a) n = 100, p = 500, rho = 0.5
  list1 <- readRDS(file="Inference_Results_n_100_p_500_rho_0.5.rds")

  # 1.(b) n = 100, p = 500, rho = 0.9
  list2 <- readRDS(file="Inference_Results_n_100_p_500_rho_0.9.rds")

  # Extract data frames for the coverage
  df1 <- list1[["n_100_p_500_rho_0.5"]]$mCoverage
  df2 <- list2[["n_100_p_500_rho_0.9"]]$mCoverage

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$",
                     "", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$")

  # Write the note for the coverage
  mynote <- paste("Average coverage for the active and non-active set,",
                  "for nominal coverage equal to 0.95.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Average coverage. n = 100, p = 500",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="Coverage_n_100_p_500.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                                     command=c(" \\toprule",
                                                                               " \\multicolumn{1}{l}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                               " \\multicolumn{1}{c}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                               " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                               " \\bottomrule \\\\[-.3cm]",
                                                                               paste0(" \\multicolumn{6}{l}{\\parbox{12cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 2. n = 100, p = 500; Tests, high-dimensional

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_500_rho_0.5"]]$tResults
  df2 <- list2[["n_100_p_500_rho_0.9"]]$tResults

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the high-dimensional coefficient vector. \\\\ n = 100, p = 500",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestHD_n_100_p_500.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 3. n = 100, p = 500; Tests, low-dimensional component

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_500_rho_0.5"]]$tResultsLow
  df2 <- list2[["n_100_p_500_rho_0.9"]]$tResultsLow

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("Family-wise error rate (FWER) and power of multiple testing,",
                  "for nominal FWER equal to 0.05.",
                  "OracleLR denotes the logistic regression estimator",
                  "that only uses the low-dimensional active set;",
                  "InfeasibleLR is the logistic regression estimator",
                  "that uses the active set and 6 additional noise variables;",
                  "DLassoCV denotes the desparsified Lasso estimator using",
                  "the nodewise Lasso with 5-fold cross-validation to",
                  "approximate the inverse;",
                  "DLassoSqrtBCW is the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "Belloni, Chernozhukov and Wang;",
                  "DLassoSqrtVdG denotes the desparsified Lasso with the",
                  "nodewise square-root Lasso using the rule proposed by",
                  "van de Geer.")

  # Print it to latex
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the low-dimensional component. \\\\ n = 100, p = 500",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestLD_n_100_p_500.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  #############################################################################################

  ### Additional simulations with "desparsified" outer product as in the paper ###
  # 1. n = 100, p = 150; Coverage

  # Load the data
  setwd(inPath)

  # 1.(a) n = 100, p = 150, rho = 0.5
  list1 <- readRDS(file="Inference_Results_Paper_n_100_p_150_rho_0.5.rds")

  # 1.(b) n = 100, p = 150, rho = 0.9
  list2 <- readRDS(file="Inference_Results_Paper_n_100_p_150_rho_0.9.rds")

  # Extract data frames for the coverage
  df1 <- list1[["n_100_p_150_rho_0.5"]]$mCoverage
  df2 <- list2[["n_100_p_150_rho_0.9"]]$mCoverage

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$",
                     "", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$")

  # Write the note for the coverage
  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Average coverage. n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="Coverage_Paper_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                                     command=c(" \\toprule",
                                                                               " \\multicolumn{1}{l}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                               " \\multicolumn{1}{c}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                               " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                               " \\bottomrule \\\\[-.3cm]",
                                                                               paste0(" \\multicolumn{6}{l}{\\parbox{12cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 2. n = 100, p = 150; Tests, high-dimensional

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_150_rho_0.5"]]$tResults
  df2 <- list2[["n_100_p_150_rho_0.9"]]$tResults

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the high-dimensional coefficient vector. \\\\ n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestHD_Paper_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 3. n = 100, p = 150; Tests, low-dimensional component

  # Extract data frames for the tests
  df1 <- list1[["n_100_p_150_rho_0.5"]]$tResultsLow
  df2 <- list2[["n_100_p_150_rho_0.9"]]$tResultsLow

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the low-dimensional component. \\\\ n = 100, p = 150",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestLD_Paper_n_100_p_150.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  #################################################################################################################################

  # 1. n = 200, p = 250; Coverage

  # Load the data
  setwd(inPath)

  # 1.(a) n = 200, p = 250, rho = 0.5
  list1 <- readRDS(file="Inference_Results_Paper_n_200_p_250_rho_0.5.rds")

  # 1.(b) n = 100, p = 150, rho = 0.9
  list2 <- readRDS(file="Inference_Results_Paper_n_200_p_250_rho_0.9.rds")

  # Extract data frames for the coverage
  df1 <- list1[["n_200_p_250_rho_0.5"]]$mCoverage
  df2 <- list2[["n_200_p_250_rho_0.9"]]$mCoverage

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$",
                     "", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}}^0$", "Avgcov $\\boldsymbol{\\beta}_{\\mathcal{S}^c}^0$")

  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  setwd(outPath)
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Average coverage. n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="Coverage_Paper_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames = F, replace=T,
                       booktabs = T, hline.after=0,  add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                                     command=c(" \\toprule",
                                                                               " \\multicolumn{1}{l}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                               " \\multicolumn{1}{c}{} &",
                                                                               " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                               " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                               " \\bottomrule \\\\[-.3cm]",
                                                                               paste0(" \\multicolumn{6}{l}{\\parbox{12cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 2. n = 200, p = 250; Tests, high-dimensional

  # Extract data frames for the tests
  df1 <- list1[["n_200_p_250_rho_0.5"]]$tResults
  df2 <- list2[["n_200_p_250_rho_0.9"]]$tResults

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  setwd(outPath)
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the high-dimensional coefficient vector. \\\\ n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestHD_Paper_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # 3. n = 200, p = 250; Tests, low-dimensional component

  # Extract data frames for the tests
  df1 <- list1[["n_200_p_250_rho_0.5"]]$tResultsLow
  df2 <- list2[["n_200_p_250_rho_0.9"]]$tResultsLow

  # Combine both results
  dfComb <- rbind(df1, emptyRow, df2)

  # Prepare it for the latex export
  dfComb <- as.data.frame(t(dfComb))
  dfComb$Methods <- row.names(dfComb)
  dfComb <- dfComb[c(6, seq(1, 5))]
  names(dfComb) <- c("Method", "FWER", "Power", "", "FWER", "Power")

  mynote <- paste("All specification as in Table X but now using the 'desparsifid outer product'")

  # Print it to latex
  outTexTable <- xtable::xtable(dfComb, type="latex", caption="Testing for the low-dimensional component. \\\\ n = 200, p = 250",
                                align="llccp{0.1cm}cc")
  xtable::print.xtable(outTexTable, type="latex", file="TestLD_Paper_n_200_p_250.tex",
                       caption.placement = "top", sanitize.text.function = function(x) {x}, include.rownames=F, replace=T, booktabs=T,
                       hline.after=0, add.to.row=list(pos=list(-1, -1, -1, -1, -1, -1, nrow(dfComb), nrow(dfComb)),
                                                      command=c(" \\toprule",
                                                                " \\multicolumn{1}{l}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.5$} &",
                                                                " \\multicolumn{1}{c}{} &",
                                                                " \\multicolumn{2}{c}{$\\rho = 0.9$} \\\\",
                                                                " \\cmidrule(l r){2-3} \\cmidrule(l r){5-6}",
                                                                " \\bottomrule \\\\[-.3cm]",
                                                                paste0(" \\multicolumn{6}{l}{\\parbox{9cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))


}
