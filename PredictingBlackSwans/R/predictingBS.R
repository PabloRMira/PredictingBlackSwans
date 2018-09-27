#' Predicting Black Swans
#'
#' Main function to reproduce the results of our paper regarding the
#' prediction analysis.
#' @param path Path to export the results
#' @param mccvnumber Number of Monte Carlo Cross-Validations.
#' @param cvfolds Number of folds for the Lasso.
#' @param seed Seed for reproducibility.
#' @param parallel Should parallel computing be used? (Only for UNIX computers)
#' @param ncores Number of cores for parallel computing.
#' @details A new directory "Prediction_Analysis" will be created
#' where the results will be placed in.
#' @export
#' @keywords PredictingBlackSwans
predictingBS <- function(path=getwd(),
                         mccvnumber=100,
                         cvfolds=5,
                         seed=813,
                         parallel=T,
                         ncores=2L) {
  # Get current directory
  curDir <- getwd()

  # Create folder structure
  predPath <- file.path(path, "Prediction_Analysis")
  mccvPath <- file.path(predPath, "MCCV_Analysis")
  testPath <- file.path(predPath, "Test_Analysis")

  dir.create(predPath, showWarnings = FALSE)
  dir.create(mccvPath, showWarnings = FALSE)
  dir.create(testPath, showWarnings = FALSE)

  # Set working directory for exporting the graphs
  setwd(mccvPath)

  # Performance Comparison
  outputList <- PredictingBlackSwans::MCCV_Analysis(mccvnumber = mccvnumber,
                                                    cvfolds    = cvfolds,
                                                    seed       = seed,
                                                    parallel   = parallel,
                                                    ncores     = ncores)
  # Table

  # Note
  mynote <- paste("Dependent variable: 2-year horizon before crisis (pre-crisis).",
                  "Out-of-sample AUC and confidence interval estimates",
                  "based on Monte Carlo Cross-Validation (Picard and Cook (1984);",
                  "Arlot and Celisse (2010) for the Lasso and LR and based on",
                  "out-of-bag data for the random forest. \\\\",
                  "$^\\S H_0: AUC_{Lasso} = AUC_{LR}$ as well as $AUC_{RF} = AUC_{Lasso}$",
                  "and $H_1: AUC_{Lasso} > AUC_{LR}$ as well as $AUC_{RF} > AUC_{Lasso}$ \\\\",
                  "$^\\star$ p-value $< 0.10$; $^{\\star\\star}$ p-value $< 0.05$; $^{\\star\\star\\star}$",
                  "p-value $< 0.01$ for $H_0: AUC = 0.5$ and $H_1: AUC > 0.5$")

  outputXtable <- xtable::xtable(outputList$outputLatex, type="latex",
                                 caption="Performance Comparison: Monte Carlo Cross Validation (MCCV)",
                                 display=c("s", "s", "s", "s", "d", "d", "d"),
                                 align="llccrrr")
  xtable::print.xtable(outputXtable, type="latex",  file="MCCV_Performance.tex", caption.placement = "top",
                       sanitize.text.function = function(x) {x}, include.rownames = F,
                       booktabs=T, replace=T, hline.after=0,
                       add.to.row=list(pos=list(-1, nrow(outputList$outputLatex), nrow(outputList$outputLatex)),
                                       command=c(" \\toprule \\toprule ",
                                                 " \\bottomrule \\\\[-.3cm]",
                                                 paste0(" \\multicolumn{6}{l}{\\parbox{16cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))

  # Plot
  ggsave("MCCV_ROC.pdf", plot=outputList$ROCPlot, units="cm", height=9.5, width=20)

  myvarlabels <- c(expression(paste("L1. \U0394 log(", C[HH], ") \u2a09 \u0394 log(C / GDP)")),
                   expression(paste("L1. \u0394 log(", C[M], ") \u2a09 \u0394 log(", C[B], " / GDP)")),
                   expression(paste(Delta, "log(", C[B], " / GDP)")),
                   "L2. \u0394 log(SPI) \u2a09 \u0394 log(C / GDP)",
                   expression(paste("L1. r \u2a09 \u0394 log(", C[B], " / GDP)")),
                   expression(paste("r \u2a09 \u0394 log(", C[B], " / GDP)")))
  p <- outputList$pVarImp +
    scale_x_discrete(labels=myvarlabels)
  ggsave("Variable_Importance.pdf", plot=p, units="cm", height=8, width=15, device=grDevices::cairo_pdf)

  myvarlabels2 <- c("L2. \u0394 log(SPI) \u2a09 \u0394 log(C / GDP)",
                    expression(paste("L1. \U0394 log(", C[HH], ") \u2a09 \u0394 log(C / GDP)")),
                    expression(paste("L1. \u0394 log(", C[M], ") \u2a09 \u0394 log(", C[B], " / GDP)")),
                    expression(paste("L1. r \u2a09 \u0394 log(", C[B], " / GDP)")),
                    expression(paste("r \u2a09 \u0394 log(", C[B], " / GDP)")),
                    expression(paste(Delta, "log(", C[B], " / GDP)")))
  myvarlabels2 <- myvarlabels2[length(myvarlabels2):1]

  # Example: Crisis probabilities and linear efffects
  p1 <- outputList$pCrisProbs$Spain
  p2 <- outputList$pLinEffects$Spain
  p2 <- outputList$pLinEffects$Spain +
    scale_fill_discrete(breaks=as.character(unique(p2$data$variable)),
                        labels=myvarlabels2)
  p2 <- p2 + labs(title="")
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g <- rbind(g1, g2, size="first")
  g$width <- grid::unit.pmax(g1$widths, g2$widths)
  grid::grid.newpage()
  grDevices::cairo_pdf(file="Spain_Crisis_Linear_Effects.pdf", width=21 / 2.54, height=21 / 2.54)
  grid::grid.draw(g)
  dev.off()

  ### Case Study: Predicting the global financial crisis (2007 / 2008) ###
  setwd(testPath)
  outputList <- PredictingBlackSwans::testAnalysis(mccvnumber = mccvnumber,
                                                   cvfolds    = cvfolds,
                                                   seed       = seed,
                                                   parallel   = parallel,
                                                   ncores     = ncores)

  # Save plots
  # Test ROC
  ggsave("ROC_Test.pdf", plot=outputList$pRocTest, units="cm", height=10, width=20) # Test ROC

  # Predictions
  ggsave("LR_Test_Probs.pdf", plot=outputList$pLogitProbs, units="cm", height=11.5, width=20) # Logistic Regression
  ggsave("Lasso_Test_Probs.pdf", plot=outputList$pLassoProbs, units="cm", height=12, width=20) # Lasso
  ggsave("Random_Forest_Test_Probs.pdf", plot=outputList$pRfProbs, units="cm", height=12, width=20) # Random Forest

  # Table
  # Note
  mynote <- paste("Dependent variable: 2-year horizon before crisis (pre-crisis).")

  outputXtable <- xtable::xtable(outputList$outputTableTest, type="latex",
                                 caption="Performance Comparison: Test Data (1998 -- 2016)",
                                 display=c("s", "s", "s", "d", "d", "d"),
                                 align="llcrrr")
  xtable::print.xtable(outputXtable, type="latex",  file="Test_Performance.tex", caption.placement = "top",
                       sanitize.text.function = function(x) {x}, include.rownames = F,
                       booktabs=T, replace=T, hline.after=0,
                       add.to.row=list(pos=list(-1, nrow(outputList$outputTableTest), nrow(outputList$outputTableTest)),
                                       command=c(" \\toprule \\toprule ",
                                                 " \\bottomrule \\\\[-.3cm]",
                                                 paste0(" \\multicolumn{5}{l}{\\parbox{16cm}{\\footnotesize \\textit{Note}: ", mynote, "}}"))))
}
