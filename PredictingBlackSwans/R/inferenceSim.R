#' Inference Simulation Study
#'
#' Replicate the results of our simulations part regarding inference. IMPORTANT NOTE: This
#' function is thought to be run as a script. See the details.
#' @param path Path to export the results.
#' @param parallel Should parallel computing be used? Note: It only works for UNIX systems.
#' @param ncores How many cores should be used for parallel computing?
#' @param seed Seed for reproducibility purposes.
#' @details This function may take to long to run for computers with few kernels or
#' for Windows-computers. Therefore we suggest to run this function as an script to split
#' the computation of the simulations in several days.
#' @examples inferenceSim()
#' @export
#' @keywords PredictingBlackSwans
inferenceSim <- function(path        = getwd(),
                         parallel    = T,
                         ncores      = getOption("mc.cores", 2L),
                         seed        = 912
) {
  # Save current directory
  curDir <- getwd()

  # Create folder structure if it does not exist
  simPath <- file.path(path, "Inference_Simulations")
  simTex <- file.path(simPath, "TeX")
  simHtml <- file.path(simPath, "Html")

  # Create new folders
  dir.create(simPath, showWarnings = FALSE)
  dir.create(simTex, showWarnings = FALSE)
  dir.create(simHtml, showWarnings = FALSE)

  # Register kernels for parallel computing if parallel = TRUE
  if (parallel == TRUE) doMC::registerDoMC(cores=ncores)

  # n = 100, p = 150, rho = 0.5
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                   n=100, p=150, rho=0.5,
                   parallel=parallel, ncores=ncores, seed=seed)
  # n = 100, p = 150, rho = 0.9
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                   n=100, p=150, rho=0.9,
                   parallel=parallel, ncores=ncores, seed=seed)
  # n = 200, p = 250, rho = 0.5
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                   n=200, p=250, rho=0.5,
                   parallel=parallel, ncores=ncores, seed=seed)
  # n = 200, p = 250, rho = 0.9
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                   n=200, p=250, rho=0.9,
                   parallel=parallel, ncores=ncores, seed=seed)

  ########################################################################################

  # Additional Simulations with n = 100, p = 500 to assess the limitations of the
  # desparsified Lasso estimator.

  # n = 100, p = 500, rho = 0.5
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                                         n=100, p=500, rho=0.5,
                                         parallel=parallel, ncores=ncores, seed=seed)
  # n = 100, p = 500, rho = 0.9
  PredictingBlackSwans::inferenceSimMain(path=simPath,
                                         n=100, p=500, rho=0.9,
                                         parallel=parallel, ncores=ncores, seed=seed)

  ########################################################################################

  # Additional simulations to show that the proposed scaling in van de Geer (2014) using
  # the "desparsified outer product of the gradient" has not good finite sample properties
  # compared to our preferred implementation using the "desparsified Hessian".

  # n = 100, p = 150, rho = 0.5; as in the original paper (van de Geer(2014))
  PredictingBlackSwans::inferenceSimMainPaper(path=simPath,
                                         n=100, p=150, rho=0.5,
                                         parallel=parallel, ncores=ncores, seed=seed)
  # n = 100, p = 150, rho = 0.9; as in the original paper (van de Geer(2014))
  PredictingBlackSwans::inferenceSimMainPaper(path=simPath,
                                         n=100, p=150, rho=0.9,
                                         parallel=parallel, ncores=ncores, seed=seed)
  # n = 200, p = 250, rho = 0.5; as in the original paper (van de Geer(2014))
  PredictingBlackSwans::inferenceSimMainPaper(path=simPath,
                                              n=200, p=250, rho=0.5,
                                              parallel=parallel, ncores=ncores, seed=seed)
  # n = 200, p = 250, rho = 0.9; as in the original paper (van de Geer(2014))
  PredictingBlackSwans::inferenceSimMainPaper(path=simPath,
                                              n=200, p=250, rho=0.9,
                                              parallel=parallel, ncores=ncores, seed=seed)

  # Generate the extra plot with the standard error of the estimates against the estimates
  PredictingBlackSwans::inferenceMeanSePlot(path=simPath)

  # Generate the results tables in the Latex-format as in the thesis
  PredictingBlackSwans::inferenceSimPrintResults(inPath=simPath, outPath=simTex)

  setwd(curDir)
}
