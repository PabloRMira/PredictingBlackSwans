% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logDiff.R
\name{logDiffFun}
\alias{logDiffFun}
\title{Compute log-differences}
\usage{
logDiffFun(inputData, inputVars)
}
\arguments{
\item{inputData}{A data.table.}

\item{inputVars}{Variables for which the log-differences will be calculated.}
}
\value{
Additional variables in the data.table. These have the names
of the input variables with a "_lDiff"-ending.
}
\description{
Compute log-differences \eqn{log(y_t) - log(y_{t-1})}
}
\examples{
logDiff(data=dat, inputVars=c("gdp", "revenue"))
}
\keyword{PredictingBlackSwans}
