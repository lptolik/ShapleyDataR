% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fbaShapley.R
\name{fbaShapley}
\alias{fbaShapley}
\title{Shapley values for FBA model.}
\usage{
fbaShapley(
  mod,
  reId,
  tol = 0.05,
  logName = "fbaShapley.log",
  tmpSave = 500,
  perfTolerance = 1e-05,
  cacheDepth = 3
)
}
\arguments{
\item{mod}{FBA model in Sybil format.}

\item{reId}{list of reaction Ids to be used in calculations.}

\item{tol}{tolerance for Shapley value convergence}

\item{logName}{name of the file to append log records to}

\item{tmpSave}{frequency of tmp results saving}

\item{perfTolerance}{tolerance for fail early check}

\item{cacheDepth}{how many first knockout have to be cached}
}
\value{
object of class \code{shapley}
}
\description{
\code{fbaShapley} return object of class \code{shapley}, which contains
\describe{
  \item{reId}{Reaction IDs as provided by argument}
  \item{shapley}{Shapley values, vector with reId as names}
  \item{sd}{Standard deviation for the Shapley values}
  \item{p}{Probability that the \code{shapley} contains real Shapley values.}
}
}
