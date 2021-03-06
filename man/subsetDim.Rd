% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R,
%   R/SubsetRangedSummarizedExperiment.R, R/SubsetSingleCellExperiment.R,
%   R/SubsetSummarizedExperiment.R, R/SubsetSpatialExperiment.R,
%   R/SubsetTreeSummarizedExperiment.R
\name{subsetDim}
\alias{subsetDim}
\alias{subsetDim,SubsetRangedSummarizedExperiment,character-method}
\alias{subsetDim,SubsetSingleCellExperiment,character-method}
\alias{subsetDim,SubsetSummarizedExperiment,character-method}
\alias{subsetDim,SubsetSpatialExperiment,character-method}
\alias{subsetDim,SubsetTreeSummarizedExperiment,character-method}
\title{Get dimensions of subsets in ExperimentSubset objects}
\usage{
subsetDim(x, subsetName)

\S4method{subsetDim}{SubsetRangedSummarizedExperiment,character}(x, subsetName)

\S4method{subsetDim}{SubsetSingleCellExperiment,character}(x, subsetName)

\S4method{subsetDim}{SubsetSummarizedExperiment,character}(x, subsetName)

\S4method{subsetDim}{SubsetSpatialExperiment,character}(x, subsetName)

\S4method{subsetDim}{SubsetTreeSummarizedExperiment,character}(x, subsetName)
}
\arguments{
\item{x}{\code{ExperimentSubset} Input \code{ExperimentSubset} object.}

\item{subsetName}{\code{character(1)} Name of the subset to retrieve the
dimensions from.}
}
\value{
A \code{vector} containing the dimensions of the specified subset
  i.e. the number of rows and the number of columns in the subset.
}
\description{
Retrieves the dimensions of the specified subset in an
  \code{ExperimentSubset} object.
}
\examples{
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetDim(es, "subset1")
}
