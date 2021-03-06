% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R,
%   R/SubsetRangedSummarizedExperiment.R, R/SubsetSingleCellExperiment.R,
%   R/SubsetSummarizedExperiment.R, R/SubsetSpatialExperiment.R,
%   R/SubsetTreeSummarizedExperiment.R
\name{setSubsetAssay}
\alias{setSubsetAssay}
\alias{setSubsetAssay,SubsetRangedSummarizedExperiment,character,ANY,character-method}
\alias{setSubsetAssay,SubsetSingleCellExperiment,character,ANY,character-method}
\alias{setSubsetAssay,SubsetSummarizedExperiment,character,ANY,character-method}
\alias{setSubsetAssay,SubsetSpatialExperiment,character,ANY,character-method}
\alias{setSubsetAssay,SubsetTreeSummarizedExperiment,character,ANY,character-method}
\title{Method for storing new assays inside subsets in ExperimentSubset 
  objects}
\usage{
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName)

\S4method{setSubsetAssay}{SubsetRangedSummarizedExperiment,character,ANY,character}(x, subsetName, inputMatrix, subsetAssayName)

\S4method{setSubsetAssay}{SubsetSingleCellExperiment,character,ANY,character}(x, subsetName, inputMatrix, subsetAssayName)

\S4method{setSubsetAssay}{SubsetSummarizedExperiment,character,ANY,character}(x, subsetName, inputMatrix, subsetAssayName)

\S4method{setSubsetAssay}{SubsetSpatialExperiment,character,ANY,character}(x, subsetName, inputMatrix, subsetAssayName)

\S4method{setSubsetAssay}{SubsetTreeSummarizedExperiment,character,ANY,character}(x, subsetName, inputMatrix, subsetAssayName)
}
\arguments{
\item{x}{\code{ExperimentSubset} Specify the input object.}

\item{subsetName}{\code{character(1)} Specify the name of the existing subset
inside which the new subset \code{assay} should be stored.}

\item{inputMatrix}{\code{dgCMatrix} The input subset \code{assay}.}

\item{subsetAssayName}{\code{character(1)} Specify the name of the new
\code{assay} against the \code{inputMatrix} parameter.}
}
\value{
Updated \code{ExperimentSubset} object with the new \code{assay}
  stored inside the specified subset.
}
\description{
Store a new subset \code{assay} inside a specified subset in the
  input \code{ExperimentSubset} object.
}
\examples{
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
counts1p <- assay(es, "subset1")
counts1p[,] <- counts1p[,] + 1
es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
es
}
