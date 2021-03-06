% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R,
%   R/SubsetRangedSummarizedExperiment.R, R/SubsetSingleCellExperiment.R,
%   R/SubsetSummarizedExperiment.R, R/SubsetSpatialExperiment.R,
%   R/SubsetTreeSummarizedExperiment.R
\name{subsetAssayNames}
\alias{subsetAssayNames}
\alias{subsetAssayNames,SubsetRangedSummarizedExperiment-method}
\alias{subsetAssayNames,SubsetSingleCellExperiment-method}
\alias{subsetAssayNames,SubsetSummarizedExperiment-method}
\alias{subsetAssayNames,SubsetSpatialExperiment-method}
\alias{subsetAssayNames,SubsetTreeSummarizedExperiment-method}
\title{Name retrieval method for all subset assays in ExperimentSubset
  objects}
\usage{
subsetAssayNames(x)

\S4method{subsetAssayNames}{SubsetRangedSummarizedExperiment}(x)

\S4method{subsetAssayNames}{SubsetSingleCellExperiment}(x)

\S4method{subsetAssayNames}{SubsetSummarizedExperiment}(x)

\S4method{subsetAssayNames}{SubsetSpatialExperiment}(x)

\S4method{subsetAssayNames}{SubsetTreeSummarizedExperiment}(x)
}
\arguments{
\item{x}{\code{ExperimentSubset} Input \code{ExperimentSubset} object.}
}
\value{
A \code{vector} containing the names of the subsets and the subset
  assays available in the input \code{ExperimentSubset} object.
}
\description{
Retrieves the names of all the subsets as well as the subset
  assays.
}
\examples{
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
assay(es, "subset1",
subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
subsetAssayNames(es)
}
