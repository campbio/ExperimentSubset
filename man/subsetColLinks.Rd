% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R,
%   R/SubsetTreeSummarizedExperiment.R
\name{subsetColLinks}
\alias{subsetColLinks}
\alias{subsetColLinks,SubsetTreeSummarizedExperiment,character-method}
\title{Accessor method for colLinks from subsets in ExperimentSubset objects}
\usage{
subsetColLinks(x, subsetName, parentColLinkData)

\S4method{subsetColLinks}{SubsetTreeSummarizedExperiment,character}(x, subsetName, parentColLinkData)
}
\arguments{
\item{x}{\code{ExperimentSubset} Input \code{ExperimentSubset} object.}

\item{subsetName}{\code{character(1)} Name of the subset to get
\code{colLinks} from.}

\item{parentColLinkData}{\code{logical(1)} Logical value indicating if parent
colLinks should be combined or not. Default \code{FALSE}.}
}
\value{
The \code{colLinks} from input object.
}
\description{
Get \code{colLinks} from a subset.
}
