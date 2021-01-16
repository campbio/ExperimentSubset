#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SubsetSummarizedExperiment <- setClass(
  Class = "SubsetSummarizedExperiment",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SummarizedExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.SubsetRangedSummarizedExperiment <- setClass(
  Class = "SubsetRangedSummarizedExperiment",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "RangedSummarizedExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SubsetSingleCellExperiment <- setClass(
  Class = "SubsetSingleCellExperiment",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SingleCellExperiment"
)
