setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))

#' An S4 class to manage subset representation.
#'
#' @slot subsetName \code{character(1)} Name of the subset.
#' @slot rowIndices \code{vector("numeric")} Indices of the rows to include in
#'   the subset.
#' @slot colIndices \code{vector("numeric")} Indices of the columns to include
#'   in the subset.
#' @slot parentAssay \code{character(1)} Name of the parent of this subset.
#' @slot internalAssay \code{SummarizedExperiment} An internal experiment object
#'   to store additional subset data.
#' @import methods
.AssaySubset <- setClass(
  "AssaySubset",
  slots = representation(
    subsetName = "character",
    rowIndices = "NullOrNumeric",
    colIndices = "NullOrNumeric",
    parentAssay = "NullOrCharacter",
    internalAssay = "ANY"
  )
)

#' An S4 class for \code{SummarizedExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SubsetSummarizedExperiment <- setClass(
  Class = "SubsetSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = SimpleList(subsets = SimpleList()),
  contains = "SummarizedExperiment"
)

#' An S4 class for \code{RangedSummarizedExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.SubsetRangedSummarizedExperiment <- setClass(
  Class = "SubsetRangedSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = SimpleList(subsets = SimpleList()),
  contains = "RangedSummarizedExperiment"
)

#' An S4 class for \code{SingleCellExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SubsetSingleCellExperiment <- setClass(
  Class = "SubsetSingleCellExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = SimpleList(subsets = SimpleList()),
  contains = "SingleCellExperiment"
)

#' An S4 class for \code{TreeSummarizedExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @export
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
.SubsetTreeSummarizedExperiment <- setClass(
  Class = "SubsetTreeSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = SimpleList(subsets = SimpleList()),
  contains = "TreeSummarizedExperiment"
)

#' An S4 class for \code{SpatialExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SpatialExperiment SpatialExperiment
.SubsetSpatialExperiment <- setClass(
  Class = "SubsetSpatialExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = SimpleList(subsets = SimpleList()),
  contains = "SpatialExperiment"
)
