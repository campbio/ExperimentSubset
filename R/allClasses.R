setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))

#' An S4 class to manage subset representation.
#'
#' Used inside the \code{Subset*} classes; not for direct use.
#' @keywords internal
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
#' @exportClass SubsetSummarizedExperiment
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SubsetSummarizedExperiment <- setClass(
  Class = "SubsetSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = prototype(subsets = SimpleList()),
  contains = "SummarizedExperiment"
)

#' An S4 class for \code{RangedSummarizedExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @exportClass SubsetRangedSummarizedExperiment
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.SubsetRangedSummarizedExperiment <- setClass(
  Class = "SubsetRangedSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = prototype(subsets = SimpleList()),
  contains = "RangedSummarizedExperiment"
)

#' An S4 class for \code{SingleCellExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @exportClass SubsetSingleCellExperiment
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SubsetSingleCellExperiment <- setClass(
  Class = "SubsetSingleCellExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = prototype(subsets = SimpleList()),
  contains = "SingleCellExperiment"
)

#' An S4 class for \code{TreeSummarizedExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @exportClass SubsetTreeSummarizedExperiment
#' @import methods
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
.SubsetTreeSummarizedExperiment <- setClass(
  Class = "SubsetTreeSummarizedExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = prototype(subsets = SimpleList()),
  contains = "TreeSummarizedExperiment"
)

#' An S4 class for \code{SpatialExperiment} objects with added support for
#'   subsets.
#' @slot subsets A list of \code{AssaySubset} objects.
#' @exportClass SubsetSpatialExperiment
#' @import methods
#' @importClassesFrom SpatialExperiment SpatialExperiment
.SubsetSpatialExperiment <- setClass(
  Class = "SubsetSpatialExperiment",
  slots = representation(subsets = "SimpleList"),
  prototype = prototype(subsets = SimpleList()),
  contains = "SpatialExperiment"
)

#' The class of every ExperimentSubset object
#'
#' A virtual class that the five \code{Subset*} classes extend:
#' \code{\linkS4class{SubsetSummarizedExperiment}},
#' \code{\linkS4class{SubsetRangedSummarizedExperiment}},
#' \code{\linkS4class{SubsetSingleCellExperiment}},
#' \code{\linkS4class{SubsetTreeSummarizedExperiment}} and
#' \code{\linkS4class{SubsetSpatialExperiment}}. Use
#' \code{is(x, "ExperimentSubset")} to test whether \code{x} is an
#' \code{ExperimentSubset} object, whichever class it extends. A
#' \code{SubsetSpatialExperiment}, for example, is not a
#' \code{SubsetSingleCellExperiment}, though it extends
#' \code{SingleCellExperiment}.
#'
#' The class is virtual, so there are no objects of the class itself. Create
#' an \code{ExperimentSubset} object with \code{\link{ExperimentSubset}()}.
#' @name ExperimentSubset-class
#' @aliases ExperimentSubset-class
#' @exportClass ExperimentSubset
#' @examples
#' # simulate a small count matrix
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("spot", 1:20)
#'
#' # a SpatialExperiment, and an ExperimentSubset object made from it
#' spe <- SpatialExperiment(assays = list(counts = counts))
#' es <- ExperimentSubset(spe)
#'
#' # es is an ExperimentSubset, whichever class it extends
#' is(es, "ExperimentSubset")
#' is(es, "SubsetSingleCellExperiment")
#'
#' # the SpatialExperiment it was made from is not
#' is(spe, "ExperimentSubset")
setClassUnion(
  "ExperimentSubset",
  c(
    "SubsetSummarizedExperiment",
    "SubsetRangedSummarizedExperiment",
    "SubsetSingleCellExperiment",
    "SubsetTreeSummarizedExperiment",
    "SubsetSpatialExperiment"
  )
)
