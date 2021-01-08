#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.ExperimentSubsetSCE <- setClass(
  Class = "ExperimentSubsetSCE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SingleCellExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ExperimentSubsetSE <- setClass(
  Class = "ExperimentSubsetSE",
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
.ExperimentSubsetRSE <- setClass(
  Class = "ExperimentSubsetRSE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "RangedSummarizedExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SpatialExperiment VisiumExperiment
.ExperimentSubsetVE <- setClass(
  Class = "ExperimentSubsetVE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "VisiumExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SpatialExperiment SpatialExperiment
.ExperimentSubsetSP <- setClass(
  Class = "ExperimentSubsetSP",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SpatialExperiment"
)


#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the \code{ExperimentSubset} object, either through manually specifying the \code{assays}, \code{rowData}, \code{colData} or directly by passing either a \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or objects inherited by these classes. A subset can also be directly created by pasing a named \code{list} to the \code{subset} parameter. This named \code{list} should have parameter values named as \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @param object A \code{SingleCellExperiment} or \code{SummarizedExperiment} object if direct conversion is required.
#' @param ... Additional parameters passed to \code{SingleCellExperiment} constructor.
#' @param subset A named \code{list} if a subset should be created from within the constructor. Named parameters in this list should be \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @return A \code{ExperimentSubset} object.
#' @export
#' @import BiocStyle
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
ExperimentSubset <- function(object,
                             ...,
                             subset = list(
                               subsetName = NA,
                               rows = NA,
                               cols = NA,
                               parentAssay = NA
                             ))
{
  if (!missing(object)) {
    # if (inherits(object, "SummarizedExperiment")) {
    #   object <- as(object, "SingleCellExperiment")
    #   es <- .ExperimentSubset(object)
    # }
    if(inherits(object, "VisiumExperiment")){
      es <- .ExperimentSubsetVE(object)
    }
    else if(inherits(object, "SpatialExperiment")){
      es <- .ExperimentSubsetSP(object)
    }
    else if(inherits(object, "SingleCellExperiment")){
      es <- .ExperimentSubsetSCE(object)
    }
    else if(inherits(object, "RangedSummarizedExperiment")){
      es <- .ExperimentSubsetRSE(object)
    }
    else if(inherits(object, "SummarizedExperiment")){
      es <- .ExperimentSubsetSE(object)
    }
  }
  else{
    se <- SingleCellExperiment::SingleCellExperiment(...)
    es <- .ExperimentSubset(se)
  }
  if (!anyNA(subset)) {
    es <- ExperimentSubset::createSubset(
      es,
      subsetName = subset$subsetName,
      rows = subset$rows,
      cols = subset$cols,
      parentAssay = subset$parentAssay
    )
  }
  es
}

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "ExperimentSubsetSE"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "ExperimentSubsetSCE"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  metadata(x@subsets[[subsetName]]@internalAssay)
}
