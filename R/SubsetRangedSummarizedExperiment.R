
#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetRangedSummarizedExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetRangedSummarizedExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname show
setMethod(
  f = "show",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(object)
  {
    .show(object)
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname setSubsetAssay
setMethod(
  f = "setSubsetAssay",
  signature = c(
    x ="SubsetRangedSummarizedExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .setSubsetAssay(x,
                    subsetName,
                    inputMatrix,
                    subsetAssayName)
  }
)

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetRangedSummarizedExperiment", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname getSubsetAssay
setMethod(
  f = "getSubsetAssay",
  signature = c(
    x ="SubsetRangedSummarizedExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName)
  {
    assay(x = x, i = subsetName)
  }
)
