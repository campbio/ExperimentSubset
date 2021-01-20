
#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetColData
setReplaceMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetColData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRowData
setReplaceMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetRowData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname getSubsetAssay
setMethod(
  f = "getSubsetAssay",
  signature = c(
    x ="SubsetSummarizedExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName)
  {
    assay(x = x, i = subsetName)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSummarizedExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetSummarizedExperiment"
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
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname show
setMethod(
  f = "show",
  signature = "SubsetSummarizedExperiment",
  definition = function(object)
  {
    .show(object)
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname setSubsetAssay
setMethod(
  f = "setSubsetAssay",
  signature = c(
    x ="SubsetSummarizedExperiment",
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
                 c("SubsetSummarizedExperiment", "character"),
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
  signature = "SubsetSummarizedExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowData)
  {
    .subsetRowData(x, subsetName, parentRowData)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColData)
  {
    .subsetColData(x, subsetName, parentColData)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)
