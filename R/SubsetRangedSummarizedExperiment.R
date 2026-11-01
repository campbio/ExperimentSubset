
#' @rdname ExperimentSubset-metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname ExperimentSubset-assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetRangedSummarizedExperiment", i = "character"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .assay(x, i, withDimnames = withDimnames, ...)
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
  definition = function(x, subsetName = NULL)
  {
    .subsetAssayNames(x, subsetName)
  }
)

#' @rdname SubsetRangedSummarizedExperiment-class
#' @param object The object to print, with its subsets and their assays.
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
                        subsetAssayName,
                        unlink = FALSE)
  {
    .setSubsetAssay(x,
                    subsetName,
                    inputMatrix,
                    subsetAssayName,
                    unlink = unlink)
  }
)

#' @rdname ExperimentSubset-assay
setReplaceMethod("assay",
                 c("SubsetRangedSummarizedExperiment", "character"),
                 function(x,
                          i,
                          withDimnames = TRUE,
                          ...,
                          value) {
                   .assay(x, i, withDimnames = withDimnames, ...) <- value
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
    .isSubset(x, subsetName)
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

#' @rdname ExperimentSubset-metadata
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
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowData)
  {
    .subsetRowData(x, subsetName, parentRowData)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColData)
  {
    .subsetColData(x, subsetName, parentColData)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, use.names = TRUE, ...) {
    .rowData(x, use.names = use.names, ...)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "NULL"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment",
                        value = "data.frame"),
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

#' @rdname subsetColData
setReplaceMethod(
  f = "subsetColData",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetColData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRowData
setReplaceMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character", value = "DataFrame"),
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
    x ="SubsetRangedSummarizedExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName,
                        subsetAssayName = NULL)
  {
    if (is.null(subsetAssayName))
      return(assay(x = x, i = subsetName))
    assay(x = x, i = subsetName, subsetAssayName = subsetAssayName)
  }
)

#' @rdname getSubsetExperiment
setMethod(
  f = "getSubsetExperiment",
  signature = c(
    x = "SubsetRangedSummarizedExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName,
                        parentData = TRUE)
  {
    .getSubsetExperiment(x, subsetName, parentData = parentData)
  }
)

#' @rdname setSubsetExperiment
setMethod(
  f = "setSubsetExperiment",
  signature = c(
    x = "SubsetRangedSummarizedExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName,
                        value,
                        unlink = FALSE)
  {
    .setSubsetExperiment(x, subsetName, value, unlink = unlink)
  }
)

#' @rdname ExperimentSubset-subsetting
setMethod(
  f = "[",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, i, j, ..., drop = TRUE) {
    out <- callNextMethod()
    .remapSubsets(x, out, i, j)
  }
)

#' @rdname ExperimentSubset-combine
setMethod(
  f = "cbind",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(..., deparse.level = 1) {
    objects <- list(...)
    .checkCombine(objects, "cols")
    out <- callNextMethod()
    .combineSubsets(out, objects, "cols")
  }
)

#' @rdname ExperimentSubset-combine
setMethod(
  f = "rbind",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(..., deparse.level = 1) {
    objects <- list(...)
    .checkCombine(objects, "rows")
    out <- callNextMethod()
    .combineSubsets(out, objects, "rows")
  }
)
