
#' @rdname ExperimentSubset-metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname ExperimentSubset-assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSingleCellExperiment", i = "character"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .assay(x, i, withDimnames = withDimnames, ...)
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetSingleCellExperiment"
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
  signature = "SubsetSingleCellExperiment",
  definition = function(x, subsetName = NULL)
  {
    .subsetAssayNames(x, subsetName)
  }
)

#' @rdname SubsetSingleCellExperiment-class
#' @param object The object to print, with its subsets and their assays.
setMethod(
  f = "show",
  signature = "SubsetSingleCellExperiment",
  definition = function(object)
  {
    .show(object)
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname setSubsetAssay
setMethod(
  f = "setSubsetAssay",
  signature = c(
    x ="SubsetSingleCellExperiment",
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
                 c("SubsetSingleCellExperiment", "character"),
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
  signature = "SubsetSingleCellExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSingleCellExperiment", "character"),
  definition = function(x, subsetName)
  {
    .isSubset(x, subsetName)
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname ExperimentSubset-altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname ExperimentSubset-altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname ExperimentSubset-altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname ExperimentSubset-altExp
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

#' @rdname ExperimentSubset-altExp
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

#' @rdname ExperimentSubset-altExp
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSingleCellExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-altExp
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "character"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-altExp
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "missing"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-altExp
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "numeric"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-altExp
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExps(x, withDimnames = withDimnames, withColData = withColData,
             ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-metadata
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowData)
  {
    .subsetRowData(x, subsetName, parentRowData)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColData)
  {
    .subsetColData(x, subsetName, parentColData)
  }
)

#' @rdname subsetColData
setReplaceMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetColData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRowData
setReplaceMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetRowData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname ExperimentSubset-reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "missing"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname ExperimentSubset-reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname ExperimentSubset-reducedDim
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "character"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-reducedDim
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "missing"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-reducedDim
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "numeric"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-reducedDim
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, use.names = TRUE, ...) {
    .rowData(x, use.names = use.names, ...)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment", value = "NULL"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname ExperimentSubset-rowData-colData
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment", value = "data.frame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
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
    x ="SubsetSingleCellExperiment",
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
    x = "SubsetSingleCellExperiment",
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
    x = "SubsetSingleCellExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName,
                        value,
                        unlink = FALSE)
  {
    .setSubsetExperiment(x, subsetName, value, unlink = unlink)
  }
)

#' @rdname reducedDimNames-set
setReplaceMethod(
  f = "reducedDimNames",
  signature = "ANY",
  definition = function(x, subsetName, value)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      SingleCellExperiment::reducedDimNames(.internalAssay(.subsets(x)[[subsetName]])) <-
        value
    }
    else{
      SingleCellExperiment::reducedDimNames(x) <- value
    }
    x
  }
)

#' @rdname reducedDimNames
setMethod(
  f = "reducedDimNames",
  signature = c("ANY"),
  definition = function(x, ...)
  {
    arglist <- list(...)
    if(!"subsetName" %in% names(arglist))
      return(SingleCellExperiment::reducedDimNames(x))
    subsetName = arglist[["subsetName"]]
    .isSubset(x, subsetName)
    SingleCellExperiment::reducedDimNames(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname ExperimentSubset-subsetting
setMethod(
  f = "[",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, i, j, ..., drop = TRUE) {
    out <- callNextMethod()
    .remapSubsets(x, out, i, j)
  }
)

#' @rdname ExperimentSubset-combine
setMethod(
  f = "cbind",
  signature = "SubsetSingleCellExperiment",
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
  signature = "SubsetSingleCellExperiment",
  definition = function(..., deparse.level = 1) {
    objects <- list(...)
    .checkCombine(objects, "rows")
    out <- callNextMethod()
    .combineSubsets(out, objects, "rows")
  }
)
