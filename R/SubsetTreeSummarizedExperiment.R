#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetTreeSummarizedExperiment", i = "character"),
  definition = function(x, i, withDimnames = TRUE, ...) {
    .assay(x, i, withDimnames = withDimnames, ...)
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetTreeSummarizedExperiment"
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
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x, subsetName = NULL)
  {
    .subsetAssayNames(x, subsetName)
  }
)

#' @rdname show
setMethod(
  f = "show",
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(object)
  {
    .show(object)
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname setSubsetAssay
setMethod(
  f = "setSubsetAssay",
  signature = c(
    x ="SubsetTreeSummarizedExperiment",
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

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetTreeSummarizedExperiment", "character"),
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
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetTreeSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .isSubset(x, subsetName)
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

#' @rdname altExpNames<-
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetTreeSummarizedExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "character"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "missing"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetTreeSummarizedExperiment", e = "numeric"),
  definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
            ...) <- value
    return(x)
  }
)

#' @rdname altExps<-
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, withDimnames = TRUE, withColData = FALSE, ...,
                        value) {
    .altExps(x, withDimnames = withDimnames, withColData = withColData,
             ...) <- value
    return(x)
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowData)
  {
    .subsetRowData(x, subsetName, parentRowData)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColData)
  {
    .subsetColData(x, subsetName, parentColData)
  }
)

#' @rdname subsetColData
#' @usage NULL
setReplaceMethod(
  f = "subsetColData",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetColData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRowData
#' @usage NULL
setReplaceMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetRowData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "missing"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "character"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "missing"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetTreeSummarizedExperiment", type = "numeric"),
  definition = function(x, type, withDimnames = TRUE, ..., value) {
    .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
    return(x)
  }
)

#' @rdname reducedDims
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, use.names = TRUE, ...) {
    .rowData(x, use.names = use.names, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetTreeSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetTreeSummarizedExperiment", value = "NULL"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetTreeSummarizedExperiment",
                        value = "data.frame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
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
    x ="SubsetTreeSummarizedExperiment",
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
    x = "SubsetTreeSummarizedExperiment",
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
    x = "SubsetTreeSummarizedExperiment",
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

#' @rdname subsetRowLinks
setMethod(
  f = "subsetRowLinks",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowLinkData)
  {
    .subsetRowLinkData(x, subsetName, parentRowLinkData)
  }
)

#' @rdname subsetColLinks
setMethod(
  f = "subsetColLinks",
  signature = c(x = "SubsetTreeSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColLinkData)
  {
    .subsetColLinkData(x, subsetName, parentColLinkData)
  }
)

#' @rdname ExperimentSubset-subsetting
setMethod(
  f = "[",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, i, j, ..., drop = TRUE) {
    out <- callNextMethod()
    .remapSubsets(x, out, i, j)
  }
)

#' @rdname ExperimentSubset-combine
setMethod(
  f = "cbind",
  signature = "SubsetTreeSummarizedExperiment",
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
  signature = "SubsetTreeSummarizedExperiment",
  definition = function(..., deparse.level = 1) {
    objects <- list(...)
    .checkCombine(objects, "rows")
    out <- callNextMethod()
    .combineSubsets(out, objects, "rows")
  }
)
