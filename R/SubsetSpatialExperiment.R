
#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSpatialExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetSpatialExperiment"
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
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname show
setMethod(
  f = "show",
  signature = "SubsetSpatialExperiment",
  definition = function(object)
  {
    .show(object)
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname setSubsetAssay
setMethod(
  f = "setSubsetAssay",
  signature = c(
    x ="SubsetSpatialExperiment",
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
                 c("SubsetSpatialExperiment", "character"),
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
  signature = "SubsetSpatialExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSpatialExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

#' @rdname altExpNames<-
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSpatialExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "character"),
  definition = function(x, e,  ..., value) {
    .altExp(x, e, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "missing"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "numeric"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExps<-
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ..., value) {
    .altExps(x, ...) <- value
    return(x)
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentRowData)
  {
    .subsetRowData(x, subsetName, parentRowData)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentColData)
  {
    .subsetColData(x, subsetName, parentColData)
  }
)

#' @rdname subsetColData
setReplaceMethod(
  f = "subsetColData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetColData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRowData
setReplaceMethod(
  f = "subsetRowData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character", value = "DataFrame"),
  definition = function(x, subsetName, value)
  {
    .subsetRowData(x, subsetName) <- value
    return(x)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "missing"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "character"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "missing"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "numeric"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDims
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSpatialExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSpatialExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
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
    x ="SubsetSpatialExperiment",
    subsetName = "character"),
  definition = function(x,
                        subsetName)
  {
    assay(x = x, i = subsetName)
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

#' @rdname subsetSpatialCoords
setMethod(
  f = "subsetSpatialCoords",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetSpatialCoords(x, subsetName)
  }
)

#' @rdname spatialCoords
setMethod(
  f = "spatialCoords",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .spatialCoords(x, ...)
  }
)

#' @rdname subsetSpatialData
setMethod(
  f = "subsetSpatialData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, parentSpatialData)
  {
    .subsetSpatialData(x, subsetName, parentSpatialData)
  }
)

#' @rdname spatialData
setMethod(
  f = "spatialData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .spatialData(x, ...)
  }
)

#' @rdname subsetSpatialData
setReplaceMethod(
  f = "subsetSpatialData",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character", value = "data.frame"),
  definition = function(x, subsetName, value)
  {
    .subsetSpatialData(x, subsetName) <- value
    return(x)
  }
)