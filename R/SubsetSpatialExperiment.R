# Methods for SubsetSpatialExperiment only. The methods it shares with the
# other Subset* classes are defined in R/SubsetMethods.R

#' @rdname subsetSpatialCoords
setMethod(
  f = "subsetSpatialCoords",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, subsetName) {
    .subsetSpatialCoords(x, subsetName)
  }
)

#' @rdname subsetSpatialCoords
setReplaceMethod(
  f = "subsetSpatialCoords",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, subsetName, value) {
    .subsetSpatialCoords(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetImgData
setMethod(
  f = "subsetImgData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, subsetName) {
    .subsetImgData(x, subsetName)
  }
)

#' @rdname subsetSpatialCoords
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
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, subsetName, parentSpatialData) {
    .subsetSpatialData(x, subsetName, parentSpatialData)
  }
)

#' @rdname subsetSpatialData
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
  signature = c(x = "SubsetSpatialExperiment", subsetName = "ANY",
                value = "ANY"),
  definition = function(x, subsetName, value) {
    .subsetSpatialData(x, subsetName) <- value
    return(x)
  }
)
