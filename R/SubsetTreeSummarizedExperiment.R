# Methods for SubsetTreeSummarizedExperiment only. The methods it shares
# with the other Subset* classes are defined in R/SubsetMethods.R

#' @rdname subsetRowLinks
setMethod(
  f = "subsetRowLinks",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, parentRowLinkData = TRUE) {
    .subsetLinks(x, subsetName, "row", parentRowLinkData)
  }
)

#' @rdname subsetColLinks
setMethod(
  f = "subsetColLinks",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, parentColLinkData = TRUE) {
    .subsetLinks(x, subsetName, "col", parentColLinkData)
  }
)

#' @rdname subsetRowTree
setMethod(
  f = "subsetRowTree",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, whichTree = 1) {
    .subsetTree(x, subsetName, "row", whichTree)
  }
)

#' @rdname subsetRowTree
setReplaceMethod(
  f = "subsetRowTree",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, value) {
    .subsetTree(x, subsetName, "row") <- value
    x
  }
)

#' @rdname subsetColTree
setMethod(
  f = "subsetColTree",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, whichTree = 1) {
    .subsetTree(x, subsetName, "col", whichTree)
  }
)

#' @rdname subsetColTree
setReplaceMethod(
  f = "subsetColTree",
  signature = signature(x = "SubsetTreeSummarizedExperiment"),
  definition = function(x, subsetName, value) {
    .subsetTree(x, subsetName, "col") <- value
    x
  }
)
