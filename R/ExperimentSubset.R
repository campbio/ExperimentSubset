#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndices = "numeric",
                                colIndices = "numeric"
                              )
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(
  subsetName = "subset",
  rowIndices = NULL,
  colIndices = NULL,
  ...)
{
  .SingleCellSubset(subsetName = subsetName,
                    rowIndices = rowIndices,
                    colIndices = colIndices)
}

#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.ExperimentSubset <- setClass("ExperimentSubset",
                          slots = representation(
                            subsets = "list"
                          ),
                          contains = "SingleCellExperiment"
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
ExperimentSubset <- function(
  subsets = list(),
  ...)
{
  se <- SingleCellExperiment::SingleCellExperiment(...)
  .ExperimentSubset(se,
                    subsets = subsets)
}

#' @title subsetAssay
#'
#' @param object \code{ExperimentSubset}, \code{SingleCellExperiment} or \code{SummarizedExperiment} object.
#'
#' @param subsetName Specify the name of the new subset.
#' @param subsetRows Specify which rows/features to subset from the input object. Either indices or names.
#' @param subsetCols Specify which columns/cells to subset from the input object. Either indices or names.
#'
#' @export
setGeneric(name = "subsetAssay",
           def = function(object, subsetName, subsetRows, subsetCols)
           {
             standardGeneric("subsetAssay")
           }
)

#' @export
setMethod(f = "subsetAssay",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, subsetRows, subsetCols)
            {
            if(is.character(subsetRows)){
              subsetRows <- match(subsetRows, rownames(assay(object, "counts")))
            }
            if(is.character(subsetCols)){
              subsetCols <- match(subsetCols, colnames(assay(object, "counts")))
            }
              scs <- SingleCellSubset(
                subsetName = subsetName,
                rowIndices = subsetRows,
                colIndices = subsetCols)
              object@subsets[[subsetName]] <- scs
              return(object)
          }
)

#' @export
setGeneric(name = "subsetNames",
           def = function(object)
           {
             standardGeneric("subsetNames")
           }
)

#' @export
setMethod(f = "subsetNames",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            return(names(object@subsets))
          }
)

#' @export
#' @importMethodsFrom SingleCellExperiment show
setMethod(f = "show",
          signature = "ExperimentSubset",
          definition = function(object)
            {
              callNextMethod()
              cat(
                  "subsets(", length(subsetNames(object)), "): ",
                  sep="")
              cat(
                  subsetNames(object)
                  )
          }
)

#' @export
#' @importMethodsFrom SummarizedExperiment assay
setMethod("assay", c("ExperimentSubset", "character"), function(x, i, ...) {
  if(i %in% subsetNames(x)){
    subsetName = i
    i = "counts"
    out <- callNextMethod()
    out <- out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
  }
  else{
    out <- callNextMethod()
  }
  out
})
