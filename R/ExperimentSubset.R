#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndex = "numeric",
                                colIndex = "numeric"
                              )
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(
  subsetName = "subset",
  rowIndex = 0,
  colIndex = 0,
  ...)
{
  .SingleCellSubset(subsetName = subsetName,
                    rowIndex = rowIndex,
                    colIndex = colIndex)
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

#' @export
setGeneric(name = "subsetAssay",
           def = function(object, subsetName, rowIndex, colIndex)
           {
             standardGeneric("subsetAssay")
           }
)

#' @export
setMethod(f = "subsetAssay",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, rowIndex, colIndex)
            {
              scs <- SingleCellSubset(
                subsetName = subsetName,
                rowIndex = rowIndex,
                colIndex = colIndex)
              object@subsets[[subsetName]] <- scs
              return(object)
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
                  "subsets(", length(names(scs@subsets)), "): ",
                  sep="")
              cat(
                  names(object@subsets))
          }
)

#' @export
#' @importMethodsFrom SummarizedExperiment assay
setMethod("assay", c("ExperimentSubset", "numeric"), function(x, i, ...) {
  out <- callNextMethod()
  out <- out[1:5, 1:5]
  out
})
