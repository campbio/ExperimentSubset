library(SingleCellExperiment)
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndex = "numeric",
                                colIndex = "numeric"
                              )
)

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


.ExperimentSubset <- setClass("ExperimentSubset",
                          slots = representation(
                            subsets = "list"
                          ),
                          contains = "SingleCellExperiment"
)

ExperimentSubset <- function(
  subsets = list(),
  ...)
{
  se <- SingleCellExperiment(...)
  .ExperimentSubset(se,
                    subsets = subsets)
}

setGeneric(name = "subsetAssay",
           def = function(object, subsetName, rowIndex, colIndex)
           {
             standardGeneric("subsetAssay")
           }
)

setMethod(f = "subsetAssay",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, rowIndex, colIndex)
          {
            scs <- SingleCellSubset(subsetName = subsetName, rowIndex = rowIndex, colIndex = colIndex)
            object@subsets[[subsetName]] <- scs
            return(object)
          }
)

setMethod(f = "show",
          signature = "ExperimentSubset",
          definition = function(object) {
  callNextMethod()
  cat(
    "subsets(0):",
    sep=""
  )
})

