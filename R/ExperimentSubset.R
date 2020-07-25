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

setGeneric(name="subsetAssay",
           def=function(obj, subsetName, rowIndex, colIndex)
           {
             standardGeneric("subsetAssay")
           }
)

setMethod(f="subsetAssay",
          signature="ExperimentSubset",
          definition=function(obj, subsetName, rowIndex, colIndex)
          {
            scs <- SingleCellSubset(subsetName = subsetName, rowIndex = rowIndex, colIndex = colIndex)
            obj@subsets[[subsetName]] <- scs
            return(obj)
          }
)

