#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndices = "numeric",
                                colIndices = "numeric",
                                useAssay = "character"
                              )
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(
  subsetName = "subset",
  rowIndices = NULL,
  colIndices = NULL,
  useAssay = "counts",
  ...)
{
  .SingleCellSubset(subsetName = subsetName,
                    rowIndices = rowIndices,
                    colIndices = colIndices,
                    useAssay = useAssay)
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
    i = paste0(subsetName, "_internal")
    out <- callNextMethod()
    out <- out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
  }
  else{
    out <- callNextMethod()
  }
  out
})

#' @export
#' @importMethodsFrom SummarizedExperiment assay<-
setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, ..., value) {
  if((nrow(value)!= nrow(x))
     || (ncol(value) != ncol(x))){
    saveSubset(
      object = x,
      subsetName = i,
      inputMatrix = value
    )
  }
  else{
    callNextMethod()
  }
})


#' @export
setGeneric(name = "saveSubset",
           def = function(object, subsetName, inputMatrix)
           {
             standardGeneric("saveSubset")
           }
)

#' @export
setMethod(f = "saveSubset",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, inputMatrix)
          {
            r <- rownames(inputMatrix)
            c <- colnames(inputMatrix)
            counts <- assay(object, "counts")

            m <- Matrix::Matrix(
              nrow = nrow(counts),
              ncol = ncol(counts),
              data = 0,
              dimnames = list(
                rownames(counts),
                colnames(counts)),
              sparse = TRUE)

            m[r,c] <- inputMatrix

            SummarizedExperiment::assay(object, paste0(subsetName, "_internal")) <- m

            object <- subsetAssay(object, subsetName, r, c) #add parameter here to store useAssay name in line below

            object@subsets[[subsetName]]@useAssay <- paste0(subsetName, "_internal") #omit this line

            return(object)
          }
)
