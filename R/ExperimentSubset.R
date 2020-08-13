#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndices = "numeric",
                                colIndices = "numeric",
                                useAssay = "character",
                                internalAssay = "SingleCellExperiment"
                              )
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(
  subsetName = "subset",
  rowIndices = NULL,
  colIndices = NULL,
  useAssay = "counts",
  internalAssay = SingleCellExperiment::SingleCellExperiment(),
  ...)
{
  .SingleCellSubset(subsetName = subsetName,
                    rowIndices = rowIndices,
                    colIndices = colIndices,
                    useAssay = useAssay,
                    internalAssay = internalAssay)
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
#' @param rows Specify which rows/features to subset from the input object. Either indices or names.
#' @param cols Specify which columns/cells to subset from the input object. Either indices or names.
#' @param useAssay Assay to use against the subset data
#'
#' @export
setGeneric(name = "subsetAssay",
           def = function(object, subsetName, rows, cols, useAssay)
           {
             standardGeneric("subsetAssay")
           }
)

#' @export
setMethod(f = "subsetAssay",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, rows, cols, useAssay)
            {
            if(is.character(rows)){
              rows <- match(rows, rownames(assay(object, useAssay)))
            }
            if(is.character(cols)){
              cols <- match(cols, colnames(assay(object, useAssay)))
            }
              scs <- SingleCellSubset(
                subsetName = subsetName,
                rowIndices = rows,
                colIndices = cols,
                useAssay = useAssay)
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
    if(!i %in% SummarizedExperiment::assayNames(x)){
      i = x@subsets[[subsetName]]@useAssay
    }
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
setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, useAssay = NULL, ..., value) {
  if((nrow(value)!= nrow(x))
     || (ncol(value) != ncol(x))){
    if(is.null(useAssay)){
      saveSubset(
        object = x,
        subsetName = i,
        inputMatrix = value
      )
    }
    else{
      subsetAssay(
        object = x,
        subsetName = i,
        rows = rownames(value),
        cols = colnames(value),
        useAssay = useAssay
      )
    }
  }
  else{
    callNextMethod()
  }
})

#' @export
setGeneric(name = "subsetRowData",
           def = function(object, subsetName)
           {
             standardGeneric("subsetRowData")
           }
)


#' @export
setMethod(f = "subsetRowData",
          signature = c("ExperimentSubset", "character"),
          definition = function(object, subsetName)
          {
            SummarizedExperiment::rowData(object)[object@subsets[[subsetName]]@rowIndices, , drop = F]
          }
)

#' @export
setGeneric(name = "subsetColData",
           def = function(object, subsetName)
           {
             standardGeneric("subsetColData")
           }
)


#' @export
setMethod(f = "subsetColData",
          signature = c("ExperimentSubset", "character"),
          definition = function(object, subsetName)
          {
            SummarizedExperiment::colData(object)[object@subsets[[subsetName]]@colIndices, , drop = F]
          }
)

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

            #analyze this further (needed for SCTK)
            rownames(counts) <- gsub("_", "-", rownames(object))
            colnames(counts) <- gsub("_", "-", colnames(object))

            m <- Matrix::Matrix(
              nrow = nrow(counts),
              ncol = ncol(counts),
              data = 0,
              dimnames = list(
                rownames(counts),
                colnames(counts)),
              sparse = TRUE)

            m[r,c] <- inputMatrix #find a faster copying mechanism

            SummarizedExperiment::assay(object, paste0(subsetName, "_internal")) <- m

             object <- subsetAssay(object, subsetName, r, c, paste0(subsetName, "_internal"))

            object@subsets[[subsetName]]@internalAssay <- SingleCellExperiment::SingleCellExperiment(list(counts = inputMatrix))

            return(object)
          }
)

#' @export
#' @importMethodsFrom SummarizedExperiment rowData
setMethod("rowData", c("ExperimentSubset"), function(x, subsetName = NULL, ...) {
  if(!is.null(subsetName)){
    out <- subsetRowData(
      object = x,
      subsetName = subsetName
    )
  }
  else{
    out <- callNextMethod()
  }
  out
})

#' @export
#' @importMethodsFrom SummarizedExperiment colData
setMethod("colData", c("ExperimentSubset"), function(x, subsetName = NULL, ...) {
  if(!is.null(subsetName)){
    out <- subsetColData(
      object = x,
      subsetName = subsetName
    )
  }
  else{
    out <- callNextMethod()
  }
  out
})
