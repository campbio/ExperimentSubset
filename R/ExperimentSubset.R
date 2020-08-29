setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))
setClassUnion("NullOrNumericOrCharacter", c("NULL", "numeric", "character"))
#add these restrictions for all functions

#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass("SingleCellSubset",
                              slots = representation(
                                subsetName = "character",
                                rowIndices = "NullOrNumeric",
                                colIndices = "NullOrNumeric",
                                parentAssay = "NullOrCharacter",
                                internalAssay = "SingleCellExperiment"
                              )
)

#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(
  subsetName = "subset",
  rowIndices = NULL,
  colIndices = NULL,
  parentAssay = "counts",
  internalAssay = SingleCellExperiment::SingleCellExperiment(),
  ...)
{
  .SingleCellSubset(subsetName = subsetName,
                    rowIndices = rowIndices,
                    colIndices = colIndices,
                    parentAssay = parentAssay,
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
  es <- .ExperimentSubset(se)
  if(!missing(subsets)){
      es <- ExperimentSubset::createSubset(es,
                         subsetName = subsets[[1]],
                         rows = subsets[[2]],
                         cols = subsets[[3]],
                         parentAssay = subsets[[4]])
  }
  es
}

#' @title createSubset
#'
#' @param object \code{ExperimentSubset}, \code{SingleCellExperiment} or \code{SummarizedExperiment} object.
#'
#' @param subsetName Specify the name of the new subset.
#' @param rows Specify which rows/features to subset from the input object. Either indices or names.
#' @param cols Specify which columns/cells to subset from the input object. Either indices or names.
#' @param parentAssay Assay to use against the subset data
#'
#' @export
setGeneric(name = "createSubset",
           def = function(object, subsetName, rows, cols, parentAssay)
           {
             standardGeneric("createSubset")
           }
)

#' @export
setMethod(f = "createSubset",
          signature = c("ExperimentSubset",
                        "character",
                        "NullOrNumericOrCharacter",
                        "NullOrNumericOrCharacter",
                        "NullOrCharacter"),
          definition = function(object,
                                subsetName,
                                rows,
                                cols,
                                parentAssay)
            {
            tempAssay <- "" #better way to use this?
            if(is.null(parentAssay)){
              tempAssay <- "counts" #get first assay instead of counts
            }
            else{
              tempAssay <- parentAssay
            }
            if(is.character(rows)){
              rows <- match(rows, rownames(assay(object, tempAssay)))
            }
            if(is.character(cols)){
              cols <- match(cols, colnames(assay(object, tempAssay)))
            }
              scs <- SingleCellSubset(
                subsetName = subsetName,
                rowIndices = rows,
                colIndices = cols,
                parentAssay = parentAssay,
                internalAssay = SingleCellExperiment::SingleCellExperiment(
                  list(
                    counts = Matrix::Matrix(
                      nrow = length(rows),
                      ncol = length(cols),
                      data = 0,
                      sparse = TRUE))))
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

#' #' @export
#' #' @importMethodsFrom SummarizedExperiment assay
#' setMethod("assay", c("ExperimentSubset", "character"), function(x, i, ...) {
#'   if(i %in% subsetNames(x)){
#'     subsetName = i
#'     if(!i %in% SummarizedExperiment::assayNames(x)){
#'       i = x@subsets[[subsetName]]@parentAssay
#'     }
#'     if(!is.null(x@subsets[[subsetName]]@parentAssay)){
#'       out <- callNextMethod()
#'       out <- out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
#'     }
#'     else{
#'       out <- assay(x@subsets[[subsetName]]@internalAssay, "counts")
#'     }
#'   }
#'   else{
#'     out <- callNextMethod()
#'   }
#'   out
#' })
#'

#' @title assay
#' @export
#' @importMethodsFrom SummarizedExperiment assay
setMethod("assay", c("ExperimentSubset", "character"), function(x, i, ...) {
  #look at main assays
  if(i %in% assayNames(x)){
    out <- callNextMethod()
  }
  #look at subsets
  else if(i %in% subsetNames(x)){
    subsetName <- i
    i <- x@subsets[[subsetName]]@parentAssay
    if(is.null(i)){
      out <- assay(x@subsets[[subsetName]]@internalAssay, "counts")
    }
    else{
      out <- assay(x, i)
      out <- out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
    }
  }
  #look inside subsets
  else{
    for(j in seq(length(x@subsets))){
      print(assayNames(x@subsets[[j]]@internalAssay))
      if(i %in% assayNames(x@subsets[[j]]@internalAssay)){
        out <- assay(x@subsets[[j]]@internalAssay, i)
        colnames(out) <- colnames(x)[x@subsets[[j]]@colIndices]
        rownames(out) <- rownames(x)[x@subsets[[j]]@rowIndices]
      }
    }
  }
  out
})

#' #' @export
#' #' @importMethodsFrom SummarizedExperiment assay<-
#' setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, parentAssay = NULL, newInternalAssay = NULL, ..., value) {
#'   if((nrow(value)!= nrow(x))
#'      || (ncol(value) != ncol(x))){
#'     if(is.null(parentAssay)){
#'       if(is.null(newInternalAssay)){
#'         storeSubset(
#'           object = x,
#'           subsetName = i,
#'           inputMatrix = value
#'         )
#'       }
#'       else{
#'         if((nrow(value)!= nrow(x@subsets[[i]]@internalAssay))
#'            || (ncol(value) != ncol(x@subsets[[i]]@internalAssay))){
#'               createSubset(
#'                 object = x,
#'                 subsetName = newInternalAssay,
#'                 rows = match(rownames(value), rownames(x)),
#'                 cols = match(colnames(value), colnames(x)),
#'                 parentAssay = i
#'               )
#'         }
#'         else{
#'           storeSubset(
#'             object = x,
#'             subsetName = i,
#'             inputMatrix = value,
#'             newInternalAssay = newInternalAssay
#'           )
#'         }
#'       }
#'     }
#'     else{
#'       createSubset(
#'         object = x,
#'         subsetName = i,
#'         rows = rownames(value),
#'         cols = colnames(value),
#'         parentAssay = parentAssay
#'       )
#'     }
#'   }
#'   else{
#'     callNextMethod()
#'   }
#' })

#' @title assay
#' @export
#' @importMethodsFrom SummarizedExperiment assay<-
setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, parentAssay = NULL, newInternalAssay = NULL, ..., value) {
  if((nrow(value)!= nrow(x))
     || (ncol(value) != ncol(x))){
    storeSubset(
                object = x,
                subsetName = i,
                inputMatrix = value,
                newInternalAssay = NULL
              )
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
            out <- SummarizedExperiment::rowData(object)[object@subsets[[subsetName]]@rowIndices, , drop = F]
            if(!is.null(object@subsets[[subsetName]]@internalAssay)){
              if(ncol(rowData(object@subsets[[subsetName]]@internalAssay)) > 0){
                out <- cbind(out, rowData(object@subsets[[subsetName]]@internalAssay))
              }
            }
            out
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
            out <- SummarizedExperiment::colData(object)[object@subsets[[subsetName]]@colIndices, , drop = F]
            if(!is.null(object@subsets[[subsetName]]@internalAssay)){
              if(ncol(colData(object@subsets[[subsetName]]@internalAssay)) > 0){
                out <- cbind(out, colData(object@subsets[[subsetName]]@internalAssay))
              }
            }
            out
          }
)

#' @export
setGeneric(name = "storeSubset",
           def = function(object, subsetName, inputMatrix, newInternalAssay)
           {
             standardGeneric("storeSubset")
           }
)

#' @export
setMethod(f = "storeSubset",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, inputMatrix, newInternalAssay = NULL)
          {
            if(is.null(newInternalAssay)){
              object <- createSubset(
                object,
                subsetName,
                rownames(inputMatrix),
                colnames(inputMatrix),
                parentAssay = NULL)

              object@subsets[[subsetName]]@internalAssay <- SingleCellExperiment::SingleCellExperiment(list(counts = inputMatrix))

            }
            else{
              assay(object@subsets[[subsetName]]@internalAssay, newInternalAssay) <- inputMatrix
            }

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


#' @export
#' @importMethodsFrom SummarizedExperiment rowData<-
setReplaceMethod("rowData", c("ExperimentSubset"), function(x, ..., subsetName, value) { #test if this needs DataFrame too
  tempValue <- rowData(x)
  rowData(x@subsets[[subsetName]]@internalAssay) <- value #DataFrame(test = c("a", "b", "c","d", "e"))
  value <- tempValue
  callNextMethod()
})

#' @export
#' @importMethodsFrom SummarizedExperiment colData<-
setReplaceMethod("colData", c("ExperimentSubset" , "DataFrame"), function(x, ..., subsetName, value) {
  tempValue <- colData(x)
  colData(x@subsets[[subsetName]]@internalAssay) <- value #DataFrame(test = c("a", "b", "c","d", "e"))
  value <- tempValue
  callNextMethod()
})
