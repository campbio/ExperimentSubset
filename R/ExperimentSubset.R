setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))
setClassUnion("NullOrNumericOrCharacter", c("NULL", "numeric", "character"))
setClassUnion("MissingOrNumericOrCharacter", c("missing", "numeric", "character"))
setClassUnion("MissingOrLogical", c("missing", "logical"))
setClassUnion("MissingOrCharacter", c("missing", "character"))
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
              assay(scs@internalAssay, "counts") <- NULL #a better way to do this
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
setGeneric(name = "subsetCount",
           def = function(object)
           {
             standardGeneric("subsetCount")
           }
)

#' @export
setMethod(f = "subsetCount",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            return(length(subsetNames(object)))
          }
)

#' @export
setGeneric(name = "subsetAssayCount",
           def = function(object)
           {
             standardGeneric("subsetAssayCount")
           }
)

#' @export
setMethod(f = "subsetAssayCount",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            return(length(subsetAssayNames(object)))
          }
)

#' @export
setGeneric(name = "showSubsetLink",
           def = function(object)
           {
             standardGeneric("showSubsetLink")
           }
)

#' @export
setMethod(f = "showSubsetLink",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            # cat("Main assay(s):\n", assayNames(es),"\n")
            # cat("Subset(s):\n")
            #   for(i in seq(length(subsetNames(object)))){
            #     cat(" ", sep = "")
            #     parent <- getFirstParent(object, subsetAssayNames(object)[i])
            #     cat(i,") ", subsetNames(object)[i], ", parent =", sep = "")
            #     for(j in seq(length(parent))){
            #       cat(" ", parent[[j]], sep = "")
            #     }
            #     cat(", assay(s) = ", assayNames(object@subsets[[i]]@internalAssay))
            #     cat("\n")
            #   }

            cat("Main assay(s):\n", assayNames(es),"\n\n")
            cat("Subset(s):")
            Name <- list()
            Parent <- list()
            Assays <- list()
            for(i in seq(length(subsetNames(object)))){
              parent <- getFirstParent(object, subsetAssayNames(object)[i])
              Name[[i]] <- subsetNames(object)[i]
              Parent[[i]] <- paste(unlist(parent), collapse = ' ')
              Assays[[i]] <- assayNames(object@subsets[[i]]@internalAssay)
            }

            Assays[lengths(Assays) == 0] <- ""

            df <- data.frame(Name = as.character(Name), Parent = as.character(Parent), Assays = as.character(Assays))
            knitr::kable(df)
          }
)

#' @export
getFirstParent <- function(object, subsetName){
  parentList <- list()
  parent <- subsetName
  while(TRUE){
    parentList <- c(parentList, parent)
    if(!is.null(object@subsets[[parent]])){
      parent <- object@subsets[[parent]]@parentAssay
    }
    else{
      for(i in seq(subsetCount(object))){
        if(parent %in% assayNames(object@subsets[[i]]@internalAssay)){
          parent <- object@subsets[[i]]@subsetName
        }
      }
      parentList <- c(parentList, parent)
      parent <- object@subsets[[parent]]@parentAssay
      #print(parent)
    }
    if(parent %in% assayNames(object)){
      parentList <- c(parentList, parent)
      break
    }
  }
  parentList[[1]] <- NULL
  return(parentList)
}


#' @export
setGeneric(name = "rownames",
           def = function(object, ...)
           {
             standardGeneric("rownames")
           }
)

#' @export
setMethod(f = "rownames",
          signature = "ANY",
          definition = function(object, subsetName)
          {
            if(missing(subsetName)){
              BiocGenerics::rownames(object)
            }
            else{
              if(subsetName %in% subsetNames(object)){
                BiocGenerics::rownames(object)[object@subsets[[subsetName]]@rowIndices]
              }
              else if(subsetName %in% subsetAssayNames(object)){
                subsetName <- .getParentAssayName(object, subsetName)
                BiocGenerics::rownames(object)[object@subsets[[subsetName]]@rowIndices]
              }
            }
          }
)

#' @export
setGeneric(name = "colnames",
           def = function(object, ...)
           {
             standardGeneric("colnames")
           }
)

#' @export
setMethod(f = "colnames",
          signature = "ANY",
          definition = function(object, subsetName)
          {
            if(missing(subsetName)){
              BiocGenerics::colnames(object)
            }
            else{
              if(subsetName %in% subsetNames(object)){
                BiocGenerics::colnames(object)[object@subsets[[subsetName]]@colIndices]
              }
              else if(subsetName %in% subsetAssayNames(object)){
                subsetName <- .getParentAssayName(object, subsetName)
                BiocGenerics::colnames(object)[object@subsets[[subsetName]]@colIndices]
              }
            }
          }
)

#' @export
setGeneric(name = "subsetAssayNames",
           def = function(object)
           {
             standardGeneric("subsetAssayNames")
           }
)

#' @export
setMethod(f = "subsetAssayNames",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            tempNames <- names(object@subsets)
            if(length(object@subsets)>0){
              for(i in seq(length(object@subsets))){
                tempNames <- c(tempNames, assayNames(object@subsets[[i]]@internalAssay))
              }
            }
            return(tempNames)
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
                  sep=""
                  )
              cat(
                  subsetNames(object)
                  )
              cat(
                "\nsubsetAssays(", length(subsetAssayNames(object)), "): ",
                sep = ""
              )
              cat(
                subsetAssayNames(object)
              )
          }
)


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
      if(i %in% assayNames(x@subsets[[j]]@internalAssay)){
        out <- assay(x@subsets[[j]]@internalAssay, i)
      }
    }
  }
  out
})


#' @title assay
#' @export
#' @importMethodsFrom SummarizedExperiment assay<-
setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, newInternalAssay = NULL, ..., value) {
  if((nrow(value)!= nrow(x))
     || (ncol(value) != ncol(x))){
    storeSubset(
                object = x,
                subsetName = i,
                inputMatrix = value,
                newInternalAssay = newInternalAssay
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
            if(subsetName %in% subsetNames(object)){
              #is a subset
              out <- SummarizedExperiment::rowData(object)[object@subsets[[subsetName]]@rowIndices, , drop = F]
              out <- cbind(out, rowData(object@subsets[[subsetName]]@internalAssay))
            }
            else if(subsetName %in% subsetAssayNames(object)){
              #is a subset assay
              subsetName <- .getParentAssayName(object, subsetName)
              out <- SummarizedExperiment::rowData(object)[object@subsets[[subsetName]]@rowIndices, , drop = F]
              out <- cbind(out, rowData(object@subsets[[subsetName]]@internalAssay))
            }
            else{
              #neither a subset nor a subset assay
              print("neither a subset nor a subset assay")
            }
            return(out)
          }
)


#' @export
setGeneric(name = "subsetColData",
           def = function(object, subsetName)
           {
             standardGeneric("subsetColData")
           }
)

.getParentAssayName <- function(object, childAssayName){
  for(i in seq(length(object@subsets))){
    if(childAssayName %in% assayNames(object@subsets[[i]]@internalAssay)){
      return(object@subsets[[i]]@subsetName)
    }
  }
}

#' @export
setMethod(f = "subsetColData",
          signature = c("ExperimentSubset", "character"),
          definition = function(object, subsetName)
          {
            if(subsetName %in% subsetNames(object)){
              #is a subset
              out <- SummarizedExperiment::colData(object)[object@subsets[[subsetName]]@colIndices, , drop = F]
              out <- cbind(out, colData(object@subsets[[subsetName]]@internalAssay))
            }
            else if(subsetName %in% subsetAssayNames(object)){
              #is a subset assay
              subsetName <- .getParentAssayName(object, subsetName)
              out <- SummarizedExperiment::colData(object)[object@subsets[[subsetName]]@colIndices, , drop = F]
              out <- cbind(out, colData(object@subsets[[subsetName]]@internalAssay))
            }
            else{
              #neither a subset nor a subset assay
              print("neither a subset nor a subset assay")
            }
            return(out)
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
              rownames(object@subsets[[subsetName]]@internalAssay) <- rownames(inputMatrix)
              colnames(object@subsets[[subsetName]]@internalAssay) <- colnames(inputMatrix)
            }

            return(object)
          }
)

#' @export
setGeneric(name = "reducedDim",
           def = function(object, type, withDimnames, subsetName)
           {
             standardGeneric("reducedDim")
           }
)

#' @export
setMethod("reducedDim", c("ExperimentSubset", "MissingOrNumericOrCharacter", "MissingOrLogical", "MissingOrCharacter"), function(object, type, withDimnames, subsetName) {
  if(missing(withDimnames)){
    withDimnames = TRUE
  }
  if(!missing(subsetName)){
    out <- SingleCellExperiment::reducedDim(object@subsets[[subsetName]]@internalAssay, type, withDimnames)
  }
  else{
    out <- SingleCellExperiment::reducedDim(object, type, withDimnames)
  }
  out
})



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
  tempValue <- NULL
  if(!missing(subsetName)){
    tempValue <- rowData(x)
    rowData(x@subsets[[subsetName]]@internalAssay) <- value #DataFrame(test = c("a", "b", "c","d", "e"))
  }
  else{
    tempValue <- value
  }
  value <- tempValue
  callNextMethod()
})

#' @export
#' @importMethodsFrom SummarizedExperiment colData<-
setReplaceMethod("colData", c("ExperimentSubset" , "DataFrame"), function(x, ..., subsetName, value) {
  tempValue <- NULL
  if(!missing(subsetName)){
    tempValue <- colData(x)
    colData(x@subsets[[subsetName]]@internalAssay) <- value #DataFrame(test = c("a", "b", "c","d", "e"))
  }
  else{
    tempValue <- value
  }
  value <- tempValue
  callNextMethod()
})
