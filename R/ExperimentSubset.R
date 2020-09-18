setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("MissingOrNullOrCharacter", c("missing","NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))
setClassUnion("NullOrMissingOrNumericOrCharacter", c("NULL", "missing", "numeric", "character"))


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
  if(grepl("\\s+", subsetName)){
    subsetName <- gsub("\\s", "", subsetName)
    warning("Removing spaces from subsetName argument.")
  }
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
                          prototype = list(
                            subsets = list()
                          ),
                          contains = "SingleCellExperiment"
)

#' @export
#' @import BiocStyle
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
ExperimentSubset <- function(
  object,
  ...,
  subset = list(subsetName = NA, rows = NA, cols = NA, parentAssay = NA))
{
  if(!missing(object)){
    if(inherits(object, "SummarizedExperiment")){
      object <- as(object, "SingleCellExperiment")
      es <- .ExperimentSubset(object)
    }
  }
  else{
    se <- SingleCellExperiment::SingleCellExperiment(...)
    es <- .ExperimentSubset(se)
  }
  if(!anyNA(subset)){
      es <- ExperimentSubset::createSubset(es,
                         subsetName = subset$subsetName,
                         rows = subset$rows,
                         cols = subset$cols,
                         parentAssay = subset$parentAssay)
  }
  es
}

#' @title createSubset
#' @description Create a subset from an already available \code{assay} in the
#' input \code{ExperimentSubset} object by specifying the rows and columns
#' to include in the subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Specify the name of the subset to create.
#' @param rows Specify the rows to include in this subset.
#' If \code{missing} or \code{NULL}, all rows are included in the subset.
#' Values can be \code{numeric} or \code{character}.
#' Default \code{NULL}.
#' @param cols Specify the columns to include in this subset.
#' If \code{missing} or \cude{NULL}, all columns are included in the subset.
#' Values can be \code{numeric} or \code{character}.
#' Default \code{NULL}.
#' @param parentAssay Specify the parent \code{assay} of the subset. This parent
#' \code{assay} must already be available in the \code{ExperimentSubset} object.
#' If \code{NULL}, the first available main \code{assay} will be marked as parent.
#' Default \code{NULL}.
#' @return An \code{ExperimentSubset} object that now contains the newly created subset.
#'
#' @export
setGeneric(name = "createSubset",
           def = function(object, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)
           {
             standardGeneric("createSubset")
           }
)

#' @export
setMethod(f = "createSubset",
          signature = c("ExperimentSubset",
                        "character",
                        "NullOrMissingOrNumericOrCharacter",
                        "NullOrMissingOrNumericOrCharacter",
                        "MissingOrNullOrCharacter"),
          definition = function(object,
                                subsetName,
                                rows,
                                cols,
                                parentAssay)
            {
            tempAssay <- ""
            if(is.null(parentAssay)){
              tempAssay <- SummarizedExperiment::assayNames(object)[1]
              parentAssay <- tempAssay
            }
            else{
              if(parentAssay %in% SummarizedExperiment::assayNames(object)
                 || parentAssay %in% subsetAssayNames(object)){
                tempAssay <- parentAssay
              }
              else{
                stop("Input parentAssay does not exist.")
              }
            }
            if(is.character(rows)){
              rows <- match(rows, rownames(assay(object, tempAssay)))
            }
            if(is.character(cols)){
              cols <- match(cols, colnames(assay(object, tempAssay)))
            }
            if(is.null(rows)){
              rows <- seq(1, dim(assay(object, tempAssay))[1])
            }
            if(is.null(cols)){
              cols <- seq(1, dim(assay(object, tempAssay))[2])
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
                      sparse = TRUE
                      )
                    )
                  )
                )

              #Check if NAs introduced in the subset
              tryCatch({
                na.fail(scs@rowIndices)
                na.fail(scs@colIndices)
              }, error = function(e){
                stop("NAs introduced in input rows or columns.
                     Some or all indicated rows or columns not found in specified parent.")
              })

              #Remove counts assay from internal SCE object of the subset to save memory
              assay(scs@internalAssay, "counts") <- NULL

              object@subsets[[subsetName]] <- scs
              return(object)
          }
)

#' @title subsetNames
#' @description Retrieves the names of the available subsets in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @return A \code{vector} of subset names.
#'
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
setGeneric(name = "altExps",
           def = function(x, withColData = FALSE, subsetName)
           {
             standardGeneric("altExps")
           }
)

#' @export
setMethod(f = "altExps",
          signature = "ANY",
          definition = function(x, withColData, subsetName)
          {
            if(!missing(subsetName)){
              if(is.null(x@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              SingleCellExperiment::altExps(x@subsets[[subsetName]]@internalAssay, withColData = withColData)
            }
            else{
              SingleCellExperiment::altExps(x, withColData = withColData)
            }
          }
)

#' @export
setGeneric(name = "altExp",
           def = function(x, e, withColData = FALSE, subsetName)
           {
             standardGeneric("altExp")
           }
)

#' @export
setMethod(f = "altExp",
          signature = "ANY",
          definition = function(x, e, withColData, subsetName)
          {
            if(!missing(subsetName)){
              if(is.null(x@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              if(missing(e)){
                SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, withColData = withColData)
              }
              else{
                SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, e, withColData = withColData)
              }
            }
            else{
              if(missing(e)){
                SingleCellExperiment::altExp(x, withColData = withColData)
              }
              else{
                SingleCellExperiment::altExp(x, e, withColData = withColData)
              }
            }
          }
)

#' @export
setGeneric(name = "altExpNames",
           def = function(x, subsetName)
           {
             standardGeneric("altExpNames")
           }
)

#' @export
setMethod(f = "altExpNames",
          signature = "ANY",
          definition = function(x, subsetName)
          {
            if(!missing(subsetName)){
              if(is.null(x@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              SingleCellExperiment::altExpNames(x@subsets[[subsetName]]@internalAssay)
            }
            else{
              SingleCellExperiment::altExpNames(x)
            }
          }
)

#' @export
setGeneric(name = "reducedDimNames",
           def = function(x, subsetName)
           {
             standardGeneric("reducedDimNames")
           }
)

#' @export
setMethod(f = "reducedDimNames",
          signature = "ANY",
          definition = function(x, subsetName)
          {
            if(!missing(subsetName)){
              if(is.null(x@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              SingleCellExperiment::reducedDimNames(x@subsets[[subsetName]]@internalAssay)
            }
            else{
              SingleCellExperiment::reducedDimNames(x)
            }
          }
)

#' @export
setGeneric(name = "altExpNames<-",
           def = function(x, subsetName, value)
           {
             standardGeneric("altExpNames<-")
           }
)

#' @export
setReplaceMethod(f = "altExpNames",
                 signature = "ANY",
                 definition = function(x, subsetName, value)
                 {
                   if(!missing(subsetName)){
                     if(is.null(x@subsets[[subsetName]])){
                       stop(paste(subsetName, "does not exist in the subsets slot of the object."))
                     }
                     SingleCellExperiment::altExpNames(x@subsets[[subsetName]]@internalAssay) <- value
                   }
                   else{
                     SingleCellExperiment::altExpNames(x) <- value
                   }
                   x
                 }
)

#' @export
setGeneric(name = "reducedDimNames<-",
           def = function(x, subsetName, value)
           {
             standardGeneric("reducedDimNames<-")
           }
)

#' @export
setReplaceMethod(f = "reducedDimNames",
                 signature = "ANY",
                 definition = function(x, subsetName, value)
                 {
                   if(!missing(subsetName)){
                     if(is.null(x@subsets[[subsetName]])){
                       stop(paste(subsetName, "does not exist in the subsets slot of the object."))
                     }
                     SingleCellExperiment::reducedDimNames(x@subsets[[subsetName]]@internalAssay) <- value
                   }
                   else{
                     SingleCellExperiment::reducedDimNames(x) <- value
                   }
                   x
                 }
)

#' @export
setGeneric(name = "altExp<-",
           def = function(x, e, withColData = FALSE, subsetName, value)
           {
             standardGeneric("altExp<-")
           }
)

#' @export
setReplaceMethod(f = "altExp",
          signature = "ANY",
          definition = function(x, e, withColData, subsetName, value)
          {
            if(!missing(subsetName)){
              if(is.null(x@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              if(missing(e)){
                SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, withColData = withColData) <- value
              }
              else{
                SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, e, withColData = withColData) <- value
              }
            }
            else{
              if(missing(e)){
                SingleCellExperiment::altExp(x, withColData = withColData) <- value
              }
              else{
                SingleCellExperiment::altExp(x, e, withColData = withColData) <- value
              }
            }
            x
          }
)

#' @export
setGeneric(name = "altExps<-",
           def = function(x, subsetName, value)
           {
             standardGeneric("altExps<-")
           }
)

#' @export
setReplaceMethod(f = "altExps",
                 signature = "ANY",
                 definition = function(x, subsetName, value)
                 {
                   if(!missing(subsetName)){
                     if(is.null(x@subsets[[subsetName]])){
                       stop(paste(subsetName, "does not exist in the subsets slot of the object."))
                     }
                     SingleCellExperiment::altExps(x@subsets[[subsetName]]@internalAssay) <- value
                   }
                   else{
                     SingleCellExperiment::altExps(x) <- value
                   }
                   x
                 }
)


#' @export
setGeneric(name = "metadata",
           def = function(object, subsetName)
           {
             standardGeneric("metadata")
           }
)

#' @export
setMethod(f = "metadata",
          signature = "ANY",
          definition = function(object, subsetName)
          {
            if(!missing(subsetName)){
              if(is.null(object@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              S4Vectors::metadata(object@subsets[[subsetName]]@internalAssay)
            }
            else{
              S4Vectors::metadata(object)
            }
          }
)

#' @title subsetDim
#' @description Retrieves the dimensions of the specified subset in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to retrieve the dimensions from.
#' @return A \code{vector} containing the dimensions of the specified
#' subset i.e. the number of rows and the number of columns in the subset.
#' @export
setGeneric(name = "subsetDim",
           def = function(object, subsetName)
           {
             standardGeneric("subsetDim")
           }
)

#' @export
setMethod(f = "subsetDim",
          signature = c("ExperimentSubset", "character"),
          definition = function(object, subsetName)
          {
            dim(object@subsets[[subsetName]]@internalAssay)
          }
)

#' @export
setGeneric(name = "metadata<-",
           def = function(object, subsetName, value)
           {
             standardGeneric("metadata<-")
           }
)

#' @export
setReplaceMethod(f = "metadata",
          signature = "ANY",
          definition = function(object, subsetName, value)
          {
            if(!missing(subsetName)){
              if(is.null(object@subsets[[subsetName]])){
                stop(paste(subsetName, "does not exist in the subsets slot of the object."))
              }
              S4Vectors::metadata(object@subsets[[subsetName]]@internalAssay) <- value
            }
            else{
              S4Vectors::metadata(object) <- value
            }
            object
          }
)

#' @title subsetCount
#' @description Get the total count of the available subsets in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @return A \code{numeric} value representing the total count of the subsets.
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

#' @title subsetAssayCount
#' @description Get the count of the total available subsets and the subset assays inside
#' these subsets in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @return A \code{numeric} value representing the sum of the subset count and subset assay count.
#' @export
setGeneric(name = "subsetAssayCount",
           def = function(object)
           {
             standardGeneric("subsetAssayCount")
           }
)

#' @param ExperimentSubset
#'
#' @export
setMethod(f = "subsetAssayCount",
          signature = "ExperimentSubset",
          definition = function(object)
          {
            return(length(subsetAssayNames(object)))
          }
)

#' @title showSubsetLink
#' @description The function displays the content of an \code{ExperimentSubset} object
#' including all available main assays, all subsets and the subset assays inside these subsets.
#' This function also depicts how and in what order the subsets in the object are linked with
#' their parents. Moreover, all supplementary data inside the subsets such as \code{reducedDims}
#' and \code{altExps} are also displayed against each subset entry.
#' @param object Input \code{ExperimentSubset} object.
#' @return Prints all the available subset information against the input \code{ExperimentSubset}
#' object.
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
            cat("Main assay(s):\n", assayNames(object),"\n\n")
            cat("Subset(s):\n")
            if(!is.null(subsetNames(object))){
              Name <- list()
              Dimensions <- list()
              Parent <- list()
              Assays <- list()
              Metadata <- list()
              ReducedDims <- list()
              AltExperiments <- list()

              for(i in seq(length(subsetNames(object)))){
                parent <- subsetParent(object, subsetAssayNames(object)[i])
                Name[[i]] <- subsetNames(object)[i]
                Parent[[i]] <- paste(unlist(parent), collapse = ' -> ')
                Assays[[i]] <- assayNames(object@subsets[[i]]@internalAssay)
                Dimensions[[i]] <- paste(unlist(subsetDim(object, subsetNames(object)[i])), collapse = ', ')
                ReducedDims[[i]] <- paste(unlist(reducedDimNames(object, subsetNames(object)[i])), collapse = ", ")
                AltExperiments[[i]] <- paste(unlist(altExpNames(object, subsetName = subsetNames(object)[i])), collapse = ", ")
              }

              Assays[lengths(Assays) == 0] <- ""
              ReducedDims[lengths(ReducedDims) == 0] <- ""
              AltExperiments[lengths(AltExperiments) == 0] <- ""

              df <- data.frame(
                Name = as.character(Name),
                Dim = as.character(Dimensions),
                Parent = as.character(Parent)
                )

              if(length(which(as.character(Assays) == "")) != subsetCount(object)){
                df <- cbind(df, Assays = as.character(Assays))
              }

              if(length(which(as.character(AltExperiments) == "")) != subsetCount(object)){
                df <- cbind(df, AltExperiments = as.character(AltExperiments))
              }

              if(length(which(as.character(ReducedDims) == "")) != subsetCount(object)){
                df <- cbind(df, ReducedDims = as.character(ReducedDims))
              }

              print(df)
            }
            else{
              cat("NULL\n")
            }
          }
)

#' @title subsetParent
#' @description Retrieves a complete subset to parent link from a specified subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Specify the name of the subset against which the subset to parent link
#' should be retrieved.
#' @return A \code{list} containing the parent link of the subset.
#' @export
setGeneric(name = "subsetParent",
           def = function(object, subsetName)
           {
             standardGeneric("subsetParent")
           }
)

#' @export
setMethod(f = "subsetParent",
          signature = "ANY",
          definition = function(object, subsetName)
          {
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
              }
              if(parent %in% assayNames(object)){
                parentList <- c(parentList, parent)
                break
              }
            }
            parentList[[1]] <- NULL
            return(parentList)
          }
)


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

#' @title subsetAssayNames
#' @description Retrieves the names of all the subsets as well as the subset assays.
#' @param object Input \code{ExperimentSubset} object.
#' @return A \code{vector} containing the names of the subsets and the subset assays
#' available in the input \code{ExperimentSubset} object.
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
      out <- SummarizedExperiment::assay(x, i)
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
setReplaceMethod("assay", c("ExperimentSubset", "character"), function(x, i, subsetAssayName = NULL, ..., value) {
  if((nrow(value)!= nrow(x))
     || (ncol(value) != ncol(x))){
    storeSubset(
                object = x,
                subsetName = i,
                inputMatrix = value,
                subsetAssayName = subsetAssayName
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
              stop("Neither a subset nor a subsetAssay.")
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
              stop("Neither a subset nor a subsetAssay.")
            }
            return(out)
          }
)


#' @title storeSubset
#' @description Store a new subset \code{assay} inside a specified subset in the input
#' \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object
#' @param subsetName Specify the name of the existing subset inside which the new subset \code{assay}
#' should be stored.
#' @param inputMatrix The input subset \code{assay}.
#' @param subsetAssayName Specify the name of the new \code{assay} against the \code{inputMatrix} parameter.
#' If \code{NULL}, a new subset is created internally using the \code{createSubset} function.
#' Default \code{NULL}.
#' @return Updated \code{ExperimentSubset} object with the new \code{assay} stored inside the specified subset.
#'
#' @export
setGeneric(name = "storeSubset",
           def = function(object, subsetName, inputMatrix, subsetAssayName)
           {
             standardGeneric("storeSubset")
           }
)

#' @export
setMethod(f = "storeSubset",
          signature = "ExperimentSubset",
          definition = function(object, subsetName, inputMatrix, subsetAssayName = NULL)
          {
            if(is.null(object@subsets[[subsetName]])){
              stop(paste(subsetName, "does not exist in the subsets slot of the object."))
            }

            if(!all(dim(object@subsets[[subsetName]]@internalAssay) == dim(inputMatrix))){
              stop("Dimensions of the inputMatrix not equal to the subset. You need to create a new subset with createSubset() function.")
            }

            if(is.null(subsetAssayName)){
              object <- createSubset(
                object,
                subsetName,
                rownames(inputMatrix),
                colnames(inputMatrix),
                parentAssay = NULL)

              object@subsets[[subsetName]]@internalAssay <- SingleCellExperiment::SingleCellExperiment(
                list(counts = inputMatrix))

            }
            else{
              assay(object@subsets[[subsetName]]@internalAssay, subsetAssayName) <- inputMatrix
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
setMethod("reducedDim", "ANY", function(object, type, withDimnames, subsetName) {
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
setGeneric(name = "reducedDims",
           def = function(object, withDimnames, subsetName)
           {
             standardGeneric("reducedDims")
           }
)

#' @export
setMethod("reducedDims", "ANY", function(object, withDimnames, subsetName) {
  if(missing(withDimnames)){
    withDimnames = TRUE
  }
  if(!missing(subsetName)){
    out <- SingleCellExperiment::reducedDims(object@subsets[[subsetName]]@internalAssay, withDimnames)
  }
  else{
    out <- SingleCellExperiment::reducedDims(object, withDimnames)
  }
  out
})

#' @export
setGeneric(name = "reducedDim<-",
           def = function(object, type, subsetName, value)
           {
             standardGeneric("reducedDim<-")
           }
)

#' @export
setReplaceMethod("reducedDim", "ANY", function(object, type, subsetName, value) {
  if(!missing(subsetName)){
    SingleCellExperiment::reducedDim(object@subsets[[subsetName]]@internalAssay, type) <- value
  }
  else{
    SingleCellExperiment::reducedDim(object, type) <- value
  }
  return(object)
})

#' @export
setGeneric(name = "reducedDims<-",
           def = function(object, subsetName, value)
           {
             standardGeneric("reducedDims<-")
           }
)

#' @export
setReplaceMethod("reducedDims", "ANY", function(object, subsetName, value) {
  if(!missing(subsetName)){
    SingleCellExperiment::reducedDims(object@subsets[[subsetName]]@internalAssay) <- value
  }
  else{
    SingleCellExperiment::reducedDims(object) <- value
  }
  return(object)
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
