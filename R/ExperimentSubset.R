
#' @title AssaySubset constructor
#' @description Constructor for creating a experiment object internally by the
#'   \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset.
#' @param rowIndices \code{vector("numeric")} Indices of the rows to include in
#'   the subset.
#' @param colIndices \code{vector("numeric")} Indices of the columns to include
#'   in the subset.
#' @param parentAssay \code{character(1)} Name of the parent of this subset.
#' @param internalAssay \code{SummarizedExperiment} An internal object to store
#'   additional subset data.
#' @return A \code{AssaySubset} object.
AssaySubset <- function(subsetName = "subset",
                        rowIndices = NULL,
                        colIndices = NULL,
                        parentAssay = "counts",
                        internalAssay = NULL)
{
  if (grepl("\\s+", subsetName)) {
    subsetName <- gsub("\\s", "", subsetName)
    warning("Removing spaces from the specified subsetName.")
  }
  
  x <- .AssaySubset(
    subsetName = subsetName,
    rowIndices = rowIndices,
    colIndices = colIndices,
    parentAssay = parentAssay,
    internalAssay = internalAssay
  )
  
  return(x)
}

#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the \code{ExperimentSubset} object, either through manually specifying the \code{assays}, \code{rowData}, \code{colData} or directly by passing either a \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or objects inherited by these classes. A subset can also be directly created by pasing a named \code{list} to the \code{subset} parameter. This named \code{list} should have parameter values named as \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @param x A \code{SingleCellExperiment} or \code{SummarizedExperiment} object if direct conversion is required.
#' @param subset A named \code{list} if a subset should be created from within the constructor. Named parameters in this list should be \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @return A \code{ExperimentSubset} object.
#' @export
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
ExperimentSubset <- function(x,
                             subset = list(
                               subsetName = NA,
                               rows = NA,
                               cols = NA,
                               parentAssay = NA
                             ))
{
  if (!is.list(x)) {
    if(is.null(assayNames(x))
       || "" %in% assayNames(x))
      stop("Experiment objects with unnamed assays are not supported. Use assayNames<- setter method to set assay names before creating ES object.")
    es <- as(x, paste0("Subset", class(x)))
  }
  else{
    sce <- SingleCellExperiment::SingleCellExperiment(x)
    es <- .SubsetSingleCellExperiment(sce)
  }
  
  if (!anyNA(subset)) {
    es <- createSubset(
      es,
      subsetName = subset$subsetName,
      rows = subset$rows,
      cols = subset$cols,
      parentAssay = subset$parentAssay
    )
  }
  es
}

#' @importMethodsFrom S4Vectors metadata
.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]]))
}

.assay <- function(x, i, ...){
  out <- NULL
  #look at main assays
  if (i %in% assayNames(x)) {
    out <- callNextMethod()
  }
  #look at subsets
  else if (i %in% subsetNames(x)) {
    subsetName <- i
    i <- .parentAssay(.subsets(x)[[subsetName]]) 
    if (is.null(i)) {
      out <- assay(x = .internalAssay(.subsets(x)[[subsetName]]) , i = "counts", ...)
    }
    else{
      out <- assay(x = x, i = i, ...)
      out <- out[.rowIndices(.subsets(x)[[subsetName]]), .colIndices(.subsets(x)[[subsetName]])]
    }
  }
  #look inside subsets
  else{
    for (j in seq(length(x@subsets))) {
      if (i %in% assayNames(.internalAssay(.subsets(x)[[j]]))) {
        out <- assay(.internalAssay(.subsets(x)[[j]]), withDimnames = FALSE,  i)
      }
    }
  }
  out
}

.createSubset <- function(x,
                          subsetName,
                          rows,
                          cols,
                          parentAssay){
  
  .subsetParamsValidity(x,
                 subsetName,
                 rows,
                 cols,
                 parentAssay)
  
  tempAssay <- ""
  if (is.null(parentAssay)) {
    tempAssay <- assayNames(x)[1]
    parentAssay <- tempAssay
  }
  else{
    test <- parentAssay %in% assayNames(x) || 
      parentAssay %in% subsetAssayNames(x)
    if (test) {
      tempAssay <- parentAssay
    }
    else{
      stop("Input parentAssay does not exist.")
    }
  }

  rows <- .rows(x, rows, tempAssay)
  cols <- .cols(x, cols, tempAssay)
  
  a <- list(Matrix(
    nrow = length(rows),
    ncol = length(cols),
    data = 0,
    sparse = TRUE))
  names(a) <- "temp"
  
  internalAssay <- SingleCellExperiment(assays = a)
  internalAssay <- as(internalAssay, gsub("Subset", "", class(x)))
  
  scs <- AssaySubset(
    subsetName = subsetName,
    rowIndices = rows,
    colIndices = cols,
    parentAssay = parentAssay,
    internalAssay = internalAssay
  )
  
  assay(.internalAssay(scs),
        withDimnames = FALSE, "temp") <- NULL
  
  .subsets(x)[[subsetName]] <- scs
  return(x)
}

.subsetAssayNames <- function(x){
  tempNames <- names(.subsets(x))
  if (length(.subsets(x)) > 0) {
    for (i in seq(length(.subsets(x)))) {
      tempNames <-
        c(
          tempNames,
          assayNames(.internalAssay(.subsets(x)[[i]]))
        )
    }
  }
  return(tempNames)
}

.show <- function(object){
  callNextMethod()
  cat("subsets(", length(subsetNames(object)), "): ",
      sep = "")
  cat(subsetNames(object))
  cat("\nsubsetAssays(", length(subsetAssayNames(object)), "): ",
      sep = "")
  cat(subsetAssayNames(object))
}

.setSubsetAssay <- function(x,
                         subsetName,
                         inputMatrix,
                         subsetAssayName = NULL){

  if (!subsetName %in% subsetNames(x))
    stop(subsetName, " does not exist in the subsets slot of the object. 
         You need to create a new subset with createSubset() function. ")

  assay(
    .internalAssay(.subsets(x)[[subsetName]]),
    withDimnames = FALSE,
    subsetAssayName) <- inputMatrix
  rownames(.internalAssay(.subsets(x)[[subsetName]])) <- rownames(inputMatrix)
  colnames(.internalAssay(.subsets(x)[[subsetName]])) <- colnames(inputMatrix)
  
  return(x)
}

'.assay<-' <- function(x, i, ..., subsetAssayName = NULL, value){
  if ((nrow(value) != nrow(x))
      || (ncol(value) != ncol(x))) {
    message("Storing assay inside subset '", i, "'.")
    setSubsetAssay(
      x = x,
      subsetName = i,
      inputMatrix = value,
      subsetAssayName = subsetAssayName
    )
  }
  else{
    callNextMethod(...)
  }
}

.subsetParent <- function(x, subsetName){
  parentList <- list()
  if (!subsetName %in% subsetAssayNames(x)) {
    stop(subsetName,
         " does not exist in the subsets slot of the object.")
  }
  test <- !is.null(.subsets(x)[[subsetName]]) &&
    is.null(.parentAssay(.subsets(x)[[subsetName]]))
  if (test) {
    return(NULL)
  }
  parent <- subsetName
  while (TRUE) {
    parentList <- c(parentList, parent)
    if (!is.null(.subsets(x)[[parent]])) {
      parent <- .parentAssay(.subsets(x)[[parent]])
    }
    else{
      for (i in seq(subsetCount(x))) {
        if (parent %in% assayNames(.internalAssay(.subsets(x)[[i]]))) {
          parent <- .subsetName(.subsets(x)[[i]])
        }
      }
      parentList <- c(parentList, parent)
      parent <- .parentAssay(.subsets(x)[[parent]])
    }
    if (parent %in% assayNames(x)) {
      parentList <- c(parentList, parent)
      break
    }
  }
  parentList[[1]] <- NULL
  return(parentList)
}

.subsetSummary <- function(x){
  cat("Main assay(s):\n",
      assayNames(x),
      "\n\n")
  cat("Subset(s):\n")
  if (!is.null(subsetNames(x))) {
    Name <- list()
    Dimensions <- list()
    Parent <- list()
    Assays <- list()
    Metadata <- list()
    ReducedDims <- list()
    AltExperiments <- list()
    
    for (i in seq(length(subsetNames(x)))) {
      parent <- subsetParent(x, subsetAssayNames(x)[i])
      Name[[i]] <- subsetNames(x)[i]
      Parent[[i]] <-
        paste(unlist(parent), collapse = ' -> ')
      if (is.null(assayNames(.internalAssay(.subsets(x)[[i]])))) {
        Assays[[i]] <- ""
      }
      else{
        Assays[[i]] <-
          paste(unlist(assayNames(.internalAssay(.subsets(x)[[i]]))), collapse = ", ")
      }
      Dimensions[[i]] <-
        paste(unlist(subsetDim(x, subsetNames(x)[i])), collapse = ', ')
      if(inherits(x, "SingleCellExperiment")){
        ReducedDims[[i]] <-
          paste(unlist(reducedDimNames(x, subsetName = subsetNames(x)[i])), collapse = ", ")
        AltExperiments[[i]] <-
          paste(unlist(altExpNames(x, subsetName = subsetNames(x)[i])), collapse = ", ")
      }
    }
    
    df <- data.frame(
      Name = as.character(Name),
      Dim = as.character(Dimensions),
      Parent = as.character(Parent)
    )
    
    if (length(which(as.character(Assays) == "")) != subsetCount(x)) {
      df <- cbind(df, Assays = as.character(Assays))
    }
    
    if(inherits(x, "SingleCellExperiment")){
      if (length(which(as.character(AltExperiments) == "")) != subsetCount(x)) {
        df <- cbind(df, AltExperiments = as.character(AltExperiments))
      }
      
      if (length(which(as.character(ReducedDims) == "")) != subsetCount(x)) {
        df <- cbind(df, ReducedDims = as.character(ReducedDims))
      }
    }
    
    print(df)
  }
  else{
    cat("NULL\n")
  }
}

#' @importMethodsFrom SingleCellExperiment altExp
.altExp <- function(x, e, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod())
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExp(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @importMethodsFrom SingleCellExperiment altExps
.altExps <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExps(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @importMethodsFrom SingleCellExperiment altExpNames
.altExpNames <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @rdname rownames
setMethod(
  f = "rownames",
  signature = "ANY",
  definition = function(x, subsetName)
  {
    if (missing(subsetName)) {
      BiocGenerics::rownames(x)
    }
    else{
      if (subsetName %in% subsetNames(x)) {
        BiocGenerics::rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
      }
      else if (subsetName %in% subsetAssayNames(x)) {
        subsetName <- .getParentAssayName(x, subsetName)
        BiocGenerics::rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
      }
    }
  }
)

#' @rdname colnames
setMethod(
  f = "colnames",
  signature = "ANY",
  definition = function(x, subsetName)
  {
    if (missing(subsetName)) {
      BiocGenerics::colnames(x)
    }
    else{
      if (subsetName %in% subsetNames(x)) {
        BiocGenerics::colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
      }
      else if (subsetName %in% subsetAssayNames(x)) {
        subsetName <- .getParentAssayName(x, subsetName)
        BiocGenerics::colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
      }
    }
  }
)

#' @importMethodsFrom SingleCellExperiment altExpNames<-
'.altExpNames<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @importMethodsFrom SingleCellExperiment altExp<-
'.altExp<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod())
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExp(.internalAssay(.subsets(x)[[subsetName]]), ...) <- value
  return(x)
}

#' @importMethodsFrom SingleCellExperiment altExps<-
'.altExps<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExps(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @importMethodsFrom S4Vectors metadata<-
'.metadata<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

.subsetRowData <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    #is a subset
    out <-
      rowData(x)[.rowIndices(.subsets(x)[[subsetName]]), , drop = FALSE]
    out <-
      cbind(out, rowData(.internalAssay(.subsets(x)[[subsetName]])))
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    #is a subset assay
    subsetName <- .getParentAssayName(x, subsetName)
    out <-
      rowData(x)[.rowIndices(.subsets(x)[[subsetName]]), , drop = FALSE]
    out <-
      cbind(out, rowData(.internalAssay(.subsets(x)[[subsetName]])))
  }
  else{
    #neither a subset nor a subset assay
    stop("Neither a subset nor a subsetAssay.")
  }
  return(out)
}

.getParentAssayName <- function(x, childAssayName) {
  for (i in seq(length(.subsets(x)))) {
    if (childAssayName %in% SummarizedExperiment::assayNames(.internalAssay(.subsets(x)[[i]]))) {
      return(.subsetName(.subsets(x)[[i]]))
    }
  }
}

.subsetColData <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    #is a subset
    out <-
      colData(x)[.colIndices(.subsets(x)[[subsetName]]), , drop = FALSE]
    out <-
      cbind(out, colData(.internalAssay(.subsets(x)[[subsetName]])))
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    #is a subset assay
    subsetName <- .getParentAssayName(x, subsetName)
    out <-
      colData(x)[.colIndices(.subsets(x)[[subsetName]]), , drop = FALSE]
    out <-
      cbind(out, colData(.internalAssay(.subsets(x)[[subsetName]])))
  }
  else{
    #neither a subset nor a subset assay
    stop("Neither a subset nor a subsetAssay.")
  }
  return(out)
}

#' @importMethodsFrom SingleCellExperiment reducedDim
.reducedDim <- function(x, type, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  reducedDim(.internalAssay(.subsets(x)[[subsetName]]), type)
}

#' @importMethodsFrom SingleCellExperiment reducedDim<-
'.reducedDim<-' <- function(x, type, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  reducedDim(.internalAssay(.subsets(x)[[subsetName]]), type, ...) <- value
  return(x)
}

#' @importMethodsFrom SingleCellExperiment reducedDims
.reducedDims <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  reducedDims(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @title reducedDims<-
#' @description A wrapper to the \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the \code{reducedDims} should be set to. If \code{missing}, \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method is called on the main object.
#' @param value A \code{list} of values to set to \code{reducedDims}.
#' @return Updated input object with \code{reducedDims} set.
#' @export
setReplaceMethod("reducedDims", "ANY", function(x, subsetName, value) {
  if (!missing(subsetName)) {
    SingleCellExperiment::reducedDims(.internalAssay(.subsets(x)[[subsetName]])) <-
      value
  }
  else{
    SingleCellExperiment::reducedDims(x) <- value
  }
  return(x)
})

.rowData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  subsetRowData(x, subsetName)
}

.colData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  subsetColData(x, subsetName)
}

#' @importMethodsFrom SummarizedExperiment rowData<-
'.rowData<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  rowData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @importMethodsFrom SummarizedExperiment colData<-
'.colData<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  colData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

.subsetColnames <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
  }
}

'.subsetColnames<-' <- function(x, subsetName, value){
  if (subsetName %in% subsetNames(x)) {
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])] <- value
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])] <- value
  }
  return(x)
}

.subsetRownames <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
  }
}

'.subsetRownames<-' <- function(x, subsetName, value){
  if (subsetName %in% subsetNames(x)) {
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])] <- value
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])] <- value
  }
  return(x)
}

#helpers

#checks if subset
.isSubset <- function(x, subsetName){
  if(!subsetName %in% subsetNames(x))
    stop(subsetName," subset does not exist.")
}

#subsets accessor (ExperimentSubset)
.subsets <- function(x) x@subsets

#subsets setter (ExperimentSubset)
'.subsets<-' <- function(x, value){
  x@subsets <- value
  return(x)
}

#subsetName accessor (AssaySubset)
.subsetName <- function(x) x@subsetName

#rowIndices accessor (AssaySubset)
.rowIndices <- function(x) x@rowIndices

#colIndices accessor (AssaySubset)
.colIndices <- function(x) x@colIndices

#parentAssay accessor (AssaySubset)
.parentAssay <- function(x) x@parentAssay

#internalAssay accessor (AssaySubset)
.internalAssay <- function(x) x@internalAssay

#internalAssay setter (AssaySubset)
'.internalAssay<-' <- function(x, value){
  x@internalAssay <- value
  return(x)
}

.subsetParamsValidity <- function(x,
                                  subsetName,
                                  rows,
                                  cols,
                                  parentAssay){
  if(subsetName %in% subsetNames(x))
    stop("A subset with the specified 'subsetName' parameter already exists
         in the object. Subset names must be unique.")
  
  if(!is.character(subsetName))
    stop("'subsetName' parameter must be a unique a character value.")
  
  testForRows <- is.null(rows) || is.numeric(rows) || is.character(rows)
  if(!testForRows)
    stop("'rows' parameter must be either a 'NULL' to include all rows,
         or a numeric vector or a character vector that specify the rows to 
         include in the subset.")
  
  testForCols <- is.null(cols) || is.numeric(cols) || is.character(cols)
  if(!testForCols)
    stop("'cols' parameter must be either a 'NULL' to include all rows,
         or a numeric vector or a character vector that specify the columns to 
         include in the subset.")
  
  testForParentAssay <- is.null(parentAssay) || is.character(parentAssay)
  if(!testForParentAssay)
    stop("'parentAssay' parameter can either be 'NULL' to use the default
         assay in the input object or a character value that specifies the
         parentAssay to use from parent object.")
  
}

.rows <- function(x, rows, assayName){
  tempRows <- NULL
  if (is.character(rows)) {
    tempRows <-
      match(rows, rownames(
        assay(x, withDimnames = TRUE, assayName)
      ))
  }
  else if (is.null(rows)) {
    tempRows <-
      seq(1, dim(
        assay(x, withDimnames = FALSE, assayName)
      )[1])
  }
  else{
    tempRows <- rows
  }
  
  dimR <- dim(assay(x, withDimnames = FALSE, assayName))[1]
  testRows <- any(tempRows > dimR)
  if(is.na(testRows))
    stop("Selected rows not available in the specified assay.")
  if(testRows)
    stop("Selected rows not available in the specified assay.")
  if(length(tempRows) > dimR)
    stop("Selected rows not available in the specified assay.")
  
  return(tempRows)
}

.cols <- function(x, cols, assayName){
  tempCols <- NULL
  if (is.character(cols)) {
    tempCols <-
      match(cols, colnames(
        assay(x, withDimnames = TRUE, assayName)
      ))
  }
  else if (is.null(cols)) {
    tempCols <-
      seq(1, dim(
        assay(x, withDimnames = FALSE, assayName)
      )[2])
  }
  else{
    tempCols <- cols
  }
  
  dimC <- dim(assay(x, withDimnames = FALSE, assayName))[2]
  testCols <- any(tempCols > dimC)
  if(is.na(testCols))
    stop("Selected columns not available in the specified assay.")
  if(testCols)
    stop("Selected columns not available in the specified assay.")
  if(length(tempCols) > dimC)
    stop("Selected columns not available in the specified assay.")
  
  return(tempCols)
}