
#' @title AssaySubset constructor
#' @description Constructor for creating a experiment object internally by the
#'   \code{ExperimentSubset} object. Should not be used directly by the user.
#' @param subsetName \code{character(1)} Name of the subset.
#' @param rowIndices \code{vector("numeric")} Indices of the rows to include in
#'   the subset.
#' @param colIndices \code{vector("numeric")} Indices of the columns to include
#'   in the subset.
#' @param parentAssay \code{character(1)} Name of the parent of this subset.
#' @param internalAssay An internal Experiment object to store additional
#'   subset data.
#' @return A \code{AssaySubset} object.
#' @keywords internal
AssaySubset <- function(subsetName = "subset",
                        rowIndices = NULL,
                        colIndices = NULL,
                        parentAssay = "counts",
                        internalAssay = NULL)
{
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
#' @description This constructor function sets up an \code{ExperimentSubset}
#'   object from a \code{SummarizedExperiment},
#'   \code{RangedSummarizedExperiment}, \code{SingleCellExperiment},
#'   \code{TreeSummarizedExperiment} or \code{SpatialExperiment} object, or
#'   from a named \code{list} of assay matrices. A subset can also be created
#'   directly by passing a named \code{list} to the \code{subset} parameter,
#'   with elements named \code{subsetName}, \code{rows}, \code{cols} and
#'   \code{parentAssay}.
#' @param x An experiment object of one of the five supported classes
#'   (subclasses of them are not supported), or a named \code{list} of assay
#'   matrices, which is passed to the \code{SingleCellExperiment} constructor
#'   as its assays.
#' @param subset A named \code{list} to create a subset from within the
#'   constructor, with elements \code{subsetName} and, optionally,
#'   \code{rows}, \code{cols} and \code{parentAssay}. Leave out
#'   \code{rows}, \code{cols} or \code{parentAssay}, or set them to
#'   \code{NULL}, to use all rows, all columns or the first assay. The
#'   default, a list of \code{NA} values, creates no subset; an \code{NA}
#'   in a list you supply is an error.
#' @return An object of the matching \code{Subset*} class, for example a
#'   \code{SubsetSingleCellExperiment} for a \code{SingleCellExperiment} or
#'   a list of assays.
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
    .checkAssayNames(x)
    es <- as(x, paste0("Subset", class(x)))
  }
  else{
    sce <- SingleCellExperiment::SingleCellExperiment(x)
    .checkAssayNames(sce)
    es <- .SubsetSingleCellExperiment(sce)
  }
  
  # the default subset is all NA, which means no subset
  if (!all(is.na(subset))) {
    if (anyNA(subset))
      stop("'subset' must not contain NA. Leave out rows, cols or ",
           "parentAssay, or set them to NULL, to use all rows, all columns ",
           "or the first assay.")
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

# Stops if any assay of x has no name
.checkAssayNames <- function(x) {
  if (is.null(assayNames(x))
      || "" %in% assayNames(x))
    stop("Experiment objects with unnamed assays are not supported. ",
         "Use assayNames<- setter method to set assay names before ",
         "creating ES object.")
}

#' @importMethodsFrom S4Vectors metadata
.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]]))
}

.assay <- function(x, i, ...){
  out <- NULL
  args <- list(...)
  withDimnames <- !isFALSE(args$withDimnames)
  # assay(x, "subset", subsetAssayName = "a") means "subset/a"
  if (!is.null(args$subsetAssayName)) {
    .isSubset(x, i)
    i <- paste0(i, "/", args$subsetAssayName)
  }
  ref <- .resolveAssay(x, i)
  if (is.null(ref))
    stop("'", i, "' is not the name of an assay, a subset or an assay ",
         "stored in a subset.")
  #look at main assays
  if (ref$kind == "main") {
    out <- callNextMethod(x, i, ...)
  } else if (ref$kind == "subset") {
    #look at subsets
    subsetName <- ref$subset
    i <- .parentAssay(.subsets(x)[[subsetName]]) 
    if (is.null(i)) {
      out <- assay(x = .internalAssay(.subsets(x)[[subsetName]]) , i = "counts", ...)
    }
    else{
      out <- assay(x = x, i = i, ...)
      out <- out[.rowIndices(.subsets(x)[[subsetName]]),
                 .colIndices(.subsets(x)[[subsetName]]), drop = FALSE]
      # use the subset's own names, which subsetRownames<- may have changed
      internal <- .internalAssay(.subsets(x)[[subsetName]])
      if (withDimnames && !is.null(rownames(internal)))
        rownames(out) <- rownames(internal)
      if (withDimnames && !is.null(colnames(internal)))
        colnames(out) <- colnames(internal)
    }
  }
  #look inside subsets
  else{
    # stored subset assays take the subset's row and column names
    out <- assay(.internalAssay(.subsets(x)[[ref$subset]]), ref$assay,
                 withDimnames = withDimnames)
  }
  out
}

.createSubset <- function(x,
                          subsetName,
                          rows,
                          cols,
                          parentAssay){
  # clean the name before the uniqueness check and use it as the key
  if (is.character(subsetName) && any(grepl("\\s", subsetName))) {
    subsetName <- gsub("\\s", "", subsetName)
    warning("Removing spaces from the specified subsetName.")
  }

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
    # stops if parentAssay is a name stored in more than one subset
    if (!is.null(.resolveAssay(x, parentAssay))) {
      tempAssay <- parentAssay
    }
    else{
      stop("Input parentAssay does not exist.")
    }
  }

  rows <- .rows(x, rows, tempAssay)
  cols <- .cols(x, cols, tempAssay)
  # rows and cols are relative to the parent; names come from the main object
  rootRows <- .toRootIndices(x, parentAssay, rows, .rowIndices)
  rootCols <- .toRootIndices(x, parentAssay, cols, .colIndices)
  
  a <- list(Matrix(
    nrow = length(rows),
    ncol = length(cols),
    data = 0,
    sparse = TRUE))
  names(a) <- "temp"
  
  internalAssay <- SingleCellExperiment(assays = a)
  if(inherits(x, "SpatialExperiment"))
    internalAssay$sample_id <- x$sample_id[rootCols]
  internalAssay <- as(internalAssay, gsub("Subset", "", class(x)))
  
  rownames(internalAssay) <- rownames(x)[rootRows]
  colnames(internalAssay) <- colnames(x)[rootCols]
  
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
  if (is.character(subsetAssayName)
      && grepl("/", subsetAssayName, fixed = TRUE))
    stop("'subsetAssayName' can't contain '/', which separates a subset ",
         "from an assay name.")
  if (is.character(subsetAssayName)
      && subsetAssayName %in% c(assayNames(x), subsetNames(x)))
    stop("'subsetAssayName' can't be the name of a main assay or a subset; ",
         "'", subsetAssayName, "' is one.")

  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (!identical(as.integer(dim(inputMatrix)), as.integer(dim(internal))))
    stop("'inputMatrix' must have the same dimensions as subset '",
         subsetName, "' (", nrow(internal), " x ", ncol(internal), ").")
  # the subset keeps its own names; inputMatrix may have none or the same
  .checkInputNames(rownames(inputMatrix), rownames(internal),
                   "rownames", subsetName)
  .checkInputNames(colnames(inputMatrix), colnames(internal),
                   "colnames", subsetName)

  # a subset of an object without names takes the input's names
  if (is.null(rownames(internal)))
    rownames(internal) <- rownames(inputMatrix)
  if (is.null(colnames(internal)))
    colnames(internal) <- colnames(inputMatrix)
  assay(internal, withDimnames = FALSE, subsetAssayName) <- inputMatrix
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  
  return(x)
}

# Stops unless inputNames are NULL or the same as the subset's current
# names, in the same order. what is "rownames" or "colnames".
.checkInputNames <- function(inputNames, currentNames, what, subsetName) {
  if (is.null(inputNames) || is.null(currentNames))
    return(invisible(NULL))
  if (!identical(as.character(inputNames), as.character(currentNames)))
    stop("The ", what, " of 'inputMatrix' must be NULL or the same as the ",
         what, " of subset '", subsetName, "', in the same order.")
}

'.assay<-' <- function(x, i, ..., subsetAssayName = NULL, value){
  # a subset name stores the assay inside that subset
  if (i %in% subsetNames(x)) {
    if (is.null(subsetAssayName))
      stop("'subsetAssayName' is needed to store an assay inside subset '",
           i, "'.")
    message("Storing assay inside subset '", i, "'.")
    setSubsetAssay(
      x = x,
      subsetName = i,
      inputMatrix = value,
      subsetAssayName = subsetAssayName
    )
  }
  else{
    if (!is.null(subsetAssayName))
      stop("'subsetAssayName' is only used when 'i' is a subset name; '",
           i, "' is not a subset.")
    callNextMethod(x, i, ..., value = value)
  }
}

.subsetParent <- function(x, subsetName){
  ref <- .resolveAssay(x, subsetName)
  if (is.null(ref) || ref$kind == "main") {
    stop(subsetName,
         " does not exist in the subsets slot of the object.")
  }
  test <- !is.null(.subsets(x)[[subsetName]]) &&
    is.null(.parentAssay(.subsets(x)[[subsetName]]))
  if (test) {
    return(NULL)
  }
  parentList <- list()
  current <- subsetName
  # an assay stored in a subset has that subset as its first parent
  if (is.null(.subsets(x)[[current]])) {
    current <- .parentSubsetName(x, current)
    parentList <- c(parentList, current)
  }
  # names resolve as in .assay(); each step moves up one subset, so a
  # longer walk means the chain has a cycle
  for (step in seq_len(subsetCount(x) + 1)) {
    parent <- .parentAssay(.subsets(x)[[current]])
    if (is.null(parent))
      return(parentList)
    if (parent %in% assayNames(x))
      return(c(parentList, parent))
    holder <- .parentSubsetName(x, parent)
    if (is.null(holder))
      stop("The parent '", parent, "' of '", current, "' cannot be found.")
    if (parent != holder)
      parentList <- c(parentList, parent)
    parentList <- c(parentList, holder)
    current <- holder
  }
  stop("The parent chain of '", subsetName, "' does not lead to a main ",
       "assay; an assay name stored in more than one subset can cause ",
       "this.")
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
      parent <- subsetParent(x, subsetNames(x)[i])
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
  if (!"subsetName" %in% names(arglist)) {
    if (missing(e))
      return(callNextMethod(x, ...))
    return(callNextMethod(x, e, ...))
  }
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  # a missing e can't be passed into a second S4 dispatch
  if (missing(e))
    return(do.call(altExp, c(list(internal), .otherArgs(...))))
  do.call(altExp, c(list(internal, e), .otherArgs(...)))
}

#' @importMethodsFrom SingleCellExperiment altExps
.altExps <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  do.call(altExps, c(list(internal), .otherArgs(...)))
}

#' @importMethodsFrom SingleCellExperiment altExpNames
.altExpNames <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @importMethodsFrom SingleCellExperiment altExpNames<-
'.altExpNames<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @importMethodsFrom SingleCellExperiment altExp<-
'.altExp<-' <- function(x, e, ..., value){
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist)) {
    if (missing(e))
      return(callNextMethod(x, ..., value = value))
    return(callNextMethod(x, e, ..., value = value))
  }
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (missing(e))
    args <- c(list(internal), .otherArgs(...), list(value = value))
  else
    args <- c(list(internal, e), .otherArgs(...), list(value = value))
  .internalAssay(.subsets(x)[[subsetName]]) <- do.call(`altExp<-`, args)
  return(x)
}

#' @importMethodsFrom SingleCellExperiment altExps<-
'.altExps<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  args <- c(list(internal), .otherArgs(...), list(value = value))
  .internalAssay(.subsets(x)[[subsetName]]) <- do.call(`altExps<-`, args)
  return(x)
}

#' @importMethodsFrom S4Vectors metadata<-
'.metadata<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

# Where an assay reference points: a main assay (kind "main"), a subset
# (kind "subset"), or an assay stored in a subset (kind "stored"). A stored
# assay can be named as "subset/assay". Returns NULL if nothing matches, and
# stops if an unqualified name is stored in more than one subset.
.resolveAssay <- function(x, name) {
  if (grepl("/", name, fixed = TRUE)) {
    parts <- strsplit(name, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2 || !parts[1] %in% subsetNames(x))
      return(NULL)
    stored <- assayNames(.internalAssay(.subsets(x)[[parts[1]]]))
    if (!parts[2] %in% stored)
      return(NULL)
    return(list(kind = "stored", subset = parts[1], assay = parts[2]))
  }
  if (name %in% assayNames(x))
    return(list(kind = "main", subset = NULL, assay = name))
  if (name %in% subsetNames(x))
    return(list(kind = "subset", subset = name, assay = NULL))
  holders <- .assayHolders(x, name)
  if (length(holders) == 0)
    return(NULL)
  if (length(holders) > 1)
    stop("'", name, "' is stored in subsets ",
         paste0("'", holders, "'", collapse = ", "),
         "; use a qualified name such as \"", holders[1], "/", name, "\".")
  list(kind = "stored", subset = holders, assay = name)
}

# Names of the subsets that store an assay called name
.assayHolders <- function(x, name) {
  holders <- character(0)
  for (subset in subsetNames(x)) {
    if (name %in% assayNames(.internalAssay(.subsets(x)[[subset]])))
      holders <- c(holders, subset)
  }
  holders
}

# Name of the subset that parentAssay refers to: the subset itself, or the
# subset storing that assay. NULL if parentAssay is a main assay or cannot
# be found.
.parentSubsetName <- function(x, parentAssay) {
  if (is.null(parentAssay))
    return(NULL)
  ref <- .resolveAssay(x, parentAssay)
  if (is.null(ref) || ref$kind == "main")
    return(NULL)
  ref$subset
}

# Converts indices relative to parentAssay into positions in the main
# object by following the chain of parent subsets. getIndices is
# .rowIndices or .colIndices.
.toRootIndices <- function(x, parentAssay, indices, getIndices) {
  for (step in seq_len(subsetCount(x) + 1)) {
    parentName <- .parentSubsetName(x, parentAssay)
    if (is.null(parentName))
      return(indices)
    parent <- .subsets(x)[[parentName]]
    indices <- getIndices(parent)[indices]
    parentAssay <- .parentAssay(parent)
  }
  stop("The parent chain of the subset does not lead to a main assay.")
}

# Row or column indices of a subset as positions in the main object
.rootIndices <- function(x, subsetName, getIndices) {
  subset <- .subsets(x)[[subsetName]]
  .toRootIndices(x, .parentAssay(subset), getIndices(subset), getIndices)
}

.subsetRowData <- function(x, subsetName, parentRowData){
  .isSubset(x, subsetName)
  if(missing(parentRowData)
     || is.null(parentRowData))
    parentRowData = FALSE
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  out <- rowData(internal)
  if (parentRowData) {
    rootRows <- .rootIndices(x, subsetName, .rowIndices)
    parentRows <- rowData(x)[rootRows, , drop = FALSE]
    out <- .tagParentColumns(cbind(parentRows, out), ncol(parentRows))
    # the subset's own row names, which subsetRownames<- may have changed
    rownames(out) <- rownames(internal)
  }
  return(out)
}

.subsetColData <- function(x, subsetName, parentColData){
  .isSubset(x, subsetName)
  if(missing(parentColData)
     || is.null(parentColData))
    parentColData = FALSE
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  out <- colData(internal)
  if (parentColData) {
    rootCols <- .rootIndices(x, subsetName, .colIndices)
    parentCols <- colData(x)[rootCols, , drop = FALSE]
    out <- .tagParentColumns(cbind(parentCols, out), ncol(parentCols))
    # the subset's own column names, which subsetColnames<- may have changed
    rownames(out) <- colnames(internal)
  }
  return(out)
}

'.subsetColData<-' <- function(x, subsetName, value){
  .isSubset(x, subsetName)
  value <- .dropParentColumns(value, "colData")
  # colData<- takes the subset's column names from value's rownames
  value <- .keepSubsetColnames(value, .subsetColnames(x, subsetName),
                               subsetName)
  colData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

'.subsetRowData<-' <- function(x, subsetName, value){
  .isSubset(x, subsetName)
  value <- .dropParentColumns(value, "rowData")
  rowData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

# Marks the first nParent columns of a subset's rowData or colData as
# coming from the main object, in the DataFrame's column metadata, so the
# setters can leave them out when the result is written back
.tagParentColumns <- function(out, nParent) {
  fromParent <- seq_len(ncol(out)) <= nParent
  columnInfo <- S4Vectors::mcols(out)
  if (is.null(columnInfo))
    columnInfo <- S4Vectors::DataFrame(fromParentData = fromParent)
  else
    columnInfo$fromParentData <- fromParent
  S4Vectors::mcols(out) <- columnInfo
  out
}

# Leaves out the columns .tagParentColumns() marked, so writing back a
# result read with the parent's data stores no copies of it
.dropParentColumns <- function(value, what) {
  if (!methods::is(value, "DataFrame"))
    return(value)
  columnInfo <- S4Vectors::mcols(value)
  if (is.null(columnInfo) || !"fromParentData" %in% colnames(columnInfo))
    return(value)
  fromParent <- columnInfo$fromParentData %in% TRUE
  if (any(fromParent))
    message("Not storing columns that come from the main object's ", what,
            ": ", paste(colnames(value)[fromParent], collapse = ", "),
            ". Change them on the main object with ", what, "(x).")
  .untagParentColumns(value[, !fromParent, drop = FALSE])
}

# Removes the .tagParentColumns() mark, which is only for the setters
.untagParentColumns <- function(value) {
  if (!methods::is(value, "DataFrame"))
    return(value)
  columnInfo <- S4Vectors::mcols(value)
  if (is.null(columnInfo) || !"fromParentData" %in% colnames(columnInfo))
    return(value)
  columnInfo$fromParentData <- NULL
  if (ncol(columnInfo) == 0)
    columnInfo <- NULL
  S4Vectors::mcols(value) <- columnInfo
  value
}

# Setting a subset's colData must not rename its columns: value's rownames
# must be NULL, and are then filled in, or the same as the subset's names.
# A value of the wrong size is left for colData<- to reject.
.keepSubsetColnames <- function(value, currentNames, subsetName) {
  if (is.null(value))
    return(value)
  if (is.null(rownames(value))) {
    if (nrow(value) == length(currentNames))
      rownames(value) <- currentNames
  } else if (!is.null(currentNames)
             && !identical(rownames(value), currentNames)) {
    stop("The rownames of 'value' must be NULL or the column names of ",
         "subset '", subsetName, "', in the same order.")
  }
  value
}

#' @importMethodsFrom SingleCellExperiment reducedDim
.reducedDim <- function(x, type, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, type, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  # a missing type can't be passed into a second S4 dispatch
  if (missing(type))
    return(do.call(reducedDim, c(list(internal), .otherArgs(...))))
  do.call(reducedDim, c(list(internal, type), .otherArgs(...)))
}

#' @importMethodsFrom SingleCellExperiment reducedDim<-
'.reducedDim<-' <- function(x, type, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, type, ..., value = value))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (missing(type))
    args <- c(list(internal), .otherArgs(...), list(value = value))
  else
    args <- c(list(internal, type), .otherArgs(...), list(value = value))
  .internalAssay(.subsets(x)[[subsetName]]) <- do.call(`reducedDim<-`, args)
  return(x)
}

#' @importMethodsFrom SingleCellExperiment reducedDims
.reducedDims <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  do.call(reducedDims, c(list(internal), .otherArgs(...)))
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
    .isSubset(x, subsetName)
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
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  parentRowData = arglist[["parentRowData"]]
  .isSubset(x, subsetName)
  out <- .subsetRowData(x, subsetName, parentRowData)
  if (isFALSE(arglist$use.names))
    rownames(out) <- NULL
  out
}

.colData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  parentColData = arglist[["parentColData"]]
  .isSubset(x, subsetName)
  .subsetColData(x, subsetName, parentColData)
}

#' @importMethodsFrom SummarizedExperiment rowData<-
'.rowData<-' <- function(x, ..., value){
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist)) {
    value <- .untagParentColumns(value)
    return(callNextMethod(x, ..., value = value))
  }
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  if (is.data.frame(value))
    value <- S4Vectors::DataFrame(value, check.names = FALSE)
  .subsetRowData(x, subsetName) <- value
  return(x)
}

#' @importMethodsFrom SummarizedExperiment colData<-
'.colData<-' <- function(x, ..., value){
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist)) {
    if (is.data.frame(value))
      stop("'value' must be a DataFrame or NULL.")
    value <- .untagParentColumns(value)
    return(callNextMethod(x, ..., value = value))
  }
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  if (is.data.frame(value))
    value <- S4Vectors::DataFrame(value, check.names = FALSE)
  .subsetColData(x, subsetName) <- value
  return(x)
}

.subsetColnames <- function(x, subsetName){
  subsetName <- .owningSubset(x, subsetName)
  colnames(.internalAssay(.subsets(x)[[subsetName]]))
}

'.subsetColnames<-' <- function(x, subsetName, value){
  subsetName <- .owningSubset(x, subsetName)
  colnames(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

.subsetRownames <- function(x, subsetName){
  subsetName <- .owningSubset(x, subsetName)
  rownames(.internalAssay(.subsets(x)[[subsetName]]))
}

'.subsetRownames<-' <- function(x, subsetName, value){
  subsetName <- .owningSubset(x, subsetName)
  rownames(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

.subsetRowLinkData <- function(x, subsetName, parentLinkRowData){
  .isSubset(x, subsetName)
  if(missing(parentLinkRowData)
     || is.null(parentLinkRowData))
    parentLinkRowData = FALSE
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  out <- rowLinks(internal)
  # a main object without a row tree has no links to add
  if (parentLinkRowData && !is.null(rowLinks(x))) {
    rootIndices <- .rootIndices(x, subsetName, .rowIndices)
    parentLinks <- rowLinks(x)[rootIndices, , drop = FALSE]
    if (is.null(out))
      out <- parentLinks
    else
      out <- cbind(parentLinks, out)
  }
  if (!is.null(out))
    rownames(out) <- rownames(internal)
  return(out)
}

.subsetColLinkData <- function(x, subsetName, parentLinkColData){
  .isSubset(x, subsetName)
  if(missing(parentLinkColData)
     || is.null(parentLinkColData))
    parentLinkColData = FALSE
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  out <- colLinks(internal)
  # a main object without a column tree has no links to add
  if (parentLinkColData && !is.null(colLinks(x))) {
    rootIndices <- .rootIndices(x, subsetName, .colIndices)
    parentLinks <- colLinks(x)[rootIndices, , drop = FALSE]
    if (is.null(out))
      out <- parentLinks
    else
      out <- cbind(parentLinks, out)
  }
  if (!is.null(out))
    rownames(out) <- colnames(internal)
  return(out)
}

.subsetSpatialCoords <- function(x, subsetName){
  .isSubset(x, subsetName)
  rootCols <- .rootIndices(x, subsetName, .colIndices)
  out <- spatialCoords(x)[rootCols, , drop = FALSE]
  return(out)
}

#' @importMethodsFrom SpatialExperiment spatialCoords
.spatialCoords <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  .subsetSpatialCoords(x, subsetName)
}

.subsetSpatialData <- function(x, subsetName, parentSpatialData){
  .isSubset(x, subsetName)
  if(missing(parentSpatialData)
     || is.null(parentSpatialData))
    parentSpatialData = FALSE
  # a DataFrame, since spot names can repeat across samples
  nSpots <- subsetDim(x, subsetName)[2]
  out <- S4Vectors::DataFrame(row.names = as.character(seq_len(nSpots)))
  if(parentSpatialData) {
    rootCols <- .rootIndices(x, subsetName, .colIndices)
    out <- spatialData(x)[rootCols, , drop = FALSE]
  }
  subsetSpatialData <- spatialData(.internalAssay(.subsets(x)[[subsetName]]))
  if(ncol(subsetSpatialData) > 0)
    out <- cbind(out, subsetSpatialData)
  # the subset's own column names, which subsetColnames<- may have changed
  spots <- colnames(.internalAssay(.subsets(x)[[subsetName]]))
  if (!is.null(spots))
    rownames(out) <- spots
  return(out)
}

#' @importMethodsFrom SpatialExperiment spatialData
.spatialData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  parentSpatialData = arglist[["parentSpatialData"]]
  .isSubset(x, subsetName)
  .subsetSpatialData(x, subsetName, parentSpatialData)
}

#' @importMethodsFrom SpatialExperiment spatialData<-
'.subsetSpatialData<-' <- function(x, subsetName, value){
  .isSubset(x, subsetName)
  spatialData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#helpers

# Arguments in ... other than subsetName, to pass to a subset's internal
# object, whose methods may not accept subsetName
.otherArgs <- function(...) {
  args <- list(...)
  if (is.null(names(args)))
    return(args)
  args[names(args) != "subsetName"]
}

#checks if is a subset
.isSubset <- function(x, subsetName){
  if(!subsetName %in% subsetNames(x))
    stop(subsetName, " subset does not exist. 'subsetName' must be the ",
         "name of a subset.")
}

# Name of the subset that holds subsetName: the subset itself, or the
# subset storing an assay of that name. Stops if there is none.
.owningSubset <- function(x, subsetName) {
  if (subsetName %in% subsetNames(x))
    return(subsetName)
  ref <- .resolveAssay(x, subsetName)
  if (!is.null(ref) && ref$kind == "stored")
    return(ref$subset)
  .isSubset(x, subsetName)
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

#checks for validity of subset parameters before creating one
.subsetParamsValidity <- function(x,
                                  subsetName,
                                  rows,
                                  cols,
                                  parentAssay){
  if (!is.character(subsetName) || length(subsetName) != 1
      || is.na(subsetName) || !nzchar(subsetName))
    stop("'subsetName' must be a single, non-empty character string.")

  if(subsetName %in% subsetNames(x))
    stop("A subset with the specified 'subsetName' parameter already exists
         in the object. Subset names must be unique.")

  if (subsetName %in% c(assayNames(x), subsetAssayNames(x)))
    stop("'subsetName' must differ from the names of existing assays; '",
         subsetName, "' is already an assay name.")

  if (grepl("/", subsetName, fixed = TRUE))
    stop("'subsetName' can't contain '/', which separates a subset from an ",
         "assay name.")
  
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
  
  testForParentAssay <- is.null(parentAssay) ||
    (is.character(parentAssay) && length(parentAssay) == 1 &&
       !is.na(parentAssay))
  if(!testForParentAssay)
    stop("'parentAssay' parameter can either be 'NULL' to use the default
         assay in the input object or a character value that specifies the
         parentAssay to use from parent object.")
  
}

.rows <- function(x, rows, assayName){
  parent <- assay(x, withDimnames = TRUE, assayName)
  dimR <- nrow(parent)
  if (is.character(rows)) {
    if (is.null(rownames(parent)))
      stop("'rows' are given as names, but the parent assay has no row ",
           "names.")
    tempRows <- match(rows, rownames(parent))
  }
  else if (is.null(rows)) {
    tempRows <- seq_len(dimR)
  }
  else{
    tempRows <- rows
  }
  
  testRows <- any(tempRows > dimR)
  if(is.na(testRows))
    stop("Selected rows not available in the specified assay.")
  if(testRows)
    stop("Selected rows not available in the specified assay.")
  if(length(tempRows) > dimR)
    stop("Selected rows not available in the specified assay.")
  if (length(tempRows) == 0)
    stop("'rows' must select at least one row.")
  if (any(tempRows < 1) || any(tempRows != round(tempRows)))
    stop("'rows' must be positive whole numbers or row names.")
  if (anyDuplicated(tempRows))
    stop("'rows' must not contain duplicates.")
  
  return(tempRows)
}

.cols <- function(x, cols, assayName){
  parent <- assay(x, withDimnames = TRUE, assayName)
  dimC <- ncol(parent)
  if (is.character(cols)) {
    if (is.null(colnames(parent)))
      stop("'cols' are given as names, but the parent assay has no column ",
           "names.")
    tempCols <- match(cols, colnames(parent))
  }
  else if (is.null(cols)) {
    tempCols <- seq_len(dimC)
  }
  else{
    tempCols <- cols
  }
  
  testCols <- any(tempCols > dimC)
  if(is.na(testCols))
    stop("Selected columns not available in the specified assay.")
  if(testCols)
    stop("Selected columns not available in the specified assay.")
  if(length(tempCols) > dimC)
    stop("Selected columns not available in the specified assay.")
  if (length(tempCols) == 0)
    stop("'cols' must select at least one column.")
  if (any(tempCols < 1) || any(tempCols != round(tempCols)))
    stop("'cols' must be positive whole numbers or column names.")
  if (anyDuplicated(tempCols))
    stop("'cols' must not contain duplicates.")
  
  return(tempCols)
}
