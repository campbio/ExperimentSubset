
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
  args <- list(...)
  withDimnames <- !isFALSE(args$withDimnames)
  # assay(x, "subset", subsetAssayName = "a") means "subset/a"
  if (!is.null(args$subsetAssayName)) {
    if (!is.character(args$subsetAssayName)
        || length(args$subsetAssayName) != 1)
      stop("'subsetAssayName' must be a single character string.")
    .isSubset(x, i)
    i <- paste0(i, "/", args$subsetAssayName)
  }
  ref <- .resolveAssay(x, i)
  if (is.null(ref))
    stop("'", i, "' is not the name of an assay, a subset or an assay ",
         "seen in a subset.")
  #look at main assays
  if (ref$kind == "main")
    return(callNextMethod(x, i, ...))
  #an assay as seen in a subset, its own or linked
  if (ref$kind == "visible")
    return(.subsetAssayView(x, ref$subset, ref$assay, withDimnames))
  #a subset: its view of the assay its parentAssay refers to, which is
  #its own if it, or a subset on the way, stored one with unlink = TRUE
  subsetName <- ref$subset
  if (is.null(.parentAssay(.subsets(x)[[subsetName]]))) {
    internal <- .internalAssay(.subsets(x)[[subsetName]])
    return(assay(x = internal, i = "counts", ...))
  }
  .subsetAssayView(x, subsetName, .primaryAssayName(x, subsetName),
                   withDimnames)
}

.createSubset <- function(x,
                          subsetName,
                          rows,
                          cols,
                          parentAssay){
  x <- .qualifyParents(x)
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
    # stops if parentAssay is a bare name of an assay stored in a subset
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

.subsetAssayNames <- function(x, subsetName = NULL) {
  # the assays seen in one subset, own and linked
  if (!is.null(subsetName)) {
    .isSubset(x, subsetName)
    return(names(.assayOwners(x, subsetName)))
  }
  # the subsets, then every subset's own assays as "subset/assay"
  tempNames <- names(.subsets(x))
  for (subset in names(.subsets(x))) {
    own <- assayNames(.internalAssay(.subsets(x)[[subset]]))
    if (length(own) > 0)
      tempNames <- c(tempNames, paste0(subset, "/", own))
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
                         subsetAssayName = NULL,
                         unlink = FALSE) {

  if (!subsetName %in% subsetNames(x))
    stop(subsetName, " does not exist in the subsets slot of the object. 
         You need to create a new subset with createSubset() function. ")
  x <- .qualifyParents(x)
  if (!is.character(subsetAssayName) || length(subsetAssayName) != 1
      || is.na(subsetAssayName) || !nzchar(subsetAssayName))
    stop("'subsetAssayName' must be a single, non-empty character string.")
  if (grepl("/", subsetAssayName, fixed = TRUE))
    stop("'subsetAssayName' can't contain '/', which separates a subset ",
         "from an assay name.")
  # a name linked from the parent is protected unless unlink = TRUE
  owners <- .assayOwners(x, subsetName)
  if (subsetAssayName %in% names(owners)
      && owners[[subsetAssayName]] != subsetName && !isTRUE(unlink)) {
    owner <- owners[[subsetAssayName]]
    source <- paste0("'", owner, "/", subsetAssayName, "'")
    if (owner == "")
      source <- paste0("the main assay '", subsetAssayName, "'")
    stop("'", subsetAssayName, "' in subset '", subsetName, "' is linked ",
         "to ", source, ". Set unlink = TRUE to store a separate '",
         subsetAssayName, "' computed on the cells of '", subsetName, "'.")
  }

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

# A subset as a standalone object of the main object's class: the main
# object cut to the subset's rows and columns, which keeps trees, ranges,
# spatial coordinates and sample_id, with the subset's own names, assays
# and data
#' @importMethodsFrom SummarizedExperiment assays assays<-
#' @importMethodsFrom TreeSummarizedExperiment rownames<- colnames<-
.getSubsetExperiment <- function(x, subsetName, parentData = TRUE) {
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  out <- as(x, class(internal))
  # the main object's assays and reduced data are replaced, so not cut;
  # rownames<- and colnames<- are the generics, which keep tree links
  assays(out, withDimnames = FALSE) <- S4Vectors::SimpleList()
  if (methods::is(out, "SingleCellExperiment")) {
    SingleCellExperiment::reducedDims(out) <- list()
    altExps(out) <- list()
  }
  out <- out[.rootIndices(x, subsetName, .rowIndices),
             .rootIndices(x, subsetName, .colIndices)]
  rownames(out) <- rownames(internal)
  colnames(out) <- colnames(internal)

  # every assay the subset sees, its primary assay first
  assayOrder <- names(.assayOwners(x, subsetName))
  primary <- intersect(.primaryAssayName(x, subsetName), assayOrder)
  assayOrder <- c(primary, setdiff(assayOrder, primary))
  views <- lapply(assayOrder, function(name) {
    .subsetAssayView(x, subsetName, name, withDimnames = FALSE)
  })
  names(views) <- assayOrder
  assays(out, withDimnames = FALSE) <- S4Vectors::SimpleList(views)

  rowData(out) <- .subsetRowData(x, subsetName, isTRUE(parentData))
  colData(out) <- .subsetColData(x, subsetName, isTRUE(parentData))
  if (methods::is(out, "SingleCellExperiment")) {
    SingleCellExperiment::reducedDims(out) <- reducedDims(internal)
    altExps(out) <- altExps(internal)
  }
  metadata(out) <- metadata(internal)
  out
}

# Stores what a tool added to an object from .getSubsetExperiment(); value
# replaces the subset's own data, and unchanged linked assays and the main
# object's tagged columns are left out
.setSubsetExperiment <- function(x, subsetName, value, unlink = FALSE) {
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (!identical(as.character(class(value)), as.character(class(internal))))
    stop("'value' must be a ", class(internal), ", the class of the main ",
         "object, as returned by getSubsetExperiment().")
  if (!identical(as.integer(dim(value)), as.integer(dim(internal))))
    stop("'value' must have the same dimensions as subset '", subsetName,
         "' (", nrow(internal), " x ", ncol(internal), "). To keep fewer ",
         "rows or columns, create a subset with createSubset().")
  # matching names show value's rows and columns weren't reordered
  for (what in c("rownames", "colnames")) {
    getNames <- match.fun(what)
    if (!is.null(getNames(internal))
        && !identical(getNames(value), getNames(internal)))
      stop("The ", what, " of 'value' must be the ", what, " of subset '",
           subsetName, "', in the same order.")
  }

  # own assays that value no longer has are removed, unless a subset was
  # made from them
  x <- .qualifyParents(x)
  owners <- .assayOwners(x, subsetName)
  kept <- assayNames(internal) %in% assayNames(value)
  .checkNotParents(x, subsetName, assayNames(internal)[!kept])
  assays(internal, withDimnames = FALSE) <-
    assays(internal, withDimnames = FALSE)[kept]
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  for (name in assayNames(value)) {
    input <- assay(value, name)
    linked <- name %in% names(owners) && owners[[name]] != subsetName
    if (linked && .sameAssay(input, .subsetAssayView(x, subsetName, name)))
      next
    x <- .setSubsetAssay(x, subsetName, input, name, unlink = unlink)
  }

  .subsetRowData(x, subsetName) <- .withoutParentColumns(rowData(value))
  .subsetColData(x, subsetName) <- .withoutParentColumns(colData(value))
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (methods::is(internal, "SingleCellExperiment")) {
    SingleCellExperiment::reducedDims(internal) <- reducedDims(value)
    altExps(internal) <- altExps(value)
  }
  metadata(internal) <- metadata(value)
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  x
}

# Stops if removing subsetName's own assays in assayNames would change or
# break another subset: one whose primary assay is one of them, reached
# directly ("s/scaled"), through the subset ("s") or through a descendant
# ("child/scaled")
.checkNotParents <- function(x, subsetName, assayNames) {
  used <- character(0)
  for (name in setdiff(subsetNames(x), subsetName)) {
    primary <- .primaryAssayName(x, name)
    owners <- .assayOwners(x, name)
    if (primary %in% assayNames && primary %in% names(owners)
        && owners[[primary]] == subsetName)
      used[name] <- primary
  }
  if (length(used) > 0)
    stop("'value' has no assay ",
         paste0("'", unique(used), "'", collapse = ", "),
         ", which subsets ", paste0("'", names(used), "'", collapse = ", "),
         " were made from. Keep it in 'value'.")
}

# TRUE if two assays hold the same values, whatever their dimnames, type
# (integer or double) or class (for example dense or sparse)
.sameAssay <- function(a, b) {
  dimnames(a) <- NULL
  dimnames(b) <- NULL
  if (identical(a, b))
    return(TRUE)
  same <- tryCatch(all.equal(a, b, tolerance = 0, check.attributes = FALSE),
                   error = function(e) FALSE)
  isTRUE(same)
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
  unlink <- isTRUE(list(...)$unlink)
  # "subset/assay" stores into that subset
  if (grepl("/", i, fixed = TRUE) && !i %in% assayNames(x)) {
    parts <- .splitQualified(i)
    if (is.null(parts) || !parts[1] %in% subsetNames(x))
      stop("'", i, "' is not \"subset/assay\" for an existing subset.")
    if (!is.null(subsetAssayName))
      stop("'subsetAssayName' can't be used with \"subset/assay\"; '", i,
           "' already names the assay.")
    return(setSubsetAssay(x, parts[1], value, parts[2], unlink = unlink))
  }
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
      subsetAssayName = subsetAssayName,
      unlink = unlink
    )
  }
  else{
    if (!is.null(subsetAssayName))
      stop("'subsetAssayName' is only used when 'i' is a subset name; '",
           i, "' is not a subset.")
    if (unlink)
      stop("'unlink' only applies when storing an assay in a subset.")
    callNextMethod(x, i, ..., value = value)
  }
}

.subsetParent <- function(x, subsetName){
  ref <- .resolveAssay(x, subsetName)
  if (is.null(ref) || ref$kind == "main") {
    stop(subsetName,
         " does not exist in the subsets slot of the object.")
  }
  current <- ref$subset
  parentList <- list()
  # an assay seen in a subset has that subset as its first parent
  if (ref$kind == "visible")
    parentList <- c(parentList, current)
  else if (is.null(.parentAssay(.subsets(x)[[current]])))
    return(NULL)
  # each step moves up one subset, so a longer walk means a cycle
  for (step in seq_len(subsetCount(x) + 1)) {
    parent <- .parentAssay(.subsets(x)[[current]])
    if (is.null(parent))
      return(parentList)
    container <- .parentSubsetName(x, parent)
    if (is.null(container))
      return(c(parentList, parent))
    if (parent != container)
      parentList <- c(parentList, parent)
    parentList <- c(parentList, container)
    current <- container
  }
  stop("The parent chain of '", subsetName, "' does not lead to a main ",
       "assay.")
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

# Splits "subset/assay" into its two names; NULL if name isn't of that form
.splitQualified <- function(name) {
  if (!grepl("/", name, fixed = TRUE))
    return(NULL)
  parts <- strsplit(name, "/", fixed = TRUE)[[1]]
  if (length(parts) != 2 || !all(nzchar(parts)) || endsWith(name, "/"))
    return(NULL)
  parts
}

# Where a reference points: a main assay (kind "main"), a subset (kind
# "subset"), or an assay as seen in a subset, named "subset/assay" (kind
# "visible", either the subset's own or linked from its parent). Returns
# NULL if nothing matches. A bare name that is only stored in a subset is
# an error that shows the qualified form.
.resolveAssay <- function(x, name) {
  if (name %in% assayNames(x))
    return(list(kind = "main", subset = NULL, assay = name))
  if (name %in% subsetNames(x))
    return(list(kind = "subset", subset = name, assay = NULL))
  if (grepl("/", name, fixed = TRUE)) {
    parts <- .splitQualified(name)
    found <- !is.null(parts) && parts[1] %in% subsetNames(x)
    if (found && parts[2] %in% names(.assayOwners(x, parts[1])))
      return(list(kind = "visible", subset = parts[1], assay = parts[2]))
    return(NULL)
  }
  holders <- .assayHolders(x, name)
  if (length(holders) > 0)
    stop("'", name, "' is stored in subset ",
         paste0("'", holders, "'", collapse = ", "), ", not in the main ",
         "object; name it as \"", holders[1], "/", name, "\" (assay() and ",
         "getSubsetAssay() also take subsetAssayName =).")
  NULL
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

# The container a parentAssay points into: NULL for the main object, or a
# subset name for a subset ("P") or an assay seen in a subset ("P/assay").
# Objects saved before parents were qualified may hold a bare stored-assay
# name; it resolves when one subset stores it. NULL if nothing matches.
.parentSubsetName <- function(x, parentAssay) {
  if (is.null(parentAssay) || parentAssay %in% assayNames(x))
    return(NULL)
  if (parentAssay %in% subsetNames(x))
    return(parentAssay)
  if (grepl("/", parentAssay, fixed = TRUE)) {
    parts <- .splitQualified(parentAssay)
    if (!is.null(parts) && parts[1] %in% subsetNames(x))
      return(parts[1])
    return(NULL)
  }
  holders <- .assayHolders(x, parentAssay)
  if (length(holders) > 1)
    stop("The parent '", parentAssay, "' is stored in subsets ",
         paste0("'", holders, "'", collapse = ", "), ".")
  if (length(holders) == 1)
    return(holders)
  NULL
}

# Subset names from the outermost container down to subsetName
.subsetChain <- function(x, subsetName) {
  chain <- subsetName
  for (step in seq_len(subsetCount(x) + 1)) {
    parent <- .parentSubsetName(x, .parentAssay(.subsets(x)[[chain[1]]]))
    if (is.null(parent))
      return(chain)
    chain <- c(parent, chain)
  }
  stop("The parent chain of '", subsetName, "' does not lead to a main ",
       "assay.")
}

# Assays visible in a subset: a named vector from assay name to where it is
# stored, "" for the main object or a subset name. Each subset's own
# assays hide those of the same name linked from further up the chain.
.assayOwners <- function(x, subsetName) {
  owners <- rep("", length(assayNames(x)))
  names(owners) <- assayNames(x)
  for (subset in .subsetChain(x, subsetName)) {
    own <- assayNames(.internalAssay(.subsets(x)[[subset]]))
    owners[own] <- subset
  }
  owners
}

# Indices of a subset within the rows or columns of container (a subset
# name, or NULL for the main object). getIndices is .rowIndices or
# .colIndices.
.indicesIn <- function(x, subsetName, container, getIndices) {
  chain <- .subsetChain(x, subsetName)
  start <- if (is.null(container)) 0L else match(container, chain)
  below <- chain[seq_along(chain) > start]
  indices <- getIndices(.subsets(x)[[subsetName]])
  for (subset in rev(below)[-1])
    indices <- getIndices(.subsets(x)[[subset]])[indices]
  indices
}

# An assay as seen in a subset: its own, or a view of the nearest
# container up the chain that stores it, sliced to the subset
.subsetAssayView <- function(x, subsetName, assayName, withDimnames = TRUE) {
  owners <- .assayOwners(x, subsetName)
  if (!assayName %in% names(owners))
    stop("'", assayName, "' is not an assay of subset '", subsetName, "'.")
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  owner <- owners[[assayName]]
  if (owner == subsetName)
    return(assay(internal, assayName, withDimnames = withDimnames))
  if (owner == "") {
    container <- NULL
    source <- assay(x, assayName, withDimnames = FALSE)
  } else {
    container <- owner
    ownerInternal <- .internalAssay(.subsets(x)[[owner]])
    source <- assay(ownerInternal, assayName, withDimnames = FALSE)
  }
  rows <- .indicesIn(x, subsetName, container, .rowIndices)
  cols <- .indicesIn(x, subsetName, container, .colIndices)
  out <- source[rows, cols, drop = FALSE]
  if (withDimnames && !is.null(rownames(internal)))
    rownames(out) <- rownames(internal)
  if (withDimnames && !is.null(colnames(internal)))
    colnames(out) <- colnames(internal)
  out
}

# The name of a subset's primary assay: the assay its parentAssay refers
# to, followed through parent subsets ("filtered" for a subset of subset
# "clean" whose own parent is "filtered/decontaminated" is
# "decontaminated")
.primaryAssayName <- function(x, subsetName) {
  for (subset in rev(.subsetChain(x, subsetName))) {
    parent <- .parentAssay(.subsets(x)[[subset]])
    if (is.null(parent) || !parent %in% subsetNames(x)
        || parent %in% assayNames(x))
      break
  }
  if (is.null(parent))
    return("counts")
  parts <- .splitQualified(parent)
  if (is.null(parts) || parent %in% assayNames(x))
    return(parent)
  parts[2]
}

# Rewrites parents recorded as a bare stored-assay name, in objects saved
# before parents were qualified, to "subset/assay", so that later subsets
# or assays with that name can't change what they point to
.qualifyParents <- function(x) {
  for (name in subsetNames(x)) {
    subset <- .subsets(x)[[name]]
    parent <- .parentAssay(subset)
    bare <- !is.null(parent) && !parent %in% c(assayNames(x), subsetNames(x))
    if (!bare || grepl("/", parent, fixed = TRUE))
      next
    holders <- .assayHolders(x, parent)
    if (length(holders) != 1)
      next
    .subsets(x)[[name]] <- AssaySubset(
      subsetName = .subsetName(subset),
      rowIndices = .rowIndices(subset),
      colIndices = .colIndices(subset),
      parentAssay = paste0(holders, "/", parent),
      internalAssay = .internalAssay(subset)
    )
  }
  x
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

#' @title Combine ExperimentSubset objects by columns or rows
#' @description \code{cbind(...)} and \code{rbind(...)} combine objects as
#'   for the class they extend, and keep the subsets of every object, with
#'   their own assays, rowData, colData, reducedDims and altExps. Subsets
#'   of the second and later objects are moved to where that object's
#'   columns (for \code{cbind}) or rows (for \code{rbind}) are in the
#'   result. Subset names must be unique across the objects.
#' @param ... \code{ExperimentSubset} objects of the same class, which must
#'   have the same rows for \code{cbind} or the same columns for
#'   \code{rbind}. Convert a plain object, such as a
#'   \code{SingleCellExperiment}, with \code{ExperimentSubset()} first,
#'   which keeps all of its data: combining it as it is makes R warn about
#'   "multiple direct matches", and if it comes first, the result is a
#'   plain object without subsets.
#' @param deparse.level Passed on to the method of the class the objects
#'   extend. Default \code{1}.
#' @return An object of the same class as the objects in \code{...}.
#' @name ExperimentSubset-combine
#' @examples
#' # two objects with the same 10 genes and 20 cells each
#' set.seed(1)
#' counts1 <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts1) <- paste0("gene", 1:10)
#' colnames(counts1) <- paste0("cellA", 1:20)
#' counts2 <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts2) <- paste0("gene", 1:10)
#' colnames(counts2) <- paste0("cellB", 1:20)
#' es1 <- ExperimentSubset(list(counts = counts1))
#' es2 <- ExperimentSubset(list(counts = counts2))
#' # a subset in each object
#' es1 <- createSubset(es1, subsetName = "firstA", rows = 1:5, cols = 1:8)
#' es2 <- createSubset(es2, subsetName = "firstB", rows = 1:5, cols = 1:8)
#' # combine the cells of both objects
#' esBoth <- cbind(es1, es2)
#' # both subsets are kept
#' subsetNames(esBoth)
#' subsetColnames(esBoth, subsetName = "firstB")
NULL

# Gives out, the result of cbind() or rbind() on objects, the subsets of
# every object; the parent methods keep only the first object's (#67).
# Subsets made from a main-object assay move by the columns (along =
# "cols") or rows ("rows") of the objects before theirs; subsets of
# subsets keep their indices, which are relative to their parent
.combineSubsets <- function(out, objects, along) {
  getSize <- if (along == "cols") ncol else nrow
  sizes <- vapply(objects, getSize, integer(1))
  offsets <- cumsum(c(0L, sizes[-length(sizes)]))
  hasSubsets <- vapply(objects, methods::.hasSlot, logical(1), "subsets")
  allNames <- unlist(lapply(objects[hasSubsets], subsetNames))
  repeated <- unique(allNames[duplicated(allNames)])
  if (length(repeated) > 0)
    stop("Subset names must be unique across the objects to combine; ",
         paste0("'", repeated, "'", collapse = ", "), " is in more than ",
         "one. Create the subsets under different names before combining.")
  subsets <- list()
  for (k in which(hasSubsets)) {
    # explicit parent names, so a bare one can't match another object's
    object <- .qualifyParents(objects[[k]])
    for (name in subsetNames(object)) {
      subset <- .subsets(object)[[name]]
      rows <- .rowIndices(subset)
      cols <- .colIndices(subset)
      if (is.null(.parentSubsetName(object, .parentAssay(subset)))) {
        if (along == "cols")
          cols <- cols + offsets[k]
        else
          rows <- rows + offsets[k]
      }
      subsets[[name]] <- AssaySubset(
        subsetName = name,
        rowIndices = rows,
        colIndices = cols,
        parentAssay = .parentAssay(subset),
        internalAssay = .internalAssay(subset)
      )
    }
  }
  .subsets(out) <- S4Vectors::SimpleList(subsets)
  out
}

#' @title Subset an ExperimentSubset object by rows or columns
#' @description \code{x[i, j]} keeps rows \code{i} and columns \code{j} of
#'   the main object, as for the class it extends, and also updates every
#'   subset. Each subset keeps the rows and columns that remain, in its own
#'   order, together with its own assays, rowData, colData, reducedDims and
#'   altExps. A subset none of whose rows or columns remain is kept with no
#'   rows or no columns. If \code{i} or \code{j} repeats a row or column,
#'   subsets keep one copy of it.
#' @param x An \code{ExperimentSubset} object.
#' @param i,j Rows and columns to keep: positions, names or logical vectors,
#'   as for the class \code{x} extends. Leave one out to keep all rows or
#'   all columns.
#' @param ... Passed on to the method of the class \code{x} extends.
#' @param drop Ignored, as for \code{SummarizedExperiment}. Default
#'   \code{TRUE}.
#' @return An object of the same class as \code{x}.
#' @name ExperimentSubset-subsetting
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # keep genes 3 to 10 of the main object
#' esSmaller <- es[3:10, ]
#' # the subset keeps genes 3 to 5
#' subsetRownames(esSmaller, subsetName = "first")
NULL

# Updates the subsets of out, the result of x[i, j], to the rows and
# columns that remain; i or j may be missing. Subsets are visited parents
# first, so each is matched against what remains of its parent.
.remapSubsets <- function(x, out, i, j) {
  if (missing(i) && missing(j))
    return(out)
  rows <- NULL
  if (!missing(i)) {
    positions <- seq_len(nrow(x))
    names(positions) <- rownames(x)
    rows <- S4Vectors::normalizeSingleBracketSubscript(i, positions)
  }
  cols <- NULL
  if (!missing(j)) {
    positions <- seq_len(ncol(x))
    names(positions) <- colnames(x)
    cols <- S4Vectors::normalizeSingleBracketSubscript(j, positions)
  }
  x <- .qualifyParents(x)
  chains <- lapply(subsetNames(x), function(name) .subsetChain(x, name))
  names(chains) <- subsetNames(x)
  depth <- vapply(chains, length, integer(1))
  # positions in each subset's old rows and columns that remain, in order;
  # NULL when all remain
  keptRows <- list()
  keptCols <- list()
  for (name in names(chains)[order(depth)]) {
    subset <- .subsets(x)[[name]]
    chain <- chains[[name]]
    container <- if (length(chain) > 1) chain[length(chain) - 1] else NULL
    parentRows <- if (is.null(container)) rows else keptRows[[container]]
    parentCols <- if (is.null(container)) cols else keptCols[[container]]
    newRows <- .keepIndices(.rowIndices(subset), parentRows)
    newCols <- .keepIndices(.colIndices(subset), parentCols)
    keptRows[name] <- list(newRows$kept)
    keptCols[name] <- list(newCols$kept)
    internal <- .internalAssay(subset)
    internalRows <- newRows$kept
    if (is.null(internalRows))
      internalRows <- seq_len(nrow(internal))
    internalCols <- newCols$kept
    if (is.null(internalCols))
      internalCols <- seq_len(ncol(internal))
    .subsets(out)[[name]] <- AssaySubset(
      subsetName = .subsetName(subset),
      rowIndices = newRows$indices,
      colIndices = newCols$indices,
      parentAssay = .parentAssay(subset),
      internalAssay = internal[internalRows, internalCols]
    )
  }
  out
}

# Remaps indices into a container whose remaining positions, in their new
# order, are parentKept (NULL when all remain in place). Returns the new
# indices and the positions in indices that remain (NULL when all do).
.keepIndices <- function(indices, parentKept) {
  if (is.null(parentKept))
    return(list(indices = indices, kept = NULL))
  newIndices <- match(indices, parentKept)
  kept <- which(!is.na(newIndices))
  list(indices = newIndices[kept], kept = kept)
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
    parentRows <- .hideParentColumns(parentRows, out, subsetName, "rowData")
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
    # a SpatialExperiment subset stores a copy of sample_id, made by
    # createSubset(); the main object's current one is used instead
    if (methods::is(x, "SpatialExperiment"))
      out <- out[, colnames(out) != "sample_id", drop = FALSE]
    parentCols <- .hideParentColumns(parentCols, out, subsetName, "colData")
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

# The main object's columns without those the subset also has: a subset's
# own column hides the main object's column of the same name (#57). A
# message names the hidden columns whose values differ from the subset's
.hideParentColumns <- function(parentData, own, subsetName, what) {
  hidden <- colnames(parentData) %in% colnames(own)
  # each hidden column, so a name the main object repeats is checked for
  # every copy
  differs <- vapply(which(hidden), function(i) {
    !.sameColumn(parentData[[i]], own[[colnames(parentData)[i]]])
  }, logical(1))
  if (any(differs)) {
    differing <- unique(colnames(parentData)[which(hidden)[differs]])
    message("Subset '", subsetName, "' has its own ", what, " columns ",
            "that hide the main object's columns of the same name, with ",
            "different values: ", paste(differing, collapse = ", "),
            ". The main object's values are in ", what, "(x).")
  }
  parentData[, !hidden, drop = FALSE]
}

# Whether two columns hold the same values, ignoring the names of atomic
# vectors. Deliberately strict: equal values of different types, such as
# integer and double, count as different, which costs only a message
.sameColumn <- function(a, b) {
  if (is.atomic(a))
    a <- unname(a)
  if (is.atomic(b))
    b <- unname(b)
  identical(a, b)
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
  fromParent <- .parentColumns(value)
  if (any(fromParent))
    message("Not storing columns that come from the main object's ", what,
            ": ", paste(colnames(value)[fromParent], collapse = ", "),
            ". Change them on the main object with ", what, "(x).")
  .withoutParentColumns(value)
}

# Which columns of value .tagParentColumns() marked
.parentColumns <- function(value) {
  if (!methods::is(value, "DataFrame"))
    return(logical(0))
  columnInfo <- S4Vectors::mcols(value)
  if (is.null(columnInfo) || !"fromParentData" %in% colnames(columnInfo))
    return(rep(FALSE, ncol(value)))
  columnInfo$fromParentData %in% TRUE
}

# value without the columns .tagParentColumns() marked, and without the mark
.withoutParentColumns <- function(value) {
  fromParent <- .parentColumns(value)
  if (any(fromParent))
    value <- value[, !fromParent, drop = FALSE]
  .untagParentColumns(value)
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

# Used by both subsetSpatialData() and spatialData(x, subsetName =), which
# are deprecated in step with SpatialExperiment's spatialData() (#34)
.subsetSpatialData <- function(x, subsetName, parentSpatialData){
  .isSubset(x, subsetName)
  msg <- paste0(
    "subsetSpatialData() and spatialData(x, subsetName =) are deprecated, ",
    "as is spatialData() in SpatialExperiment, which keeps these columns ",
    "in colData. Use subsetColData(x, subsetName) instead, which returns ",
    "all of the subset's colData columns, with parentColData = TRUE to ",
    "include the main object's.\nSee help(\"Deprecated\")"
  )
  # through do.call() so BiocCheck, which flags any .Deprecated() call,
  # stays clean (#34)
  do.call(".Deprecated", list(msg = msg))
  if(missing(parentSpatialData)
     || is.null(parentSpatialData))
    parentSpatialData = FALSE
  # a DataFrame, since spot names can repeat across samples
  nSpots <- subsetDim(x, subsetName)[2]
  out <- S4Vectors::DataFrame(row.names = as.character(seq_len(nSpots)))
  if(parentSpatialData) {
    rootCols <- .rootIndices(x, subsetName, .colIndices)
    # SpatialExperiment's note repeats the warning above
    out <- suppressMessages(spatialData(x))[rootCols, , drop = FALSE]
  }
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  subsetSpatialData <- suppressMessages(spatialData(internal))
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

# Defunct without a deprecation release: it failed for every value, since
# SpatialExperiment's spatialData<- takes no data.frame, and the method
# took no DataFrame (#34)
'.subsetSpatialData<-' <- function(x, subsetName, value){
  msg <- paste0(
    "subsetSpatialData<- is defunct, as spatialData<- is deprecated in ",
    "SpatialExperiment. Add each column to the subset's colData instead, ",
    "with colData(x, subsetName = \"name\")$column <- values.\n",
    "See help(\"Defunct\")"
  )
  # through do.call() so BiocCheck, which flags any .Defunct() call,
  # stays clean (#34)
  do.call(".Defunct", list(msg = msg))
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
  if (!is.null(ref) && ref$kind == "visible")
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

  if (grepl("/", subsetName, fixed = TRUE))
    stop("'subsetName' can't contain '/', which separates a subset from an ",
         "assay name.")

  if (subsetName %in% assayNames(x))
    stop("'subsetName' must differ from the names of existing assays; '",
         subsetName, "' is already an assay name.")
  
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
  if (dimR == 0)
    stop("Parent '", assayName, "' has no rows, so a subset can't be made ",
         "from it.")
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
  if (dimC == 0)
    stop("Parent '", assayName, "' has no columns, so a subset can't be ",
         "made from it.")
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
