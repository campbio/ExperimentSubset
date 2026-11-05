
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
#' @param x An experiment object of one of the five supported classes, or a
#'   named \code{list} of assay matrices, which is passed to the
#'   \code{SingleCellExperiment} constructor as its assays. An object of a
#'   subclass of a supported class is converted as the class it extends,
#'   with a warning: the slots and methods that the subclass adds are not
#'   kept. An object that is already an \code{ExperimentSubset} is
#'   returned as is, with the subset from \code{subset}, if given.
#' @param subset A named \code{list} to create a subset from within the
#'   constructor, with elements \code{subsetName} and, optionally,
#'   \code{rows}, \code{cols} and \code{parentAssay}. Leave out
#'   \code{rows}, \code{cols} or \code{parentAssay}, or set them to
#'   \code{NULL}, to use all rows, all columns or the first assay. The
#'   default, a list of \code{NA} values, creates no subset; an \code{NA}
#'   in a list you supply is an error.
#' @return An object of the matching \code{Subset*} class, for example a
#'   \code{SubsetSingleCellExperiment} for a \code{SingleCellExperiment} or
#'   a list of assays. For an \code{ExperimentSubset} object, the same
#'   object, with the new subset if \code{subset} is given.
#' @seealso \code{\linkS4class{ExperimentSubset}}, the class every
#'   \code{Subset*} class extends, for testing whether an object is an
#'   \code{ExperimentSubset}.
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
  if (methods::is(x, "ExperimentSubset")) {
    es <- x
  } else if (!is.list(x)) {
    baseClass <- .baseClass(x)
    if (is.na(baseClass))
      stop("'x' must be a SummarizedExperiment, RangedSummarizedExperiment, ",
           "SingleCellExperiment, TreeSummarizedExperiment or ",
           "SpatialExperiment object, or a named list of assays; it is an ",
           "object of class ", class(x)[1], ".")
    .checkAssayNames(x)
    if (!identical(as.character(class(x)), baseClass)) {
      warning("'x' is an object of class ", class(x)[1], ", a subclass of ",
              baseClass,
              ", so it is converted as a ", baseClass, ": the result is a ",
              "Subset", baseClass, ", without the slots and methods that ",
              class(x)[1], " adds.")
      x <- methods::as(x, baseClass, strict = TRUE)
    }
    es <- as(x, paste0("Subset", baseClass))
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

#' @title Metadata of subsets in ExperimentSubset objects
#' @description Get or set the \code{metadata} list of a subset by naming
#'   it with \code{subsetName}. Without \code{subsetName}, \code{metadata()}
#'   and \code{metadata<-} act on the main object, as for the class
#'   \code{x} extends.
#' @param x An \code{ExperimentSubset} object.
#' @param ... \describe{
#'     \item{\code{subsetName}}{\code{character(1)} Name of the subset whose
#'       metadata to get or set. Leave it out for the main object's.}
#'   }
#'   Without \code{subsetName}, all arguments are passed on to the method
#'   of the class \code{x} extends.
#' @param value A \code{list} to store as the metadata.
#' @return \code{metadata()} returns a \code{list}. \code{metadata<-}
#'   returns \code{x} with the metadata set.
#' @seealso \code{\link[S4Vectors]{metadata}}
#' @name ExperimentSubset-metadata
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # store a note in the subset's metadata
#' metadata(es, subsetName = "first") <- list(note = "the first 8 cells")
#' # read it back
#' metadata(es, subsetName = "first")
#' # the main object's metadata is unchanged
#' metadata(es)
NULL

#' @importMethodsFrom S4Vectors metadata
.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @title Assays of subsets in ExperimentSubset objects
#' @description Get or set an assay of the main object or of a subset.
#'   Unlike the other methods, \code{assay()} and \code{assay<-} name the
#'   subset with \code{i}, not with \code{subsetName}.
#'
#'   \code{assay(x, i)} returns:
#'   \itemize{
#'     \item a main assay, when \code{i} is its name, as for the class
#'       \code{x} extends;
#'     \item the subset's rows and columns of the assay it was made from,
#'       when \code{i} is a subset name. If the subset, or a subset it was
#'       made from, stored its own assay of that name with
#'       \code{unlink}, that assay is returned instead;
#'     \item an assay as a subset sees it, when \code{i} is
#'       \code{"subset/assay"}, or a subset name with
#'       \code{subsetAssayName = "assay"}: either one stored in the subset,
#'       or one linked from its parent, cut to the subset's rows and
#'       columns.
#'   }
#'   \code{assay(x, i) <- value} stores \code{value} in a subset when
#'   \code{i} is \code{"subset/assay"}, or a subset name with
#'   \code{subsetAssayName}, as \code{\link{setSubsetAssay}} does. A name
#'   with a \code{"/"} is read as \code{"subset/assay"} unless it is the
#'   name of a main assay, and a subset name without
#'   \code{subsetAssayName} is an error. Any other name stores a main
#'   assay.
#' @param x An \code{ExperimentSubset} object.
#' @param i \code{character(1)} The name of a main assay, a subset name, or
#'   \code{"subset/assay"}.
#' @param withDimnames \code{logical(1)} As for \code{SummarizedExperiment}:
#'   \code{TRUE} gives the assay returned the row and column names of the
#'   main object or the subset, and \code{FALSE} keeps the names stored
#'   with the assay, if any. When setting a main assay, passed on to the
#'   method of the class \code{x} extends. Default \code{TRUE}.
#' @param ... \describe{
#'     \item{\code{subsetAssayName}}{\code{character(1)} With a subset name
#'       in \code{i}, the name of the assay to get from or store in that
#'       subset. It is an error with any other \code{i}, including
#'       \code{"subset/assay"}. Default \code{NULL}.}
#'     \item{\code{unlink}}{\code{TRUE}, \code{FALSE} or a character
#'       vector of assay names. When storing in a subset, store the assay
#'       even though the subset sees an assay of that name linked from its
#'       parent, with \code{TRUE} or a vector that contains its name, as for
#'       \code{\link{setSubsetAssay}}. It is an error when storing a main
#'       assay. Default \code{FALSE}.}
#'   }
#'   For a main assay, other arguments are passed on to the method of the
#'   class \code{x} extends.
#' @param value A matrix-like assay with the main object's dimensions, or,
#'   to store it in a subset, the subset's.
#' @return \code{assay()} returns a matrix-like assay. \code{assay<-}
#'   returns \code{x} with the assay stored.
#' @seealso \code{\link{getSubsetAssay}}, \code{\link{setSubsetAssay}},
#'   \code{\link[SummarizedExperiment]{assay}}
#' @name ExperimentSubset-assay
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # the subset's 5 x 8 part of the counts
#' firstCounts <- assay(es, "first")
#' # store a new assay in the subset
#' assay(es, "first", subsetAssayName = "logcounts") <- log1p(firstCounts)
#' # read it back
#' assay(es, "first/logcounts")
NULL

# Why assay() or assay<- can't take subsetAssayName with i, which is not a
# subset name: i already names an assay as "subset/assay" (qualified), or
# it names something other than a subset. The caller says which, from the
# test it already made. Shared so both give the same reason (#76)
.subsetAssayNameMisuse <- function(i, qualified) {
  if (qualified)
    return(paste0("'subsetAssayName' can't be used with \"subset/assay\"; '",
                  i, "' already names the assay."))
  paste0("'subsetAssayName' is only used when 'i' is a subset name; '", i,
         "' is not a subset.")
}

# The two names in i, "subset/assay"; stops unless the subset exists.
# Used by assay() and assay<- so both say the same (#76)
.qualifiedParts <- function(x, i) {
  parts <- .splitQualified(i)
  if (is.null(parts) || !parts[1] %in% subsetNames(x))
    stop("'", i, "' is not \"subset/assay\" for an existing subset.")
  parts
}

.assay <- function(x, i, ...){
  .checkAssayI(i)
  args <- list(...)
  withDimnames <- !isFALSE(args$withDimnames)
  # assay(x, "subset", subsetAssayName = "a") means "subset/a"
  if (!is.null(args$subsetAssayName)) {
    # NA would become "subset/NA" (#85)
    if (!.isOneString(args$subsetAssayName))
      stop("'subsetAssayName' must be a single character string.")
    if (!i %in% subsetNames(x)) {
      qualified <- grepl("/", i, fixed = TRUE) && !i %in% assayNames(x)
      # as for the setter, a "/" name must first be "subset/assay" for an
      # existing subset
      if (qualified)
        .qualifiedParts(x, i)
      stop(.subsetAssayNameMisuse(i, qualified))
    }
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
  #its own if it, or a subset on the way, stored one with unlink
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
  internalAssay <- as(internalAssay, .baseClass(x))
  
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
                         unlink = FALSE,
                         quiet = FALSE) {

  .checkSubsetName(subsetName)
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
  .checkUnlink(unlink)
  # a name linked from the parent is protected unless unlink allows it
  source <- .checkLinked(.assayOwners(x, subsetName), subsetName,
                         subsetAssayName, unlink)

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
  if (length(source) > 0 && !quiet)
    .unlinkMessage(subsetName, subsetAssayName, source)

  return(x)
}

# Stops unless unlink is TRUE, FALSE or a character vector of assay names
.checkUnlink <- function(unlink) {
  validLogical <- is.logical(unlink) && length(unlink) == 1 && !is.na(unlink)
  validNames <- is.character(unlink) && !anyNA(unlink) && all(nzchar(unlink))
  if (!validLogical && !validNames)
    stop("'unlink' must be TRUE, FALSE or a character vector of assay ",
         "names.")
  invisible(NULL)
}

# TRUE if unlink allows storing the assay assayName over a linked one
.unlinks <- function(unlink, assayName) {
  isTRUE(unlink) || (is.character(unlink) && assayName %in% unlink)
}

# Says that subsetName now stores its own assayName, linked from source
.unlinkMessage <- function(subsetName, assayName, source) {
  message("Subset '", subsetName, "' now stores its own '", assayName,
          "', unlinked from ", source, ".")
}

# For the assays in assayNames that subsetName sees linked from its
# parent (owners is .assayOwners(x, subsetName)), stops unless unlink
# allows storing its own, naming every one it doesn't, and returns where
# each is linked from, as .linkedSources() does
.checkLinked <- function(owners, subsetName, assayNames, unlink) {
  .checkUnlinked(.linkedSources(owners, subsetName, assayNames),
                 subsetName, unlink)
}

# Where each of the assays in assayNames that subsetName sees linked from
# its parent is linked from (for example "the main assay 'counts'"), named
# by assay; character(0) when none is linked
.linkedSources <- function(owners, subsetName, assayNames) {
  assayNames <- assayNames[assayNames %in% names(owners)]
  assayNames <- assayNames[owners[assayNames] != subsetName]
  if (length(assayNames) == 0)
    return(character(0))
  sources <- paste0("'", owners[assayNames], "/", assayNames, "'")
  main <- owners[assayNames] == ""
  sources[main] <- paste0("the main assay '", assayNames[main], "'")
  names(sources) <- assayNames
  sources
}

# Stops unless unlink allows storing every linked assay in sources (named
# by assay, or "altExp/assay" for an assay of an altExp, as from
# .linkedSources()), naming every one it doesn't, and returns sources
.checkUnlinked <- function(sources, subsetName, unlink) {
  if (length(sources) == 0)
    return(character(0))
  assayNames <- names(sources)
  blocked <- assayNames[!vapply(assayNames, .unlinks, logical(1),
                                unlink = unlink)]
  if (length(blocked) == 0)
    return(sources)
  # with names in unlink, suggest adding these to them
  suggested <- blocked
  if (is.character(unlink))
    suggested <- c(unname(unlink), blocked)
  suggested <- paste(deparse(suggested), collapse = "")
  quoted <- paste0("'", blocked, "'")
  linkedTo <- paste0(quoted, " to ", sources[blocked])
  linkedTo[1] <- paste0(quoted[1], " in subset '", subsetName,
                        "' is linked to ", sources[[blocked[1]]])
  separate <- if (length(blocked) == 1) "a separate" else "separate"
  what <- paste(separate, .andList(quoted))
  stop(.andList(linkedTo), ". Set unlink = TRUE or unlink = ", suggested,
       " to store ", what, " computed on the cells of '", subsetName,
       "'.")
}

# Joins words as "a", "a and b" or "a, b and c"
.andList <- function(words) {
  n <- length(words)
  if (n < 2)
    return(paste(words, collapse = ""))
  paste(paste(words[-n], collapse = ", "), "and", words[n])
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
    # every altExp the subset sees, its own and linked
    altExps(out, withDimnames = FALSE) <-
      altExps(.withSubsetAltExps(x, subsetName), withDimnames = FALSE)
  }
  metadata(out) <- metadata(internal)
  out
}

# Stores what a tool added to an object from .getSubsetExperiment(); value
# replaces the subset's own data, and unchanged linked assays and the main
# object's tagged columns are left out
.setSubsetExperiment <- function(x, subsetName, value, unlink = FALSE) {
  .isSubset(x, subsetName)
  .checkUnlink(unlink)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  .checkSubsetValue(internal, subsetName, value)

  # unchanged linked assays are not stored; changed ones only if unlink
  # allows it, checked before anything is stored
  x <- .qualifyParents(x)
  owners <- .assayOwners(x, subsetName)
  linked <- names(owners)[owners != subsetName]
  stored <- character(0)
  for (name in assayNames(value)) {
    if (!name %in% linked
        || !.sameAssay(assay(value, name),
                       .subsetAssayView(x, subsetName, name)))
      stored <- c(stored, name)
  }
  # the same for the altExps' assays; one check for both, so one error
  # names every assay unlink doesn't allow
  isSce <- methods::is(internal, "SingleCellExperiment")
  if (isSce)
    altExpPlan <- .planAltExps(x, subsetName, altExps(value))
  # where each unlinked assay was linked from, for the messages once
  # everything is stored
  sources <- .linkedSources(owners, subsetName, stored)
  if (isSce)
    sources <- c(sources, altExpPlan$sources)
  sources <- .checkUnlinked(sources, subsetName, unlink)

  # own assays that value no longer has are removed, unless a subset was
  # made from them
  kept <- assayNames(internal) %in% assayNames(value)
  .checkNotParents(x, subsetName, assayNames(internal)[!kept])
  assays(internal, withDimnames = FALSE) <-
    assays(internal, withDimnames = FALSE)[kept]
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  for (name in stored)
    x <- .setSubsetAssay(x, subsetName, assay(value, name), name,
                         unlink = unlink, quiet = TRUE)

  .subsetRowData(x, subsetName) <- .withoutParentColumns(rowData(value))
  .subsetColData(x, subsetName) <- .withoutParentColumns(colData(value))
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  if (isSce)
    SingleCellExperiment::reducedDims(internal) <- reducedDims(value)
  metadata(internal) <- metadata(value)
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  if (isSce)
    x <- .applyAltExps(x, subsetName, altExpPlan$parts)
  for (name in names(sources))
    .unlinkMessage(subsetName, name, sources[[name]])
  x
}

# Stops unless value, for setSubsetExperiment(), has the class, dimensions
# and row and column names of subsetName's internal object
.checkSubsetValue <- function(internal, subsetName, value) {
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
  invisible(NULL)
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
  .checkAssayI(i)
  unlink <- list(...)[["unlink"]]
  if (is.null(unlink))
    unlink <- FALSE
  .checkUnlink(unlink)
  # "subset/assay" stores into that subset
  if (grepl("/", i, fixed = TRUE) && !i %in% assayNames(x)) {
    parts <- .qualifiedParts(x, i)
    if (!is.null(subsetAssayName))
      stop(.subsetAssayNameMisuse(i, qualified = TRUE))
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
      stop(.subsetAssayNameMisuse(i, qualified = FALSE))
    if (!isFALSE(unlink) && length(unlink) > 0)
      stop("'unlink' only applies when storing an assay in a subset.")
    callNextMethod(x, i, ..., value = value)
  }
}

.subsetParent <- function(x, subsetName){
  .checkSubsetName(subsetName)
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

#' @title Alternative experiments of subsets in ExperimentSubset objects
#' @description Get or set the alternative experiments (\code{altExps}) of
#'   a subset by naming it with \code{subsetName}. Without
#'   \code{subsetName}, these methods act on the main object, as for
#'   \code{SingleCellExperiment}. They exist for the classes that extend
#'   \code{SingleCellExperiment}: \code{SubsetSingleCellExperiment},
#'   \code{SubsetSpatialExperiment} and \code{SubsetTreeSummarizedExperiment}.
#' @param x An \code{ExperimentSubset} object.
#' @param e \code{character(1)} or \code{integer(1)} The name or position
#'   of the alternative experiment. Leave it out for the first one.
#' @param withDimnames \code{logical(1)} Passed on to
#'   \code{\link[SingleCellExperiment]{altExp<-}} or \code{altExps<-}.
#'   Default \code{TRUE}.
#' @param withColData \code{logical(1)} Passed on to
#'   \code{\link[SingleCellExperiment]{altExp<-}} or \code{altExps<-}.
#'   Default \code{FALSE}.
#' @param ... \describe{
#'     \item{\code{subsetName}}{\code{character(1)} Name of the subset whose
#'       alternative experiments to get or set. Leave it out for the main
#'       object's.}
#'   }
#'   Other arguments are passed on to the \code{SingleCellExperiment}
#'   method. With \code{subsetName}, \code{altExpNames()} and
#'   \code{altExpNames<-} ignore other arguments.
#' @param value For \code{altExp<-}, a \code{SummarizedExperiment} with one
#'   column per column of the subset (or of the main object). For
#'   \code{altExps<-}, a list of them. For \code{altExpNames<-}, a
#'   \code{character} vector of new names.
#' @return \code{altExp()} returns a \code{SummarizedExperiment},
#'   \code{altExps()} a list of them and \code{altExpNames()} a
#'   \code{character} vector. The setters return \code{x} with the change
#'   made.
#' @seealso \code{\link[SingleCellExperiment]{altExps}}
#' @name ExperimentSubset-altExp
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # counts of 2 spike-in genes in the subset's 8 cells
#' spikeCounts <- matrix(rpois(16, lambda = 5), nrow = 2)
#' spikes <- SingleCellExperiment(list(counts = spikeCounts))
#' # store them as an alternative experiment of the subset
#' altExp(es, "spikes", subsetName = "first") <- spikes
#' # list the subset's alternative experiments
#' altExpNames(es, subsetName = "first")
#' # read it back
#' altExp(es, "spikes", subsetName = "first")
NULL

#' @importMethodsFrom SingleCellExperiment altExp
.altExp <- function(x, e, ...) {
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist)) {
    if (missing(e))
      return(callNextMethod(x, ...))
    return(callNextMethod(x, e, ...))
  }
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .withSubsetAltExps(x, subsetName)
  # a missing e can't be passed into a second S4 dispatch
  if (missing(e))
    return(do.call(altExp, c(list(internal), .otherArgs(...))))
  do.call(altExp, c(list(internal, e), .otherArgs(...)))
}

#' @importMethodsFrom SingleCellExperiment altExps
.altExps <- function(x, ...) {
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .withSubsetAltExps(x, subsetName)
  do.call(altExps, c(list(internal), .otherArgs(...)))
}

#' @importMethodsFrom SingleCellExperiment altExpNames
.altExpNames <- function(x, ...) {
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ...))
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  .altExpHolders(x, subsetName)$names
}

#' @importMethodsFrom SingleCellExperiment altExpNames<-
'.altExpNames<-' <- function(x, ..., value) {
  arglist <- list(...)
  if (!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  visible <- .altExpHolders(x, subsetName)$names
  if (length(value) != length(visible))
    stop("'value' must give one name for each of the ", length(visible),
         " altExps of subset '", subsetName, "'.")
  linked <- .altExpHolders(x, subsetName, own = FALSE)
  for (i in which(visible %in% linked$names & visible != value)) {
    holders <- linked$holders[[match(visible[i], linked$names)]]
    stop("Subset '", subsetName, "' can't rename '", visible[i],
         "', which is linked to ",
         .altExpSource(visible[i], holders[1]), ".")
  }
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  altExpNames(internal) <- value[match(altExpNames(internal), visible)]
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  x
}

#' @importMethodsFrom SingleCellExperiment altExp<-
'.altExp<-' <- function(x, e, ..., value) {
  arglist <- list(...)
  unlink <- .unlinkArg(arglist)
  if (!"subsetName" %in% names(arglist)) {
    if (missing(e))
      return(callNextMethod(x, ..., value = value))
    return(callNextMethod(x, e, ..., value = value))
  }
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  # the SingleCellExperiment method checks value and handles e,
  # withDimnames and withColData on everything the subset sees
  internal <- .withSubsetAltExps(x, subsetName)
  if (missing(e))
    args <- c(list(internal), .otherArgs(...), list(value = value))
  else
    args <- c(list(internal, e), .otherArgs(...), list(value = value))
  internal <- do.call(`altExp<-`, args)
  .storeAltExps(x, subsetName, altExps(internal, withDimnames = FALSE),
                unlink)
}

#' @importMethodsFrom SingleCellExperiment altExps<-
'.altExps<-' <- function(x, ..., value) {
  arglist <- list(...)
  unlink <- .unlinkArg(arglist)
  if (!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .withSubsetAltExps(x, subsetName)
  args <- c(list(internal), .otherArgs(...), list(value = value))
  internal <- do.call(`altExps<-`, args)
  .storeAltExps(x, subsetName, altExps(internal, withDimnames = FALSE),
                unlink)
}

# The unlink argument of the altExp setters, checked; FALSE if not given.
# It only applies with subsetName.
.unlinkArg <- function(arglist) {
  unlink <- arglist[["unlink"]]
  if (is.null(unlink))
    return(FALSE)
  .checkUnlink(unlink)
  if (!"subsetName" %in% names(arglist) && !isFALSE(unlink)
      && length(unlink) > 0)
    stop("'unlink' only applies when storing in a subset.")
  unlink
}

# The altExps a subset sees: their names, the main object's first, then
# those stored only on subsets down the chain; and for each, where it is
# stored, top first: "" for the main object, or a subset name. With
# own = FALSE, what the subset sees linked, without its own altExps.
#' @importMethodsFrom SingleCellExperiment altExpNames
.altExpHolders <- function(x, subsetName, own = TRUE) {
  chain <- .subsetChain(x, subsetName)
  if (!own)
    chain <- chain[-length(chain)]
  names <- altExpNames(x)
  holders <- rep(list(""), length(names))
  for (subset in chain) {
    for (name in altExpNames(.internalAssay(.subsets(x)[[subset]]))) {
      i <- match(name, names)
      if (is.na(i)) {
        names <- c(names, name)
        holders <- c(holders, list(subset))
      } else {
        holders[[i]] <- c(holders[[i]], subset)
      }
    }
  }
  list(names = names, holders = holders)
}

# How an error or message names where an altExp is stored
.altExpSource <- function(name, holder) {
  if (holder == "")
    return(paste0("the main altExp '", name, "'"))
  paste0("altExp '", name, "' of subset '", holder, "'")
}

# The altExp called name as subsetName sees it, built from the parts
# stored in holders (from .altExpHolders()), top first, each cut to the
# subset's columns. A part of the same class, rows and rownames adds its
# assays, each hiding one of the same name, and replaces rowData, colData,
# metadata and reducedDims that it has; any other part replaces the
# altExp. Returns the altExp and, named by assay, where each is stored.
.subsetAltExpView <- function(x, subsetName, name, holders) {
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  view <- NULL
  owners <- character(0)
  for (holder in holders) {
    if (holder == "") {
      part <- .storedAltExp(x, name)
      part <- part[, .indicesIn(x, subsetName, NULL, .colIndices)]
    } else {
      part <- .storedAltExp(.internalAssay(.subsets(x)[[holder]]), name)
      if (holder != subsetName)
        part <- part[, .indicesIn(x, subsetName, holder, .colIndices)]
    }
    partOwners <- rep(holder, length(assayNames(part)))
    names(partOwners) <- assayNames(part)
    if (is.null(view) || !.sameAltExpShape(view, part)) {
      view <- part
      owners <- partOwners
      next
    }
    view <- .mergeAltExpPart(view, part)
    owners[names(partOwners)] <- partOwners
  }
  if (!is.null(colnames(internal)))
    colnames(view) <- colnames(internal)
  list(altExp = view, owners = owners)
}

# The altExp called name as stored in object, without its dimnames set
#' @importMethodsFrom SingleCellExperiment altExps
.storedAltExp <- function(object, name) {
  altExps(object, withDimnames = FALSE)[[match(name, altExpNames(object))]]
}

# TRUE if a subset's part of an altExp adds to the linked one rather than
# replacing it: the same class, number of rows and rownames
.sameAltExpShape <- function(linked, part) {
  identical(class(linked), class(part)) && nrow(linked) == nrow(part) &&
    identical(rownames(linked), rownames(part))
}

# The linked altExp with a subset's part of the same shape merged in
.mergeAltExpPart <- function(view, part) {
  for (name in assayNames(part))
    assay(view, name, withDimnames = FALSE) <-
      assay(part, name, withDimnames = FALSE)
  if (ncol(rowData(part)) > 0)
    rowData(view) <- rowData(part)
  if (ncol(colData(part)) > 0)
    colData(view) <- colData(part)
  if (length(metadata(part)) > 0)
    metadata(view) <- metadata(part)
  if (methods::is(part, "SingleCellExperiment")
      && length(SingleCellExperiment::reducedDims(part)) > 0)
    SingleCellExperiment::reducedDims(view, withDimnames = FALSE) <-
      SingleCellExperiment::reducedDims(part, withDimnames = FALSE)
  view
}

# A copy of subsetName's internal object holding every altExp the subset
# sees, for the SingleCellExperiment altExp methods to run on
.withSubsetAltExps <- function(x, subsetName) {
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  visible <- .altExpHolders(x, subsetName)
  views <- lapply(seq_along(visible$names), function(i) {
    .subsetAltExpView(x, subsetName, visible$names[i],
                      visible$holders[[i]])$altExp
  })
  names(views) <- visible$names
  SingleCellExperiment::altExps(internal, withDimnames = FALSE) <- views
  internal
}

# Stores value, every altExp subsetName should see, as the subset's own
# altExps, with the unlink rule for changed linked assays, and says which
# assays are unlinked
.storeAltExps <- function(x, subsetName, value, unlink) {
  plan <- .planAltExps(x, subsetName, value)
  sources <- .checkUnlinked(plan$sources, subsetName, unlink)
  x <- .applyAltExps(x, subsetName, plan$parts)
  for (name in names(sources))
    .unlinkMessage(subsetName, name, sources[[name]])
  x
}

# What subsetName stores to see the altExps in value (a named list): for
# each, NULL if the linked altExp is the same, the value whole if nothing
# is linked under its name or its shape differs, or else a part with only
# the assays, rowData, colData, metadata and reducedDims that differ.
# Also, named "altExp/assay", where each changed linked assay is linked
# from, for the unlink check.
.planAltExps <- function(x, subsetName, value) {
  linked <- .altExpHolders(x, subsetName, own = FALSE)
  parts <- vector("list", length(value))
  names(parts) <- names(value)
  sources <- character(0)
  for (i in seq_along(value)) {
    name <- names(value)[i]
    new <- value[[i]]
    j <- match(name, linked$names)
    if (is.na(j)) {
      parts[[i]] <- new
      next
    }
    old <- .subsetAltExpView(x, subsetName, name, linked$holders[[j]])
    changed <- names(old$owners)
    if (.sameAltExpShape(old$altExp, new)) {
      changed <- .changedAssays(old$altExp, new)
      parts[i] <- list(.altExpPart(old$altExp, new))
    } else {
      parts[[i]] <- new
    }
    changed <- intersect(changed, names(old$owners))
    if (length(changed) > 0) {
      changedSources <- vapply(old$owners[changed], .altExpSource,
                               character(1), name = name)
      names(changedSources) <- paste0(name, "/", changed)
      sources <- c(sources, changedSources)
    }
  }
  list(parts = parts, sources = sources)
}

# Names of the assays of new that old doesn't have or holds other values
# for
.changedAssays <- function(old, new) {
  changed <- character(0)
  for (name in assayNames(new)) {
    if (!name %in% assayNames(old)
        || !.sameAssay(assay(new, name, withDimnames = FALSE),
                       assay(old, name, withDimnames = FALSE)))
      changed <- c(changed, name)
  }
  changed
}

# new, of the same shape as the linked altExp old, without what old
# already has: unchanged assays, and rowData, colData, metadata and
# reducedDims that are the same. NULL if nothing differs.
.altExpPart <- function(old, new) {
  part <- new
  own <- .changedAssays(old, new)
  assays(part, withDimnames = FALSE) <-
    assays(new, withDimnames = FALSE)[own]
  if (.sameColumns(rowData(new), rowData(old)))
    rowData(part) <- rowData(part)[, 0]
  if (.sameColumns(colData(new), colData(old)))
    colData(part) <- colData(part)[, 0]
  if (identical(metadata(new), metadata(old)))
    metadata(part) <- list()
  hasReducedDims <- methods::is(part, "SingleCellExperiment")
  if (hasReducedDims && .sameAssayList(
    SingleCellExperiment::reducedDims(new, withDimnames = FALSE),
    SingleCellExperiment::reducedDims(old, withDimnames = FALSE))) {
    SingleCellExperiment::reducedDims(part) <- list()
    hasReducedDims <- FALSE
  }
  if (length(own) == 0 && ncol(rowData(part)) == 0
      && ncol(colData(part)) == 0 && length(metadata(part)) == 0
      && !hasReducedDims)
    return(NULL)
  part
}

# TRUE if two DataFrames have the same columns, whatever their rownames
.sameColumns <- function(a, b) {
  identical(colnames(a), colnames(b)) &&
    identical(lapply(as.list(a), unname), lapply(as.list(b), unname))
}

# TRUE if two lists of matrices have the same names and values, as
# .sameAssay() compares them
.sameAssayList <- function(a, b) {
  a <- as.list(a)
  b <- as.list(b)
  if (!identical(names(a), names(b)) || length(a) != length(b))
    return(FALSE)
  all(vapply(seq_along(a), function(i) .sameAssay(a[[i]], b[[i]]),
             logical(1)))
}

# Stores parts (from .planAltExps()) as subsetName's own altExps, in place
# of those it had; NULL parts are left out
.applyAltExps <- function(x, subsetName, parts) {
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  parts <- parts[!vapply(parts, is.null, logical(1))]
  SingleCellExperiment::altExps(internal, withDimnames = FALSE) <- parts
  .internalAssay(.subsets(x)[[subsetName]]) <- internal
  x
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
#'   result. Subset names must be unique across the objects, and
#'   differ from the assay names of all of them.
#'
#'   \code{combineCols(x, ...)} and \code{combineRows(x, ...)} combine
#'   objects whose rows (for \code{combineCols}) or columns (for
#'   \code{combineRows}) differ, as for the class they extend: the result
#'   has all of the rows (or columns), matched by name, and fills the
#'   gaps. They keep the subsets of every object too, matching each
#'   subset's rows (or columns) by name to the result's, or by position
#'   with \code{use.names = FALSE}.
#' @details For some classes, \code{combineCols()} and \code{combineRows()}
#'   of the class the objects extend return a different class:
#'   \code{combineCols()} of a \code{TreeSummarizedExperiment} or
#'   \code{SpatialExperiment} returns a \code{SingleCellExperiment}, and
#'   \code{combineRows()} of a \code{SingleCellExperiment},
#'   \code{TreeSummarizedExperiment} or \code{SpatialExperiment} returns a
#'   \code{RangedSummarizedExperiment}. The result is then the matching
#'   \code{ExperimentSubset} class, such as
#'   \code{SubsetSingleCellExperiment}, and each subset is converted to
#'   that class too, so it loses what the class can't hold, as the main
#'   object does: for example its reducedDims and altExps in a
#'   \code{SubsetRangedSummarizedExperiment}.
#'
#'   \code{combineCols()} and \code{combineRows()} use the method for the
#'   class of \code{x} only, so if \code{x} is a plain object, such as a
#'   \code{SingleCellExperiment}, the result is a plain object without
#'   subsets. Convert it with \code{ExperimentSubset()} first.
#' @param ... \code{ExperimentSubset} objects of the same class, which must
#'   have the same rows for \code{cbind} or the same columns for
#'   \code{rbind}. Convert a plain object, such as a
#'   \code{SingleCellExperiment}, with \code{ExperimentSubset()} first,
#'   which keeps all of its data: combining it as it is makes R warn about
#'   "multiple direct matches", and if it comes first, the result is a
#'   plain object without subsets. For \code{combineCols} and
#'   \code{combineRows}, the objects to combine with \code{x}.
#' @param deparse.level Passed on to the method of the class the objects
#'   extend. Default \code{1}.
#' @param x For \code{combineCols} and \code{combineRows}, the first
#'   \code{ExperimentSubset} object to combine.
#' @param delayed,fill Passed on to the method of the class the objects
#'   extend: whether to combine assays as \code{DelayedArray} objects, and
#'   the value for the gaps. Defaults \code{TRUE} and \code{NA}.
#' @param use.names Whether to match rows (for \code{combineCols}) or
#'   columns (for \code{combineRows}) by name, as the method of the class
#'   the objects extend does; if \code{FALSE}, they are matched by position
#'   and the objects must have the same number of them. Default
#'   \code{TRUE}.
#' @return For \code{cbind} and \code{rbind}, an object of the same class
#'   as the objects in \code{...}. For \code{combineCols} and
#'   \code{combineRows}, the \code{ExperimentSubset} class matching the
#'   result of the method of the class \code{x} extends; see Details.
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
#'
#' # a third object, with genes 6 to 10 and 5 genes the others don't have
#' counts3 <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts3) <- paste0("gene", 6:15)
#' colnames(counts3) <- paste0("cellC", 1:20)
#' es3 <- ExperimentSubset(list(counts = counts3))
#' # a subset of the genes only the third object has
#' es3 <- createSubset(es3, subsetName = "newGenes", rows = 6:10)
#' # combine the cells, keeping the genes of both; the gaps are NA
#' esAll <- combineCols(es1, es3)
#' dim(esAll)
#' # the subsets of both objects are kept
#' subsetNames(esAll)
#' # the new subset keeps its genes
#' subsetRownames(esAll, subsetName = "newGenes")
NULL

# Whether each of objects is an ExperimentSubset object, with subsets
.hasSubsets <- function(objects) {
  vapply(objects, methods::is, logical(1), "ExperimentSubset")
}

# The supported class that x is or extends, found from its class
# hierarchy, so a subclass gets the class it extends, not its own (#89);
# NA if x extends none. For an ExperimentSubset object, it is the class
# that the stored objects of its subsets have
.baseClass <- function(x) {
  supported <- c("SummarizedExperiment", "RangedSummarizedExperiment",
                 "SingleCellExperiment", "TreeSummarizedExperiment",
                 "SpatialExperiment")
  # closest class first
  classes <- methods::extends(class(x))
  classes[classes %in% supported][1]
}

# Checks that subset names are unique across the objects to combine, and
# differ from the assay names of the other objects (#67, #73)
.checkSubsetNames <- function(objects) {
  hasSubsets <- .hasSubsets(objects)
  allNames <- unlist(lapply(objects[hasSubsets], subsetNames))
  repeated <- unique(allNames[duplicated(allNames)])
  if (length(repeated) > 0)
    stop("Subset names must be unique across the objects to combine; ",
         paste0("'", repeated, "'", collapse = ", "),
         if (length(repeated) == 1) " is" else " are", " in more than ",
         "one. Create the subsets under different names before combining.")
  # combineCols() and combineRows() keep the assays of every object, so a
  # subset named as another object's assay would be hidden by it (#73)
  for (k in which(hasSubsets)) {
    otherAssays <- unlist(lapply(objects[-k],
                                 SummarizedExperiment::assayNames))
    clashes <- intersect(subsetNames(objects[[k]]), otherAssays)
    if (length(clashes) > 0)
      stop(paste0("'", clashes, "'", collapse = ", "),
           if (length(clashes) == 1) " is a subset" else " are subsets",
           " in one object and an assay in another; subset names must ",
           "differ from assay names. Rename one of them before combining.")
  }
  invisible(NULL)
}

# Checks objects before cbind() (along = "cols") or rbind() ("rows")
# combines them, so nothing is copied when they can't be: subset names
# must be unique across the objects, and for cbind() each object with
# subsets must have the first object's rows, in the same order, since its
# subsets keep their row indices. rbind() already requires the same
# column names (#67)
.checkCombine <- function(objects, along) {
  .checkSubsetNames(objects)
  if (along == "cols") {
    for (k in which(.hasSubsets(objects))) {
      if (!identical(rownames(objects[[k]]), rownames(objects[[1]])))
        stop("cbind() needs the same rows, in the same order, in each ",
             "object with subsets; object ", k, " has different rows ",
             "from object 1.")
    }
  }
  invisible(NULL)
}

# Where the rows and columns of each object with subsets are in out, the
# result of combining objects along "cols" (cbind(), combineCols()) or
# "rows" (rbind(), combineRows()): list(rows =, cols =), each a list with
# a vector of positions in out per object, or NULL where they don't
# change. Along the combined dimension, each object follows the objects
# before it; on the other, positions are matched by name when useNames is
# TRUE, as combineCols() and combineRows() do, and don't change otherwise
.combineMaps <- function(objects, out, along, useNames) {
  byCols <- along == "cols"
  getSize <- if (byCols) ncol else nrow
  getOtherNames <- if (byCols) rownames else colnames
  sizes <- vapply(objects, getSize, integer(1))
  offsets <- cumsum(c(0L, sizes[-length(sizes)]))
  shifted <- matched <- vector("list", length(objects))
  for (k in which(.hasSubsets(objects))) {
    shifted[[k]] <- seq_len(sizes[k]) + offsets[k]
    if (isTRUE(useNames))
      matched[[k]] <- match(getOtherNames(objects[[k]]),
                            getOtherNames(out))
  }
  if (byCols)
    list(rows = matched, cols = shifted)
  else
    list(rows = shifted, cols = matched)
}

# The names a subset stores on one dimension, after its object's names
# there, from, became to in the combined object: combineCols() and
# combineRows() with use.names = FALSE match by position and keep the
# first object's names. Stored names that are still the object's follow;
# names set with subsetRownames<- or subsetColnames<- stay (#73)
.followNames <- function(stored, from, to) {
  if (is.null(from))
    return(if (is.null(stored)) to else stored)
  if (is.null(stored))
    return(stored)
  if (is.null(to))
    return(if (isTRUE(all(stored == from))) NULL else stored)
  same <- which(stored == from)
  stored[same] <- to[same]
  stored
}

# Gives out, the result of combining objects along "cols" or "rows", the
# subsets of every object. cbind() and rbind() of the parent class keep
# only the first object's subsets (#67); combineCols() and combineRows()
# keep none, and return a new object of the class they extend or one it
# extends, which is wrapped as the matching Subset class (#73). Subsets
# made from a main-object assay move to where their rows and columns are
# in out (see .combineMaps()); subsets of subsets keep their indices,
# which are relative to their parent. Each subset's stored object is
# converted to the class out extends, if that differs, and its names on
# the matched dimension follow out's (see .followNames())
.combineSubsets <- function(out, objects, along, useNames = FALSE) {
  if (!methods::is(out, "ExperimentSubset"))
    out <- ExperimentSubset(out)
  maps <- .combineMaps(objects, out, along, useNames)
  baseClass <- .baseClass(out)
  byCols <- along == "cols"
  getOtherIndices <- if (byCols) .rowIndices else .colIndices
  getOtherNames <- if (byCols) rownames else colnames
  subsets <- list()
  for (k in which(.hasSubsets(objects))) {
    # explicit parent names, so a bare one can't match another object's
    object <- .qualifyParents(objects[[k]])
    # matched by position, out has the first object's names (#73)
    renamed <- !isTRUE(useNames) &&
      !identical(getOtherNames(object), getOtherNames(out))
    for (name in subsetNames(object)) {
      subset <- .subsets(object)[[name]]
      rows <- .rowIndices(subset)
      cols <- .colIndices(subset)
      if (is.null(.parentSubsetName(object, .parentAssay(subset)))) {
        if (!is.null(maps$rows[[k]]))
          rows <- maps$rows[[k]][rows]
        if (!is.null(maps$cols[[k]]))
          cols <- maps$cols[[k]][cols]
      }
      internal <- .internalAssay(subset)
      if (!identical(as.character(class(internal)), baseClass))
        internal <- methods::as(internal, baseClass)
      if (renamed) {
        positions <- .rootIndices(object, name, getOtherIndices)
        stored <- getOtherNames(internal)
        followed <- .followNames(stored,
                                 getOtherNames(object)[positions],
                                 getOtherNames(out)[positions])
        if (!identical(followed, stored)) {
          if (byCols)
            rownames(internal) <- followed
          else
            colnames(internal) <- followed
        }
      }
      subsets[[name]] <- AssaySubset(
        subsetName = name,
        rowIndices = rows,
        colIndices = cols,
        parentAssay = .parentAssay(subset),
        internalAssay = internal
      )
    }
  }
  .subsets(out) <- S4Vectors::SimpleList(subsets)
  # SpatialExperiment's cbind() renames sample_ids repeated across the
  # objects; each subset's stored copy follows the main object's
  if (methods::is(out, "SpatialExperiment")) {
    subsets <- .subsets(out)
    for (name in names(subsets)) {
      rootCols <- .rootIndices(out, name, .colIndices)
      internal <- .internalAssay(subsets[[name]])
      internal$sample_id <- out$sample_id[rootCols]
      .internalAssay(subsets[[name]]) <- internal
    }
    .subsets(out) <- subsets
  }
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

#' @title Reduced dimensions of subsets in ExperimentSubset objects
#' @description Get or set the reduced dimensions, such as a PCA or UMAP,
#'   of a subset by naming it with \code{subsetName}. Without
#'   \code{subsetName}, these methods act on the main object, as for
#'   \code{SingleCellExperiment}. They exist for the classes that extend
#'   \code{SingleCellExperiment}: \code{SubsetSingleCellExperiment},
#'   \code{SubsetSpatialExperiment} and \code{SubsetTreeSummarizedExperiment}.
#'   \code{\link{reducedDimNames}} and \code{\link{reducedDimNames<-}}
#'   have their own pages.
#' @param x An \code{ExperimentSubset} object.
#' @param type \code{character(1)} or \code{integer(1)} The name or
#'   position of the reduced dimension. Leave it out for the first one.
#' @param withDimnames \code{logical(1)} Passed on to
#'   \code{\link[SingleCellExperiment]{reducedDim<-}} or
#'   \code{\link[SingleCellExperiment]{reducedDims<-}}. Default \code{TRUE}.
#' @param ... \describe{
#'     \item{\code{subsetName}}{\code{character(1)} Name of the subset whose
#'       reduced dimensions to get or set. Leave it out for the main
#'       object's. Give it by name: for \code{reducedDims<-}, a subset name
#'       in second place is \code{withDimnames}, and an error.}
#'   }
#'   Other arguments are passed on to the \code{SingleCellExperiment}
#'   method, for the subset or the main object.
#' @param value For \code{reducedDim<-}, a matrix with one row per column
#'   of the subset (or of the main object). For \code{reducedDims<-}, a
#'   named list of such matrices, which replaces all of them.
#' @return \code{reducedDim()} returns a matrix and \code{reducedDims()} a
#'   list of them. \code{reducedDim<-} and \code{reducedDims<-} return
#'   \code{x} with the reduced dimensions stored.
#' @seealso \code{\link[SingleCellExperiment]{reducedDims}}
#' @name ExperimentSubset-reducedDim
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # a 2-dimensional embedding of the subset's 8 cells
#' embedding <- matrix(rnorm(16), ncol = 2)
#' # store it in the subset
#' reducedDim(es, "PCA", subsetName = "first") <- embedding
#' # read it back
#' reducedDim(es, "PCA", subsetName = "first")
#' # all of the subset's reduced dimensions
#' reducedDims(es, subsetName = "first")
#' # replace them all with a list of two
#' reducedDims(es, subsetName = "first") <- list(PCA = embedding,
#'                                               UMAP = embedding * 2)
#' reducedDimNames(es, subsetName = "first")
NULL

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

# Exported so ExperimentSubset::`reducedDims<-` still works; it is now
# SingleCellExperiment's generic, with the methods in R/Subset*.R (#26).
#' @importMethodsFrom SingleCellExperiment reducedDims<-
#' @exportMethod reducedDims<-
`.reducedDims<-` <- function(x, ..., value) {
  arglist <- list(...)
  # the package's old generic took subsetName second, where SCE's has
  # withDimnames, so reducedDims(x, "name") <- value now lands there (#26)
  if (is.character(arglist[["withDimnames"]]))
    stop("'withDimnames' must be TRUE or FALSE, not a subset name. To set ",
         "the reducedDims of a subset, name the subset with subsetName =, ",
         "as in reducedDims(x, subsetName = \"", arglist[["withDimnames"]][1],
         "\") <- value.")
  if (!"subsetName" %in% names(arglist))
    return(callNextMethod(x, ..., value = value))
  subsetName <- arglist[["subsetName"]]
  .isSubset(x, subsetName)
  internal <- .internalAssay(.subsets(x)[[subsetName]])
  args <- c(list(internal), .otherArgs(...), list(value = value))
  .internalAssay(.subsets(x)[[subsetName]]) <- do.call(`reducedDims<-`, args)
  x
}

# The package keeps its own reducedDimNames generics: SingleCellExperiment's
# have no ..., so subsetName can't pass through them (#26). Calls without
# subsetName go to SingleCellExperiment's unchanged.

#' @rdname reducedDimNames
setMethod(
  f = "reducedDimNames",
  signature = "ANY",
  definition = function(x, ...) {
    .checkReducedDimNamesArgs(x, ...)
    arglist <- list(...)
    if (!"subsetName" %in% names(arglist))
      return(SingleCellExperiment::reducedDimNames(x))
    subsetName <- arglist[["subsetName"]]
    .isSubset(x, subsetName)
    internal <- .internalAssay(.subsets(x)[[subsetName]])
    SingleCellExperiment::reducedDimNames(internal)
  }
)

# SingleCellExperiment's reducedDimNames() takes only x, so anything in ...
# other than a named subsetName would be dropped without a message (#82).
.checkReducedDimNamesArgs <- function(x, ...) {
  nSubsetNames <- sum(names(list(...)) == "subsetName")
  if (nSubsetNames > 1L)
    stop("reducedDimNames() takes one 'subsetName', not several.")
  others <- .otherArgs(...)
  if (length(others) == 0L)
    return(invisible(NULL))
  otherNames <- names(others)
  if (is.null(otherNames))
    otherNames <- rep("", length(others))
  first <- others[[1]]
  # a lone string by position is a subset name; say how to give it, or,
  # for an x without subsets, give .isSubset()'s error that it has none
  if (nSubsetNames == 0L && length(others) == 1L && otherNames == ""
      && is.character(first) && length(first) == 1L) {
    if (!.hasSubsets(list(x)))
      .isSubset(x, first)
    stop("To get the reducedDimNames of a subset, name the subset with ",
         "subsetName =, as in reducedDimNames(x, subsetName = \"", first,
         "\").")
  }
  otherNames[otherNames == ""] <- "an unnamed argument"
  stop("reducedDimNames() takes only 'x' and 'subsetName'. Remove: ",
       paste(unique(otherNames), collapse = ", "), ".")
}

#' @rdname reducedDimNames-set
setReplaceMethod(
  f = "reducedDimNames",
  signature = "ANY",
  definition = function(x, subsetName, value) {
    if (missing(subsetName)) {
      SingleCellExperiment::reducedDimNames(x) <- value
      return(x)
    }
    .isSubset(x, subsetName)
    internal <- .internalAssay(.subsets(x)[[subsetName]])
    SingleCellExperiment::reducedDimNames(internal) <- value
    .internalAssay(.subsets(x)[[subsetName]]) <- internal
    x
  }
)

#' @title rowData and colData of subsets in ExperimentSubset objects
#' @description Get or set the \code{rowData} or \code{colData} of a subset
#'   by naming it with \code{subsetName}. Without \code{subsetName}, these
#'   methods act on the main object, as for the class \code{x} extends.
#'
#'   A subset stores only its own columns. With \code{parentRowData = TRUE}
#'   or \code{parentColData = TRUE}, the main object's columns come first,
#'   for the subset's rows or columns, and a subset's own column hides the
#'   main object's column of the same name, with a message if their values
#'   differ. When setting a subset's data, columns that came from the main
#'   object, in data read with \code{parentRowData = TRUE} or
#'   \code{parentColData = TRUE}, are left out with a message, so only the
#'   subset's own are stored. To add one column, use
#'   \code{colData(x, subsetName = "name")$column <- values}.
#' @param x An \code{ExperimentSubset} object.
#' @param use.names \code{logical(1)} Whether the \code{rowData} returned
#'   has row names. Default \code{TRUE}.
#' @param ... \describe{
#'     \item{\code{subsetName}}{\code{character(1)} Name of the subset whose
#'       data to get or set. Leave it out for the main object's.}
#'     \item{\code{parentRowData}}{\code{logical(1)} For \code{rowData()}
#'       of a subset, whether to include the main object's columns.
#'       Default \code{FALSE}.}
#'     \item{\code{parentColData}}{\code{logical(1)} For \code{colData()}
#'       of a subset, whether to include the main object's columns.
#'       Default \code{FALSE}.}
#'   }
#'   Without \code{subsetName}, other arguments are passed on to the method
#'   of the class \code{x} extends.
#' @param value A \code{DataFrame} with one row per row (\code{rowData<-})
#'   or column (\code{colData<-}) of the subset or the main object, or
#'   \code{NULL} to remove all columns (a \code{SpatialExperiment} keeps
#'   its \code{sample_id} column). For a subset, a \code{data.frame}
#'   is also accepted and converted to a \code{DataFrame}; for the main
#'   object, \code{colData<-} with a \code{data.frame} is an error.
#' @return \code{rowData()} and \code{colData()} return a \code{DataFrame}.
#'   The setters return \code{x} with the data set.
#' @seealso \code{\link{subsetRowData}}, \code{\link{subsetColData}},
#'   \code{\link[SummarizedExperiment]{rowData}},
#'   \code{\link[SummarizedExperiment]{colData}}
#' @name ExperimentSubset-rowData-colData
#' @examples
#' # an object with 10 genes and 20 cells
#' set.seed(1)
#' counts <- matrix(rpois(200, lambda = 5), nrow = 10)
#' rownames(counts) <- paste0("gene", 1:10)
#' colnames(counts) <- paste0("cell", 1:20)
#' es <- ExperimentSubset(list(counts = counts))
#' # a sample label for every cell of the main object
#' es$sample <- rep(c("A", "B"), each = 10)
#' # a subset of the first 5 genes and the first 8 cells
#' es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
#' # a cluster label for each of the subset's 8 cells
#' colData(es, subsetName = "first")$cluster <- rep(1:2, times = 4)
#' # the subset's own columns
#' colData(es, subsetName = "first")
#' # with the main object's columns first
#' colData(es, subsetName = "first", parentColData = TRUE)
#' # the subset's rowData, which has no columns yet
#' rowData(es, subsetName = "first")
NULL

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

# Arguments in ... other than subsetName and unlink, to pass to a subset's
# internal object, whose methods may not accept them
.otherArgs <- function(...) {
  args <- list(...)
  if (is.null(names(args)))
    return(args)
  args[!names(args) %in% c("subsetName", "unlink")]
}

# TRUE if value is one non-NA string, which a name used in if() must be
.isOneString <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value)
}

# Stops unless subsetName is one non-NA string, before it reaches an if()
# that would fail with a base R error naming no argument (#83)
.checkSubsetName <- function(subsetName) {
  if (!.isOneString(subsetName))
    stop("'subsetName' must be one subset name, a single character string.")
}

# The same for assay()'s i, which names an assay or a subset (#85)
.checkAssayI <- function(i) {
  if (!.isOneString(i))
    stop("'i' must be one name of an assay, a subset or a \"subset/assay\", ",
         "a single character string.")
}

#checks if is a subset
.isSubset <- function(x, subsetName){
  # only the "ANY" reducedDimNames methods can get such an x (#82); that
  # comes first, as no subsetName would work on it
  if (!.hasSubsets(list(x)))
    stop("'x' is a ", class(x)[1], ", which has no subsets. Convert it ",
         "with ExperimentSubset() before using subsetName =.")
  .checkSubsetName(subsetName)
  if(!subsetName %in% subsetNames(x))
    stop(subsetName, " subset does not exist. 'subsetName' must be the ",
         "name of a subset.")
}

# Name of the subset that holds subsetName: the subset itself, or the
# subset storing an assay of that name. Stops if there is none.
.owningSubset <- function(x, subsetName) {
  .checkSubsetName(subsetName)
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
