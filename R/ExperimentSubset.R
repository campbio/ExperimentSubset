setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))

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

#' An S4 class to create subset objects to store inside an
#' \code{ExperimentSubset} object.
#'
#' @slot subsetName \code{character(1)} Name of the subset.
#' @slot rowIndices \code{vector("numeric")} Indices of the rows to include in
#'   the subset.
#' @slot colIndices \code{vector("numeric")} Indices of the columns to include
#'   in the subset.
#' @slot parentAssay \code{character(1)} Name of the parent of this subset.
#' @slot internalAssay \code{SummarizedExperiment} An internal experiment object
#'   to store additional subset data.
#' @import methods
.AssaySubset <- setClass(
  "AssaySubset",
  slots = representation(
    subsetName = "character",
    rowIndices = "NullOrNumeric",
    colIndices = "NullOrNumeric",
    parentAssay = "NullOrCharacter",
    internalAssay = "ANY"
  )
)

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
    warning("Removing spaces from subsetName argument.")
  }
  .AssaySubset(
    subsetName = subsetName,
    rowIndices = rowIndices,
    colIndices = colIndices,
    parentAssay = parentAssay,
    internalAssay = internalAssay
  )
}

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.ExperimentSubsetSCE <- setClass(
  Class = "ExperimentSubsetSCE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SingleCellExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ExperimentSubsetSE <- setClass(
  Class = "ExperimentSubsetSE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SummarizedExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.ExperimentSubsetRSE <- setClass(
  Class = "ExperimentSubsetRSE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "RangedSummarizedExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SpatialExperiment VisiumExperiment
.ExperimentSubsetVE <- setClass(
  Class = "ExperimentSubsetVE",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "VisiumExperiment"
)

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SpatialExperiment SpatialExperiment
.ExperimentSubsetSP <- setClass(
  Class = "ExperimentSubsetSP",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SpatialExperiment"
)


#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the \code{ExperimentSubset} object, either through manually specifying the \code{assays}, \code{rowData}, \code{colData} or directly by passing either a \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or objects inherited by these classes. A subset can also be directly created by pasing a named \code{list} to the \code{subset} parameter. This named \code{list} should have parameter values named as \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @param object A \code{SingleCellExperiment} or \code{SummarizedExperiment} object if direct conversion is required.
#' @param ... Additional parameters passed to \code{SingleCellExperiment} constructor.
#' @param subset A named \code{list} if a subset should be created from within the constructor. Named parameters in this list should be \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @return A \code{ExperimentSubset} object.
#' @export
#' @import BiocStyle
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
ExperimentSubset <- function(object,
                             ...,
                             subset = list(
                               subsetName = NA,
                               rows = NA,
                               cols = NA,
                               parentAssay = NA
                             ))
{
  if (!missing(object)) {
    # if (inherits(object, "SummarizedExperiment")) {
    #   object <- as(object, "SingleCellExperiment")
    #   es <- .ExperimentSubset(object)
    # }
    if(inherits(object, "VisiumExperiment")){
      es <- .ExperimentSubsetVE(object)
    }
    else if(inherits(object, "SpatialExperiment")){
      es <- .ExperimentSubsetSP(object)
    }
    else if(inherits(object, "SingleCellExperiment")){
      es <- .ExperimentSubsetSCE(object)
    }
    else if(inherits(object, "RangedSummarizedExperiment")){
      es <- .ExperimentSubsetRSE(object)
    }
    else if(inherits(object, "SummarizedExperiment")){
      es <- .ExperimentSubsetSE(object)
    }
  }
  else{
    sce <- SingleCellExperiment::SingleCellExperiment(...)
    es <- .ExperimentSubset(sce)
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

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "ExperimentSubsetSE"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "ExperimentSubsetSCE"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  metadata(x@subsets[[subsetName]]@internalAssay)
}

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "ExperimentSubsetSE", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "ExperimentSubsetRSE", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "ExperimentSubsetSCE", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "ExperimentSubsetSP", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "ExperimentSubsetVE", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

.assay <- function(x, i, ...){
  out <- NULL
  #look at main assays
  if (i %in% assayNames(x)) {
    out <- callNextMethod()
  }
  #look at subsets
  else if (i %in% subsetNames(x)) {
    subsetName <- i
    i <- x@subsets[[subsetName]]@parentAssay #change to .subsets
    if (is.null(i)) {
      out <- assay(x = x@subsets[[subsetName]]@internalAssay, i = "counts", ...)
    }
    else{
      out <- assay(x = x, i = i, ...)
      out <- out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
    }
  }
  #look inside subsets
  else{
    for (j in seq(length(x@subsets))) {
      if (i %in% assayNames(x@subsets[[j]]@internalAssay)) {
        out <- assay(x = x@subsets[[j]]@internalAssay, i = i)
      }
    }
  }
  out
  # arglist <- list(...)
  # if(!"subsetName" %in% names(arglist))
  #   return(callNextMethod(...))
  # subsetName = arglist[["subsetName"]]
  # assay(x@subsets[[subsetName]]@internalAssay, i)
}

#' @title Subset creation method for ExperimentSubset objects
#' @description Create a subset from an already available \code{assay} in the
#'   input \code{ExperimentSubset} object by specifying the rows and columns to
#'   include in the subset.
#' @param object \code{ExperimentSubset} Specify the object from which a subset
#'   should be created. Input can also be any object inherited from
#'   \code{SummarizedExperiment} for immediate conversion and subset formation.
#' @param subsetName \code{character(1)} Specify the name of the subset to
#'   create.
#' @param rows \code{vector("numeric")} Specify the rows to include in this
#'   subset. If \code{missing} or \code{NULL}, all rows are included in the
#'   subset. Values can be \code{numeric} or \code{character}. Default
#'   \code{NULL}.
#' @param cols \code{vector("numeric")} Specify the columns to include in this
#'   subset. If \code{missing} or \code{NULL}, all columns are included in the
#'   subset. Values can be \code{numeric} or \code{character}. Default
#'   \code{NULL}.
#' @param parentAssay \code{character(1)} Specify the parent \code{assay} of the
#'   subset. This parent \code{assay} must already be available in the
#'   \code{ExperimentSubset} object. If \code{NULL}, the first available main
#'   \code{assay} will be marked as parent. Default \code{NULL}.
#' @return An \code{ExperimentSubset} object that now contains the newly created
#'   subset.
#' @rdname createSubset
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' es
setGeneric(
  name = "createSubset",
  def = function(object,
                 subsetName,
                 rows = NULL,
                 cols = NULL,
                 parentAssay = NULL)
    standardGeneric("createSubset"),
  signature = "object"
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubsetSE"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(object,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubsetSCE"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(object,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubsetRSE"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(object,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubsetSP"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(object,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubsetVE"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(object,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

.createSubset <- function(object,
                          subsetName,
                          rows,
                          cols,
                          parentAssay){
  #checking parameters
  stopifnot(
    is.character(subsetName),
    is.null(rows) || is.numeric(rows) || is.character(rows),
    is.null(cols) || is.numeric(cols) || is.character(cols),
    is.null(parentAssay) || is.character(parentAssay)
  )
  
  tempAssay <- ""
  if (is.null(parentAssay)) {
    tempAssay <- assayNames(object)[1]
    parentAssay <- tempAssay
  }
  else{
    test <- parentAssay %in% assayNames(object) || 
      parentAssay %in% subsetAssayNames(object)
    if (test) {
      tempAssay <- parentAssay
    }
    else{
      stop("Input parentAssay does not exist.")
    }
  }
  if (is.character(rows)) {
    rows <-
      match(rows, rownames(
        assay(object, withDimnames = TRUE, tempAssay)
      ))
  }
  if (is.character(cols)) {
    cols <-
      match(cols, colnames(
        assay(object, withDimnames = TRUE, tempAssay)
      ))
  }
  if (is.null(rows)) {
    rows <-
      seq(1, dim(
        assay(object, withDimnames = FALSE, tempAssay)
      )[1])
  }
  if (is.null(cols)) {
    cols <-
      seq(1, dim(
        assay(object, withDimnames = FALSE, tempAssay)
      )[2])
  }
  
  #Check if count of stored row/column indices greater than the subset
  test <- length(rows) > dim(assay(object, withDimnames = FALSE, tempAssay))[1] || 
    length(cols) > dim(assay(object, withDimnames = FALSE, tempAssay))[2]
  if (test) {
    stop("More rows or columns selected than available in the parentAssay.")
  }
  
  # Create an initial object for the internalAssay
  a <- list(Matrix(
    nrow = length(rows),
    ncol = length(cols),
    data = 0,
    sparse = TRUE))
  names(a) <- "temp"
  internalAssay <- SummarizedExperiment(assays = a)
  
  # Convert to class of root object (e.g. SingleCellExperiment)
  if(inherits(object, "VisiumExperiment")){
    internalAssay <- as(internalAssay, "VisiumExperiment")
  }
  else if(inherits(object, "SpatialExperiment")){
    internalAssay <- as(internalAssay, "SpatialExperiment")
  }
  else if(inherits(object, "SingleCellExperiment")){
    internalAssay <- as(internalAssay, "SingleCellExperiment")
  }
  else if(inherits(object, "RangedSummarizedExperiment")){
    internalAssay <- as(internalAssay, "RangedSummarizedExperiment")
  }
  else if(inherits(object, "SummarizedExperiment")){
    internalAssay <- as(internalAssay, "SummarizedExperiment")
  }
  #add more conditions here but may need a better solution
  
  scs <- AssaySubset(
    subsetName = subsetName,
    rowIndices = rows,
    colIndices = cols,
    parentAssay = parentAssay,
    internalAssay = internalAssay
  )
  assay(.internalAssay(scs),
        withDimnames = FALSE, "temp") <- NULL
  
  #Check if NAs introduced in the subset
  tryCatch({
    na.fail(.rowIndices(scs))
    na.fail(.colIndices(scs))
  }, error = function(e) {
    stop(
      "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent."
    )
  })
  
  .subsets(object)[[subsetName]] <- scs
  return(object)
}

#' @title Name retrieval method for all subset assays in ExperimentSubset
#'   objects
#' @description Retrieves the names of all the subsets as well as the subset
#'   assays.
#' @param object \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @return A \code{vector} containing the names of the subsets and the subset
#'   assays available in the input \code{ExperimentSubset} object.
#' @rdname subsetAssayNames
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' assay(es, "subset1",
#' subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#' subsetAssayNames(es)
setGeneric(
  name = "subsetAssayNames",
  def = function(object)
  {
    standardGeneric("subsetAssayNames")
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "ExperimentSubsetSE",
  definition = function(object)
  {
    .subsetAssayNames(object)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "ExperimentSubsetRSE",
  definition = function(object)
  {
    .subsetAssayNames(object)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "ExperimentSubsetSCE",
  definition = function(object)
  {
    .subsetAssayNames(object)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "ExperimentSubsetSP",
  definition = function(object)
  {
    .subsetAssayNames(object)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "ExperimentSubsetVE",
  definition = function(object)
  {
    .subsetAssayNames(object)
  }
)

.subsetAssayNames <- function(object){
  tempNames <- names(.subsets(object))
  if (length(.subsets(object)) > 0) {
    for (i in seq(length(.subsets(object)))) {
      tempNames <-
        c(
          tempNames,
          assayNames(.internalAssay(.subsets(object)[[i]]))
        )
    }
  }
  return(tempNames)
}

#' @title show
#' @description Show the \code{ExperimentSubset} object
#' @param object Input \code{ExperimentSubset} object.
#' @return Displays the overall contents of the \code{ExperimentSubset} object.
#' @rdname show
#' @export
#' @importMethodsFrom SingleCellExperiment show
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
setMethod(
  f = "show",
  signature = "ExperimentSubsetSE",
  definition = function(object)
  {
    .show(object)
  }
)

#' @title show
#' @description Show the \code{ExperimentSubset} object
#' @param object Input \code{ExperimentSubset} object.
#' @return Displays the overall contents of the \code{ExperimentSubset} object.
#' @rdname show
#' @export
#' @importMethodsFrom SingleCellExperiment show
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
setMethod(
  f = "show",
  signature = "ExperimentSubsetRSE",
  definition = function(object)
  {
    .show(object)
  }
)

#' @title show
#' @description Show the \code{ExperimentSubset} object
#' @param object Input \code{ExperimentSubset} object.
#' @return Displays the overall contents of the \code{ExperimentSubset} object.
#' @rdname show
#' @export
#' @importMethodsFrom SingleCellExperiment show
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
setMethod(
  f = "show",
  signature = "ExperimentSubsetSCE",
  definition = function(object)
  {
    .show(object)
  }
)

#' @title show
#' @description Show the \code{ExperimentSubset} object
#' @param object Input \code{ExperimentSubset} object.
#' @return Displays the overall contents of the \code{ExperimentSubset} object.
#' @rdname show
#' @export
#' @importMethodsFrom SingleCellExperiment show
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
setMethod(
  f = "show",
  signature = "ExperimentSubsetSP",
  definition = function(object)
  {
    .show(object)
  }
)

#' @title show
#' @description Show the \code{ExperimentSubset} object
#' @param object Input \code{ExperimentSubset} object.
#' @return Displays the overall contents of the \code{ExperimentSubset} object.
#' @rdname show
#' @export
#' @importMethodsFrom SingleCellExperiment show
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
setMethod(
  f = "show",
  signature = "ExperimentSubsetVE",
  definition = function(object)
  {
    .show(object)
  }
)

.show <- function(object){
  callNextMethod()
  cat("subsets(", length(subsetNames(object)), "): ",
      sep = "")
  cat(subsetNames(object))
  cat("\nsubsetAssays(", length(subsetAssayNames(object)), "): ",
      sep = "")
  cat(subsetAssayNames(object))
}


#' @title Get names of subsets in ExperimentSubset objects
#' @description Retrieves the names of the available subsets in an
#'   \code{ExperimentSubset} object.
#' @param object \code{ExperimentSubset} Specify the input object.
#' @return A \code{vector} of subset names.
#' @rdname subsetNames
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetNames(es)
setGeneric(
  name = "subsetNames",
  def = function(object)
  {
    standardGeneric("subsetNames")
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "ExperimentSubsetSE",
  definition = function(object)
  {
    return(names(.subsets(object)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "ExperimentSubsetRSE",
  definition = function(object)
  {
    return(names(.subsets(object)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "ExperimentSubsetSCE",
  definition = function(object)
  {
    return(names(.subsets(object)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "ExperimentSubsetSP",
  definition = function(object)
  {
    return(names(.subsets(object)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "ExperimentSubsetVE",
  definition = function(object)
  {
    return(names(.subsets(object)))
  }
)

#' @title Method for storing new assays in ExperimentSubset objects
#' @description Store a new subset \code{assay} inside a specified subset in the
#'   input \code{ExperimentSubset} object.
#' @param object \code{ExperimentSubset} Specify the input object.
#' @param subsetName \code{character(1)} Specify the name of the existing subset
#'   inside which the new subset \code{assay} should be stored.
#' @param inputMatrix \code{dgCMatrix} The input subset \code{assay}.
#' @param subsetAssayName \code{character(1)} Specify the name of the new
#'   \code{assay} against the \code{inputMatrix} parameter. If \code{NULL}, a
#'   new subset is created internally using the \code{createSubset} function.
#'   Default \code{NULL}.
#' @return Updated \code{ExperimentSubset} object with the new \code{assay}
#'   stored inside the specified subset.
#' @rdname storeSubset
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' counts1p <- assay(es, "subset1")
#' counts1p[,] <- counts1p[,] + 1
#' es <- storeSubset(es, "subset1", counts1p, "scaledSubset1")
#' es
setGeneric(
  name = "storeSubset",
  def = function(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  {
    standardGeneric("storeSubset")
  }
)

#' @rdname storeSubset
setMethod(
  f = "storeSubset",
  signature = "ExperimentSubsetSE",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    .storeSubset(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubset
setMethod(
  f = "storeSubset",
  signature = "ExperimentSubsetRSE",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    .storeSubset(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubset
setMethod(
  f = "storeSubset",
  signature = "ExperimentSubsetSCE",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    .storeSubset(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubset
setMethod(
  f = "storeSubset",
  signature = "ExperimentSubsetSP",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    .storeSubset(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubset
setMethod(
  f = "storeSubset",
  signature = "ExperimentSubsetVE",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    .storeSubset(object,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

.storeSubset <- function(object,
                         subsetName,
                         inputMatrix,
                         subsetAssayName = NULL){
  test <- is.null(.subsets(object)[[subsetName]]) &&
    !is.null(subsetAssayName)
  if (test) {
    stop(subsetName, " does not exist in the subsets slot of the object.")
  }
  
  if (!is.null(.subsets(object)[[subsetName]])) {
    test <- !all(dim(.internalAssay(.subsets(object)[[subsetName]])) == dim(inputMatrix)) && 
      is.null(subsetAssayName)
    if (test) {
      stop(
        "Dimensions of the inputMatrix not equal to the subset. You need to create a new subset with createSubset() function."
      )
    }
  }
  
  if (is.null(subsetAssayName)) {
    if (subsetName %in% subsetNames(object)) {
      stop(subsetName,
           " already exists. Please choose a different subsetName parameter."
      )
    }
    
    object <- createSubset(
      object,
      subsetName,
      rownames(inputMatrix),
      colnames(inputMatrix),
      parentAssay = NULL
    )
    
    internalAssay <- SummarizedExperiment(list(counts = inputMatrix))
    if(inherits(object, "VisiumExperiment")){
      internalAssay <- as(internalAssay, "VisiumExperiment")
    }
    else if(inherits(object, "SpatialExperiment")){
      internalAssay <- as(internalAssay, "SpatialExperiment")
    }
    else if(inherits(object, "SingleCellExperiment")){
      internalAssay <- as(internalAssay, "SingleCellExperiment")
    }
    else if(inherits(object, "RangedSummarizedExperiment")){
      internalAssay <- as(internalAssay, "RangedSummarizedExperiment")
    }
    else if(inherits(object, "SummarizedExperiment")){
      internalAssay <- as(internalAssay, "SummarizedExperiment")
    }
    .internalAssay(.subsets(object)[[subsetName]]) <- internalAssay
      
    
  }
  else{
    assay(.internalAssay(.subsets(object)[[subsetName]]),
                                withDimnames = FALSE,
                                subsetAssayName) <- inputMatrix
    rownames(.internalAssay(.subsets(object)[[subsetName]])) <-
      rownames(inputMatrix)
    colnames(.internalAssay(.subsets(object)[[subsetName]])) <-
      colnames(inputMatrix)
  }
  
  return(object)
}

#' @rdname assay<-
setReplaceMethod("assay",
                 c("ExperimentSubsetSE", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("ExperimentSubsetRSE", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("ExperimentSubsetSCE", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("ExperimentSubsetSP", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("ExperimentSubsetVE", "character"),
                 function(x,
                          i,
                          ...,
                          subsetAssayName = NULL,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

'.assay<-' <- function(x, i, ..., subsetAssayName = NULL, value){
  if ((nrow(value) != nrow(x))
      || (ncol(value) != ncol(x))) {
    storeSubset(
      object = x,
      subsetName = i,
      inputMatrix = value,
      subsetAssayName = subsetAssayName
    )
  }
  else{
    callNextMethod(...)
  }
}