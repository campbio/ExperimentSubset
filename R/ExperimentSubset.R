setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))


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

#' An S4 class to create an \code{ExperimentSubset} object with support for subsets.
#'
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SubsetSingleCellExperiment <- setClass(
  Class = "SubsetSingleCellExperiment",
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
.SubsetSummarizedExperiment <- setClass(
  Class = "SubsetSummarizedExperiment",
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
.SubsetRangedSummarizedExperiment <- setClass(
  Class = "SubsetRangedSummarizedExperiment",
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
.SubsetVisiumExperiment <- setClass(
  Class = "SubsetVisiumExperiment",
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
.SubsetSpatialExperiment <- setClass(
  Class = "SubsetSpatialExperiment",
  slots = representation(subsets = "list"),
  prototype = list(subsets = list()),
  contains = "SpatialExperiment"
)


#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the \code{ExperimentSubset} object, either through manually specifying the \code{assays}, \code{rowData}, \code{colData} or directly by passing either a \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or objects inherited by these classes. A subset can also be directly created by pasing a named \code{list} to the \code{subset} parameter. This named \code{list} should have parameter values named as \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @param x A \code{SingleCellExperiment} or \code{SummarizedExperiment} object if direct conversion is required.
#' @param ... Additional parameters passed to \code{SingleCellExperiment} constructor.
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
ExperimentSubset <- function(x = NULL,
                             ...,
                             subset = list(
                               subsetName = NA,
                               rows = NA,
                               cols = NA,
                               parentAssay = NA
                             ))
{
  if (!missing(x)
      && !is.null(x)) {
    es <- as(x, paste0("Subset", class(x)))
  }
  else{
    sce <- SingleCellExperiment::SingleCellExperiment(...)
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

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ...) {
    .metadata(x, ...)
  }
)

.metadata <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSummarizedExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetRangedSummarizedExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSingleCellExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetSpatialExperiment", i = "character"),
  definition = function(x, i, ...) {
    .assay(x, i, ...)
  }
)

#' @rdname assay
setMethod(
  f = "assay",
  signature = signature(x = "SubsetVisiumExperiment", i = "character"),
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

#' @title Subset creation method for ExperimentSubset objects
#' @description Create a subset from an already available \code{assay} in the
#'   input \code{ExperimentSubset} object by specifying the rows and columns to
#'   include in the subset.
#' @param x \code{ExperimentSubset} Specify the object from which a subset
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
  def = function(x,
                 subsetName,
                 rows = NULL,
                 cols = NULL,
                 parentAssay = NULL)
    standardGeneric("createSubset"),
  signature = "x"
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "SubsetSummarizedExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
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
    "SubsetSingleCellExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
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
    "SubsetRangedSummarizedExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
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
    "SubsetSpatialExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
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
    "SubsetVisiumExperiment"
  ),
  definition = function(x,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    return(.createSubset(x,
                         subsetName,
                         rows,
                         cols,
                         parentAssay))
  }
)

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
  
  # internalAssay <- SingleCellExperiment(assays = a)
  # internalAssay <- as(internalAssay, gsub("Subset", "", class(x)))
  
  internalAssay <- SingleCellExperiment(assays = a)
  if(inherits(x, "SpatialExperiment"))
    internalAssay <- as(internalAssay, "VisiumExperiment")
  else
    internalAssay <- as(internalAssay, gsub("Subset", "", class(x)))
  internalAssay <- .VisiumSpatialSubset(x, internalAssay, cols)
  
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

.VisiumSpatialSubset <- function(x, internalAssay, cols){
  if(inherits(x, "VisiumExperiment")){
    spatialCoords(internalAssay) <- spatialCoords(x)[cols, ]
    scaleFactors(internalAssay) <- scaleFactors(x)
  }
  else if(inherits(x, "SpatialExperiment")){
    spatialCoords(internalAssay) <- spatialCoords(x)[cols, ]
    internalAssay <- as(internalAssay, "SpatialExperiment")
  }
  return(internalAssay)
}

#' @title Name retrieval method for all subset assays in ExperimentSubset
#'   objects
#' @description Retrieves the names of all the subsets as well as the subset
#'   assays.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
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
  def = function(x)
  {
    standardGeneric("subsetAssayNames")
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

#' @rdname subsetAssayNames
setMethod(
  f = "subsetAssayNames",
  signature = "SubsetVisiumExperiment",
  definition = function(x)
  {
    .subsetAssayNames(x)
  }
)

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
  signature = "SubsetSummarizedExperiment",
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
  signature = "SubsetRangedSummarizedExperiment",
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
  signature = "SubsetSingleCellExperiment",
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
  signature = "SubsetSpatialExperiment",
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
  signature = "SubsetVisiumExperiment",
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
#' @param x \code{ExperimentSubset} Specify the input object.
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
  def = function(x)
  {
    standardGeneric("subsetNames")
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @rdname subsetNames
setMethod(
  f = "subsetNames",
  signature = "SubsetVisiumExperiment",
  definition = function(x)
  {
    return(names(.subsets(x)))
  }
)

#' @title Method for storing new assays in ExperimentSubset objects
#' @description Store a new subset \code{assay} inside a specified subset in the
#'   input \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Specify the input object.
#' @param subsetName \code{character(1)} Specify the name of the existing subset
#'   inside which the new subset \code{assay} should be stored.
#' @param inputMatrix \code{dgCMatrix} The input subset \code{assay}.
#' @param subsetAssayName \code{character(1)} Specify the name of the new
#'   \code{assay} against the \code{inputMatrix} parameter. If \code{NULL}, a
#'   new subset is created internally using the \code{createSubset} function.
#'   Default \code{NULL}.
#' @return Updated \code{ExperimentSubset} object with the new \code{assay}
#'   stored inside the specified subset.
#' @rdname storeSubsetAssay
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
#' es <- storeSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
#' es
setGeneric(
  name = "storeSubsetAssay",
  def = function(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  {
    standardGeneric("storeSubsetAssay")
  }
)

#' @rdname storeSubsetAssay
setMethod(
  f = "storeSubsetAssay",
  signature = c(
    x ="SubsetSummarizedExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .storeSubsetAssay(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubsetAssay
setMethod(
  f = "storeSubsetAssay",
  signature = c(
    x ="SubsetRangedSummarizedExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .storeSubsetAssay(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubsetAssay
setMethod(
  f = "storeSubsetAssay",
  signature = c(
    x ="SubsetSingleCellExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .storeSubsetAssay(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubsetAssay
setMethod(
  f = "storeSubsetAssay",
  signature = c(
    x ="SubsetSpatialExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .storeSubsetAssay(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

#' @rdname storeSubsetAssay
setMethod(
  f = "storeSubsetAssay",
  signature = c(
    x ="SubsetVisiumExperiment",
    subsetName = "character",
    subsetAssayName = "character"),
  definition = function(x,
                        subsetName,
                        inputMatrix,
                        subsetAssayName)
  {
    .storeSubsetAssay(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName)
  }
)

.storeSubsetAssay <- function(x,
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

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetSummarizedExperiment", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetRangedSummarizedExperiment", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetSingleCellExperiment", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetSpatialExperiment", "character"),
                 function(x,
                          i,
                          ...,
                          value) {
                   .assay(x, i, ...) <- value
                   return(x)
                 })

#' @rdname assay<-
setReplaceMethod("assay",
                 c("SubsetVisiumExperiment", "character"),
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
    message("Storing assay inside subset '", i, "'.")
    storeSubsetAssay(
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

#' @title Subset parent hierarchy retrieval method for ExperimentSubset objects
#' @description Retrieves a complete subset to parent link from a specified
#'   subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Specify the name of the subset against
#'   which the subset to parent link should be retrieved.
#' @return A \code{list} containing the parent link of the subset.
#' @rdname subsetParent
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
#' subsetParent(es, "subset1pAssay")
setGeneric(
  name = "subsetParent",
  def = function(x, subsetName)
  {
    standardGeneric("subsetParent")
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetSummarizedExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetSingleCellExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetSpatialExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "SubsetVisiumExperiment",
  definition = function(x, subsetName)
  {
    .subsetParent(x, subsetName)
  }
)

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

#' @title Method for displaying child-parent link structure of subsets in
#'   ExperimentSubset objects
#' @description The function displays the content of an \code{ExperimentSubset}
#'   object including all available main assays, all subsets and the subset
#'   assays inside these subsets. This function also depicts how and in what
#'   order the subsets in the object are linked with their parents. Moreover,
#'   all supplementary data inside the subsets such as \code{reducedDims} and
#'   \code{altExps} are also displayed against each subset entry.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @return Prints all the available subset information against the input
#'   \code{ExperimentSubset} object.
#' @rdname subsetSummary
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' assay(es, "subset1",
#' subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#' subsetSummary(es)
setGeneric(
  name = "subsetSummary",
  def = function(x)
  {
    standardGeneric("subsetSummary")
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

#' @rdname subsetSummary
setMethod(
  f = "subsetSummary",
  signature = "SubsetVisiumExperiment",
  definition = function(x)
  {
    .subsetSummary(x)
  }
)

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

#' @title Get dimensions of subsets in ExperimentSubset objects
#' @description Retrieves the dimensions of the specified subset in an
#'   \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to retrieve the
#'   dimensions from.
#' @return A \code{vector} containing the dimensions of the specified subset
#'   i.e. the number of rows and the number of columns in the subset.
#' @rdname subsetDim
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetDim(es, "subset1")
setGeneric(
  name = "subsetDim",
  def = function(x, subsetName)
  {
    standardGeneric("subsetDim")
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSingleCellExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetSpatialExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("SubsetVisiumExperiment", "character"),
  definition = function(x, subsetName)
  {
    dim(.internalAssay(.subsets(x)[[subsetName]]))
  }
)


#' @title reducedDimNames
#' @description A wrapper to the \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset from which the \code{reducedDimNames} should be fetched from. If \code{missing}, \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method is called on the main object.
#' @return The \code{reducedDimNames} from the specified subset or same as \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} when \code{subsetName} is \code{missing}.
#' @rdname reducedDimNames
#' @export
setGeneric(
  name = "reducedDimNames",
  def = function(x, ...)
  {
    standardGeneric("reducedDimNames")
  }
)

#' @rdname reducedDimNames
setMethod(
  f = "reducedDimNames",
  signature = c("ANY"),
  definition = function(x, ...)
  {
    arglist <- list(...)
    if(!"subsetName" %in% names(arglist))
      return(SingleCellExperiment::reducedDimNames(x))
    subsetName = arglist[["subsetName"]]
    .isSubset(x, subsetName)
    SingleCellExperiment::reducedDimNames(.internalAssay(.subsets(x)[[subsetName]]))
  }
)

#' @title Subset count method for ExperimentSubset objects
#' @description Get the total count of the available subsets in an
#'   \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @return A \code{numeric} value representing the total count of the subsets.
#' @rdname subsetCount
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetCount(es)
setGeneric(
  name = "subsetCount",
  def = function(x)
  {
    standardGeneric("subsetCount")
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "SubsetVisiumExperiment",
  definition = function(x)
  {
    return(length(subsetNames(x)))
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "character"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "missing"),
  definition = function(x, e,  ...) {
    .altExp(x, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "numeric"),
  definition = function(x, e,  ...) {
    .altExp(x, e, ...)
  }
)

.altExp <- function(x, e, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod())
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExp(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ...) {
    .altExps(x, ...)
  }
)

.altExps <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExps(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ...) {
    .altExpNames(x, ...)
  }
)

.altExpNames <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]]))
}

#' @title rownames
#' @description Get \code{rownames} from an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{rownames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{rownames} in \code{BiocGenerics} package.
#' @param ... Additional parameters amd \code{subsetName} parameter to pass the name of the subset to get \code{rownames} from.
#' @param subsetName Name of the subset to get \code{rownames} from. If \code{missing}, \code{rownames} from main object are returned.
#' @return A \code{vector} of \code{rownames}.
#' @rdname rownames
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' rownames(es, subsetName = "subset1")
setGeneric(
  name = "rownames",
  def = function(x, ...)
  {
    standardGeneric("rownames")
  }
)

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

#' @title colnames
#' @description Get \code{colnames} from an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param ... Additional parameters amd \code{subsetName} parameter to pass the name of the subset to get \code{colnames} from.
#' @param subsetName Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @return A \code{vector} of \code{colnames}.
#' @rdname colnames
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' colnames(es, subsetName = "subset1")
setGeneric(
  name = "colnames",
  def = function(x, ...)
  {
    standardGeneric("colnames")
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

#' @rdname altExpNames<-
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSingleCellExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

#' @rdname altExpNames<-
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetSpatialExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

#' @rdname altExpNames<-
setReplaceMethod(
  f = "altExpNames",
  signature = signature(x = "SubsetVisiumExperiment", value = "character"),
  definition = function(x, ..., value) {
    .altExpNames(x, ...) <- value
    return(x)
  }
)

'.altExpNames<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExpNames(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @title reducedDimNames<-
#' @description A wrapper to the \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the \code{reducedDimNames<-} should be set to. If \code{missing}, \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method is called on the main object.
#' @param value Input value same as \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @return Input object with \code{reducedDimNames<-} set.
#' @rdname reducedDimNames-set
#' @export
setGeneric(
  name = "reducedDimNames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("reducedDimNames<-")
  }
)

#' @rdname reducedDimNames-set
setReplaceMethod(
  f = "reducedDimNames",
  signature = "ANY",
  definition = function(x, subsetName, value)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      SingleCellExperiment::reducedDimNames(.internalAssay(.subsets(x)[[subsetName]])) <-
        value
    }
    else{
      SingleCellExperiment::reducedDimNames(x) <- value
    }
    x
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "character"),
  definition = function(x, e,  ..., value) {
    .altExp(x, e, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "character"),
  definition = function(x, e,  ..., value) {
    .altExp(x, e, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "character"),
  definition = function(x, e,  ..., value) {
    .altExp(x, e, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "missing"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "missing"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "missing"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSingleCellExperiment", e = "numeric"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetSpatialExperiment", e = "numeric"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

#' @rdname altExp<-
setReplaceMethod(
  f = "altExp",
  signature = signature(x = "SubsetVisiumExperiment", e = "numeric"),
  definition = function(x, e,  ..., value) {
    .altExp(x, ...) <- value
    return(x)
  }
)

'.altExp<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod())
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExp(.internalAssay(.subsets(x)[[subsetName]]), ...) <- value
  return(x)
}

#' @rdname altExps<-
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ..., value) {
    .altExps(x, ...) <- value
    return(x)
  }
)

#' @rdname altExps<-
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ..., value) {
    .altExps(x, ...) <- value
    return(x)
  }
)

#' @rdname altExps<-
setReplaceMethod(
  f = "altExps",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ..., value) {
    .altExps(x, ...) <- value
    return(x)
  }
)

'.altExps<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  altExps(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)
#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

#' @rdname metadata<-
setReplaceMethod(
  f = "metadata",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ..., value) {
    .metadata(x, ...) <- value
    return(x)
  }
)

'.metadata<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  metadata(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @title Count method for subset assays in ExperimentSubset objects
#' @description Get the count of the total available subsets and the subset
#'   assays inside these subsets in an \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @return A \code{numeric} value representing the sum of the subset count and
#'   subset assay count.
#' @rdname subsetAssayCount
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' assay(es, "subset1",
#' subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#' subsetAssayCount(es)
setGeneric(
  name = "subsetAssayCount",
  def = function(x)
  {
    standardGeneric("subsetAssayCount")
  }
)


#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetRangedSummarizedExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSingleCellExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetSpatialExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "SubsetVisiumExperiment",
  definition = function(x)
  {
    return(length(subsetAssayNames(x)))
  }
)

#' @title Accessor method for rowData from subsets in ExperimentSubset objects
#' @description Get \code{rowData} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{rowData} from.
#' @return The \code{rowData} from input object.
#' @rdname subsetRowData
#' @export
setGeneric(
  name = "subsetRowData",
  def = function(x, subsetName)
  {
    standardGeneric("subsetRowData")
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetSingleCellExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetSpatialExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("SubsetVisiumExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetRowData(x, subsetName)
  }
)


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

#' @title Accessor method for colData from subsets in ExperimentSubset objects
#' @description Get \code{colData} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{colData} from.
#' @return The \code{colData} from input object.
#' @rdname subsetColData
#' @export
setGeneric(
  name = "subsetColData",
  def = function(x, subsetName)
  {
    standardGeneric("subsetColData")
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetRangedSummarizedExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetSingleCellExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetSpatialExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("SubsetVisiumExperiment", "character"),
  definition = function(x, subsetName)
  {
    .subsetColData(x, subsetName)
  }
)

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

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "character"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "missing"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "missing"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)
#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

#' @rdname reducedDim
setMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "numeric"),
  definition = function(x, type, ...) {
    .reducedDim(x, type, ...)
  }
)

.reducedDim <- function(x, type, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  reducedDim(.internalAssay(.subsets(x)[[subsetName]]), type)
}

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "character"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "character"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "character"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "missing"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "missing"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "missing"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSingleCellExperiment", type = "numeric"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetSpatialExperiment", type = "numeric"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

#' @rdname reducedDim<-
setReplaceMethod(
  f = "reducedDim",
  signature = signature(x = "SubsetVisiumExperiment", type = "numeric"),
  definition = function(x, type, ..., value) {
    .reducedDim(x, type, ...) <- value
    return(x)
  }
)

'.reducedDim<-' <- function(x, type, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  reducedDim(.internalAssay(.subsets(x)[[subsetName]]), type, ...) <- value
  return(x)
}

#' @rdname reducedDims
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

#' @rdname reducedDims
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

#' @rdname reducedDims
setMethod(
  f = "reducedDims",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .reducedDims(x, ...)
  }
)

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
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:1500), cols = c(1:1500),
#' parentAssay = "counts")
#' reducedDims(es, subsetName = "subset1") <- list(
#' PCA_1 = scater::calculatePCA(assay(es, "subset1")),
#' PCA_2 = scater::calculatePCA(assay(es, "subset1")))
#' reducedDims(es, subsetName = "subset1")
setGeneric(
  name = "reducedDims<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("reducedDims<-")
  }
)

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

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

#' @rdname rowData
setMethod(
  f = "rowData",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ...) {
    .rowData(x, ...)
  }
)

.rowData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  subsetRowData(x, subsetName)
}

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetSpatialExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

#' @rdname colData
setMethod(
  f = "colData",
  signature = signature(x = "SubsetVisiumExperiment"),
  definition = function(x, ...) {
    .colData(x, ...)
  }
)

.colData <- function(x, ...){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  subsetColData(x, subsetName)
}

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSingleCellExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetSpatialExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

#' @rdname rowData<-
setReplaceMethod(
  f = "rowData",
  signature = signature(x = "SubsetVisiumExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .rowData(x, ...) <- value
    return(x)
  }
)

'.rowData<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  rowData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetRangedSummarizedExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSingleCellExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetSpatialExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

#' @rdname colData<-
setReplaceMethod(
  f = "colData",
  signature = signature(x = "SubsetVisiumExperiment", value = "DataFrame"),
  definition = function(x, ..., value) {
    .colData(x, ...) <- value
    return(x)
  }
)

'.colData<-' <- function(x, ..., value){
  arglist <- list(...)
  if(!"subsetName" %in% names(arglist))
    return(callNextMethod(...))
  subsetName = arglist[["subsetName"]]
  .isSubset(x, subsetName)
  colData(.internalAssay(.subsets(x)[[subsetName]])) <- value
  return(x)
}

#' @title subsetColnames
#' @description Get \code{colnames} from an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @return A \code{vector} of \code{colnames}.
#' @rdname subsetColnames
#' @export
setGeneric(
  name = "subsetColnames",
  def = function(x, subsetName)
  {
    standardGeneric("subsetColnames")
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

#' @rdname subsetColnames
setMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetVisiumExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetColnames(x, subsetName)
  }
)

.subsetColnames <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    colnames(x)[.colIndices(.subsets(x)[[subsetName]])]
  }
}

#' @title subsetColnames
#' @description Set \code{colnames} to an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @param value Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @return A \code{vector} of \code{colnames}.
#' @rdname subsetColnames
#' @export
setGeneric(
  name = "subsetColnames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetColnames<-")
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetColnames
setReplaceMethod(
  f = "subsetColnames",
  signature = c(x = "SubsetVisiumExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetColnames(x, subsetName) <- value
    return(x)
  }
)

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

#' @title subsetRownames
#' @description Get \code{rownames} from an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @return A \code{vector} of \code{colnames}.
#' @rdname subsetRownames
#' @export
setGeneric(
  name = "subsetRownames",
  def = function(x, subsetName)
  {
    standardGeneric("subsetRownames")
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

#' @rdname subsetRownames
setMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetVisiumExperiment", subsetName = "character"),
  definition = function(x, subsetName)
  {
    .subsetRownames(x, subsetName)
  }
)

.subsetRownames <- function(x, subsetName){
  if (subsetName %in% subsetNames(x)) {
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
  }
  else if (subsetName %in% subsetAssayNames(x)) {
    subsetName <- .getParentAssayName(x, subsetName)
    rownames(x)[.rowIndices(.subsets(x)[[subsetName]])]
  }
}

#' @title subsetRownames
#' @description Set \code{colnames} to an \code{ExperimentSubset} object or a subset in the \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{colnames} in \code{BiocGenerics} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @param value Name of the subset to get \code{colnames} from. If \code{missing}, \code{colnames} from main object are returned.
#' @return A \code{vector} of \code{colnames}.
#' @rdname subsetRownames
#' @export
setGeneric(
  name = "subsetRownames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetRownames<-")
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetRangedSummarizedExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSingleCellExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetSpatialExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

#' @rdname subsetRownames
setReplaceMethod(
  f = "subsetRownames",
  signature = c(x = "SubsetVisiumExperiment", subsetName = "character"),
  definition = function(x, subsetName, value)
  {
    .subsetRownames(x, subsetName) <- value
    return(x)
  }
)

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
