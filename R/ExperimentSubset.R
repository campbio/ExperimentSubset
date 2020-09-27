setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("MissingOrNullOrCharacter", c("missing", "NULL", "character"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))
setClassUnion("NullOrMissingOrNumericOrCharacter",
              c("NULL", "missing", "numeric", "character"))

#' An S4 class to create subset objects to store inside an \code{ExperimentSubset} object.
#'
#' @slot subsetName Name of the subset.
#' @slot rowIndices Indices of the rows to include in the subset.
#' @slot colIndices Indices of the columns to include in the subset.
#' @slot parentAssay Name of the parent of this subset.
#' @slot internalAssay An internal \code{SingleCellExperiment} object to store additional subset data.
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass(
  "SingleCellSubset",
  slots = representation(
    subsetName = "character",
    rowIndices = "NullOrNumeric",
    colIndices = "NullOrNumeric",
    parentAssay = "NullOrCharacter",
    internalAssay = "SingleCellExperiment"
  )
)

#' @title SingleCellSubset Constructor
#' @description Constructor for creating a \code{SingleCellSubset} object internally by the \code{ExperimentSubset} object.
#' @param subsetName Name of the subset.
#' @param rowIndices Indices of the rows to include in the subset.
#' @param colIndices Indices of the columns to include in the subset.
#' @param parentAssay Name of the parent of this subset.
#' @param internalAssay An internal \code{SingleCellExperiment} object to store additional subset data.
#' @return A \code{SingleCellSubset} object.
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
SingleCellSubset <- function(subsetName = "subset",
                             rowIndices = NULL,
                             colIndices = NULL,
                             parentAssay = "counts",
                             internalAssay = SingleCellExperiment::SingleCellExperiment())
{
  if (grepl("\\s+", subsetName)) {
    subsetName <- gsub("\\s", "", subsetName)
    warning("Removing spaces from subsetName argument.")
  }
  .SingleCellSubset(
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
.ExperimentSubset <- setClass(
  "ExperimentSubset",
  slots = representation(
    root = "SummarizedExperiment",
    subsets = "list"),
  prototype = list(
    subsets = list())
)

#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the \code{ExperimentSubset} object, either through manually specifying the \code{assays}, \code{rowData}, \code{colData} or directly by passing either a \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or objects inherited by these classes. A subset can also be directly created by pasing a named \code{list} to the \code{subset} parameter. This named \code{list} should have parameter values named as \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @param object A \code{SingleCellExperiment} or \code{SummarizedExperiment} object if direct conversion is required.
#' @param ... Additional paramters passed to \code{SingleCellExperiment} constructor.
#' @param subset A named \code{list} if a subset should be created from within the constructor. Named parameters in this list should be \code{subsetName}, \code{rows}, \code{cols} and \code{parentAssay}.
#' @return A \code{ExperimentSubset} object.
#' @export
#' @import BiocStyle
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
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
    if (inherits(object, "SummarizedExperiment")) {
      object <- as(object, "SingleCellExperiment")
      es <- .ExperimentSubset(root = object)
    }
  }
  else{
    se <- SingleCellExperiment::SingleCellExperiment(...)
    es <- .ExperimentSubset(root = se)
  }
  if (!anyNA(subset)) {
    es <- ExperimentSubset::createSubset(
      es,
      subsetName = subset$subsetName,
      rows = subset$rows,
      cols = subset$cols,
      parentAssay = subset$parentAssay
    )
  }
  es
}

#' @title createSubset
#' @description Create a subset from an already available \code{assay} in the input \code{ExperimentSubset} object by specifying the rows and columns to include in the subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Specify the name of the subset to create.
#' @param rows Specify the rows to include in this subset. If \code{missing} or \code{NULL}, all rows are included in the subset. Values can be \code{numeric} or \code{character}. Default \code{NULL}.
#' @param cols Specify the columns to include in this subset. If \code{missing} or \code{NULL}, all columns are included in the subset. Values can be \code{numeric} or \code{character}. Default \code{NULL}.
#' @param parentAssay Specify the parent \code{assay} of the subset. This parent \code{assay} must already be available in the \code{ExperimentSubset} object. If \code{NULL}, the first available main \code{assay} will be marked as parent. Default \code{NULL}.
#' @return An \code{ExperimentSubset} object that now contains the newly created subset.
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
  {
    standardGeneric("createSubset")
  }
)

#' @rdname createSubset
setMethod(
  f = "createSubset",
  signature = c(
    "ExperimentSubset",
    "character",
    "NullOrMissingOrNumericOrCharacter",
    "NullOrMissingOrNumericOrCharacter",
    "MissingOrNullOrCharacter"
  ),
  definition = function(object,
                        subsetName,
                        rows,
                        cols,
                        parentAssay)
  {
    tempAssay <- ""
    if (is.null(parentAssay)) {
      tempAssay <- SummarizedExperiment::assayNames(object@root)[1]
      parentAssay <- tempAssay
    }
    else{
      if (parentAssay %in% SummarizedExperiment::assayNames(object@root)
          || parentAssay %in% subsetAssayNames(object@root)) {
        tempAssay <- parentAssay
      }
      else{
        stop("Input parentAssay does not exist.")
      }
    }
    if (is.character(rows)) {
      rows <-
        match(rows, rownames(
          SummarizedExperiment::assay(object@root, withDimnames = FALSE, tempAssay)
        ))
    }
    if (is.character(cols)) {
      cols <-
        match(cols, colnames(
          SummarizedExperiment::assay(object@root, withDimnames = FALSE, tempAssay)
        ))
    }
    if (is.null(rows)) {
      rows <-
        seq(1, dim(
          SummarizedExperiment::assay(object@root, withDimnames = FALSE, tempAssay)
        )[1])
    }
    if (is.null(cols)) {
      cols <-
        seq(1, dim(
          SummarizedExperiment::assay(object@root, withDimnames = FALSE, tempAssay)
        )[2])
    }
    scs <- SingleCellSubset(
      subsetName = subsetName,
      rowIndices = rows,
      colIndices = cols,
      parentAssay = parentAssay,
      internalAssay = SingleCellExperiment::SingleCellExperiment(list(
        counts = Matrix::Matrix(
          nrow = length(rows),
          ncol = length(cols),
          data = 0,
          sparse = TRUE
        )
      ))
    )

    #Check if NAs introduced in the subset
    tryCatch({
      stats::na.fail(scs@rowIndices)
      stats::na.fail(scs@colIndices)
    }, error = function(e) {
      stop(
        "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent."
      )
    })

    #Remove counts assay from internal SCE object of the subset to save memory
    SummarizedExperiment::assay(scs@internalAssay, withDimnames = FALSE, "counts") <-
      NULL

    object@subsets[[subsetName]] <- scs
    return(object)
  }
)

#' @title storeSubset
#' @description Store a new subset \code{assay} inside a specified subset in the input \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object
#' @param subsetName Specify the name of the existing subset inside which the new subset \code{assay} should be stored.
#' @param inputMatrix The input subset \code{assay}.
#' @param subsetAssayName Specify the name of the new \code{assay} against the \code{inputMatrix} parameter. If \code{NULL}, a new subset is created internally using the \code{createSubset} function. Default \code{NULL}.
#' @return Updated \code{ExperimentSubset} object with the new \code{assay} stored inside the specified subset.
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
  signature = "ExperimentSubset",
  definition = function(object,
                        subsetName,
                        inputMatrix,
                        subsetAssayName = NULL)
  {
    if (is.null(object@subsets[[subsetName]]) &&
        !is.null(subsetAssayName)) {
      stop(paste(subsetName, "does not exist in the subsets slot of the object."))
    }

    if (!is.null(object@subsets[[subsetName]])) {
      if (!all(dim(object@subsets[[subsetName]]@internalAssay) == dim(inputMatrix))
          && is.null(subsetAssayName)) {
        stop(
          "Dimensions of the inputMatrix not equal to the subset. You need to create a new subset with createSubset() function."
        )
      }
    }

    if (is.null(subsetAssayName)) {
      if (subsetName %in% subsetNames(object)) {
        stop(
          paste(
            subsetName,
            "already exists. Please choose a different subsetName parameter."
          )
        )
      }

      object <- createSubset(object,
                             subsetName,
                             rownames(inputMatrix),
                             colnames(inputMatrix),
                             parentAssay = NULL)

      object@subsets[[subsetName]]@internalAssay <-
        SingleCellExperiment::SingleCellExperiment(list(counts = inputMatrix))

    }
    else{
      SummarizedExperiment::assay(object@subsets[[subsetName]]@internalAssay,
                                  withDimnames = FALSE,
                                  subsetAssayName) <- inputMatrix
      rownames(object@subsets[[subsetName]]@internalAssay) <-
        rownames(inputMatrix)
      colnames(object@subsets[[subsetName]]@internalAssay) <-
        colnames(inputMatrix)
    }

    return(object)
  }
)

#' @title assay
#' @description Method to get an \code{assay} from an \code{ExperimentSubset} object or a subset from an \code{ExperimentSubset} object or any object supported by \code{assay} from \code{SummarizedExperiment}.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{assay} from \code{SummarizedExperiment}.
#' @param i Name of an \code{assay} or name of a subset or name of a subset \code{assay}.
#' @param withDimnames Set whether dimnames should be applied to \code{assay}. Default \code{FALSE}.
#' @param ... Additional parameters.
#' @return The \code{assay} from the input object.
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
#' es
setMethod("assay", c("ExperimentSubset", "character"), function(x, i, withDimnames = FALSE, ...) {
  #look at main assays
  if (i %in% SummarizedExperiment::assayNames(x@root)) {
    out <- SummarizedExperiment::assay(x@root, i, withDimnames = withDimnames, ... = ...)
  }
  #look at subsets
  else if (i %in% subsetNames(x)) {
    subsetName <- i
    i <- x@subsets[[subsetName]]@parentAssay
    if (is.null(i)) {
      out <-
        SummarizedExperiment::assay(x@subsets[[subsetName]]@internalAssay, withDimnames = FALSE, "counts")
    }
    else{
      out <- SummarizedExperiment::assay(x@root, withDimnames = FALSE, i)
      out <-
        out[x@subsets[[subsetName]]@rowIndices, x@subsets[[subsetName]]@colIndices]
    }
  }
  #look inside subsets
  else{
    for (j in seq(length(x@subsets))) {
      if (i %in% SummarizedExperiment::assayNames(x@subsets[[j]]@internalAssay)) {
        out <-
          SummarizedExperiment::assay(x@subsets[[j]]@internalAssay, withDimnames = FALSE,  i)
      }
    }
  }
  out
})

#' @title subsetNames
#' @description Retrieves the names of the available subsets in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
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
  signature = "ExperimentSubset",
  definition = function(object)
  {
    return(names(object@subsets))
  }
)
