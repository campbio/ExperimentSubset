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
#' @slot internalAssay An internal experiment object to store additional subset data.
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.SingleCellSubset <- setClass(
  "SingleCellSubset",
  slots = representation(
    subsetName = "character",
    rowIndices = "NullOrNumeric",
    colIndices = "NullOrNumeric",
    parentAssay = "NullOrCharacter",
    internalAssay = "ANY"
  )
)

#' @title SingleCellSubset Constructor
#' @description Constructor for creating a experiment object internally by the
#'   \code{ExperimentSubset} object.
#' @param subsetName Name of the subset.
#' @param rowIndices Indices of the rows to include in the subset.
#' @param colIndices Indices of the columns to include in the subset.
#' @param parentAssay Name of the parent of this subset.
#' @param internalAssay An internal object to store additional subset data.
#' @return A \code{SingleCellSubset} object.
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

#' An S4 class to create an \code{ExperimentSubset} object with support for
#' subsets.
#'
#' @slot root The root object from which all consequent subsets will be created.
#'   This can be any object that is inherited from \code{SummarizedExperiment}.
#' @slot subsets A \code{list} of \code{SingleCellSubset} objects.
#' @export
#' @import methods
.ExperimentSubset <- setClass(
  "ExperimentSubset",
  slots = representation(root = "ANY",
                         subsets = "list"),
  prototype = list(subsets = list()),
  validity = function(object) {
    if (is.null(object@root)) {
      return("The root object cannot be 'NULL'.")
    } else if (!inherits(object@root, "SummarizedExperiment")) {
      return("The root slot of an 'ExperimentSubset' object can only contain an object which is inherited from 'SummarizedExperiment'.")
    } else {
      return(TRUE)
    }
  }
)

#' @title ExperimentSubset constructor
#' @description This constructor function is used to setup the
#'   \code{ExperimentSubset} object by passing either a
#'   \code{SingleCellExperiment} or \code{SummarizedExperiment} objects or
#'   objects inherited by these classes. A subset can also be directly created
#'   by passing a named \code{list} to the \code{subset} parameter. This named
#'   \code{list} should have parameter values named as \code{subsetName},
#'   \code{rows}, \code{cols} and \code{parentAssay}.
#' @param object A \code{SingleCellExperiment} or \code{SummarizedExperiment}
#'   object as the root.
#' @param subset A named \code{list} if a subset should be created from within
#'   the constructor. Named parameters in this list should be \code{subsetName},
#'   \code{rows}, \code{cols} and \code{parentAssay}.
#' @return A \code{ExperimentSubset} object.
#' @export
#' @import Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es
ExperimentSubset <- function(object,
                             subset = list(
                               subsetName = NA,
                               rows = NA,
                               cols = NA,
                               parentAssay = NA
                             ))
{
  if (!missing(object)) {
    es <- .ExperimentSubset(root = object)
  }
  else{
    stop("root object cannot be empty.")
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
#' @description Create a subset from an already available \code{assay} in the
#'   input \code{ExperimentSubset} object by specifying the rows and columns to
#'   include in the subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Specify the name of the subset to create.
#' @param rows Specify the rows to include in this subset. If \code{missing} or
#'   \code{NULL}, all rows are included in the subset. Values can be
#'   \code{numeric} or \code{character}. Default \code{NULL}.
#' @param cols Specify the columns to include in this subset. If \code{missing}
#'   or \code{NULL}, all columns are included in the subset. Values can be
#'   \code{numeric} or \code{character}. Default \code{NULL}.
#' @param parentAssay Specify the parent \code{assay} of the subset. This parent
#'   \code{assay} must already be available in the \code{ExperimentSubset}
#'   object. If \code{NULL}, the first available main \code{assay} will be
#'   marked as parent. Default \code{NULL}.
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
          || parentAssay %in% subsetAssayNames(object)) {
        tempAssay <- parentAssay
      }
      else{
        stop("Input parentAssay does not exist.")
      }
    }
    if (is.character(rows)) {
      rows <-
        match(rows, base::rownames(
          ExperimentSubset::assay(object, withDimnames = TRUE, tempAssay)
        ))
    }
    if (is.character(cols)) {
      cols <-
        match(cols, base::colnames(
          ExperimentSubset::assay(object, withDimnames = TRUE, tempAssay)
        ))
    }
    if (is.null(rows)) {
      rows <-
        seq(1, dim(
          ExperimentSubset::assay(object, withDimnames = FALSE, tempAssay)
        )[1])
    }
    if (is.null(cols)) {
      cols <-
        seq(1, dim(
          ExperimentSubset::assay(object, withDimnames = FALSE, tempAssay)
        )[2])
    }

    #Check if count of stored row/column indices greater than the subset
    if (length(rows) > dim(ExperimentSubset::assay(object, withDimnames = FALSE, tempAssay))[1]
        ||
        length(cols) > dim(ExperimentSubset::assay(object, withDimnames = FALSE, tempAssay))[2]) {
      stop("More rows or columns selected than available in the parentAssay.")
    }

    internalAssay <-
      SingleCellExperiment::SingleCellExperiment(list(
        counts = Matrix::Matrix(
          nrow = length(rows),
          ncol = length(cols),
          data = 0,
          sparse = TRUE
        )
      ))

    internalAssay <- as(internalAssay, class(object@root)[1])

    scs <- SingleCellSubset(
      subsetName = subsetName,
      rowIndices = rows,
      colIndices = cols,
      parentAssay = parentAssay,
      internalAssay = internalAssay
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
#' @description Store a new subset \code{assay} inside a specified subset in the
#'   input \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object
#' @param subsetName Specify the name of the existing subset inside which the
#'   new subset \code{assay} should be stored.
#' @param inputMatrix The input subset \code{assay}.
#' @param subsetAssayName Specify the name of the new \code{assay} against the
#'   \code{inputMatrix} parameter. If \code{NULL}, a new subset is created
#'   internally using the \code{createSubset} function. Default \code{NULL}.
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

      object <- createSubset(
        object,
        subsetName,
        base::rownames(inputMatrix),
        base::colnames(inputMatrix),
        parentAssay = NULL
      )

      object@subsets[[subsetName]]@internalAssay <-
        SingleCellExperiment::SingleCellExperiment(list(counts = inputMatrix))

    }
    else{
      SummarizedExperiment::assay(object@subsets[[subsetName]]@internalAssay,
                                  withDimnames = FALSE,
                                  subsetAssayName) <- inputMatrix
      base::rownames(object@subsets[[subsetName]]@internalAssay) <-
        base::rownames(inputMatrix)
      base::colnames(object@subsets[[subsetName]]@internalAssay) <-
        base::colnames(inputMatrix)
    }

    return(object)
  }
)

#' @title assay
#' @description Method to get an \code{assay} from an \code{ExperimentSubset}
#'   object or a subset from an \code{ExperimentSubset} object or any object
#'   supported by \code{assay} from \code{SummarizedExperiment}.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{assay} from \code{SummarizedExperiment}.
#' @param i Name of an \code{assay} or name of a subset or name of a subset
#'   \code{assay}.
#' @param withDimnames Set whether dimnames should be applied to \code{assay}.
#'   Default \code{FALSE}.
#' @param ... Additional parameters.
#' @return The \code{assay} from the input object.
#' @export
#' @importMethodsFrom SummarizedExperiment assay
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
  out <- NULL
  #look at main assays
  if (i %in% SummarizedExperiment::assayNames(x@root)) {
    out <-
      SummarizedExperiment::assay(x@root, i, withDimnames = withDimnames, ... = ...)
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
      out <- ExperimentSubset::assay(x, withDimnames = FALSE, i)
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
  if (is.null(out)) {
    stop("requested assay not found")
  }
  out
})

#' @title subsetNames
#' @description Retrieves the names of the available subsets in an
#'   \code{ExperimentSubset} object.
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

#' @title altExps
#' @description A wrapper to the \link[SingleCellExperiment]{altExps} method
#'   with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \link[SingleCellExperiment]{altExps} method.
#' @param withColData Same as \link[SingleCellExperiment]{altExps}. Default
#'   \code{FALSE}.
#' @param subsetName Specify the name of the subset from which the
#'   \code{altExps} should be fetched from. If \code{missing},
#'   \link[SingleCellExperiment]{altExps} method is called on the main object.
#' @return \code{altExps} from the specified subset or same as
#'   \link[SingleCellExperiment]{altExps} when \code{subsetName} is
#'   \code{missing}.
#' @rdname altExps
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExps(es, subsetName = "subset1") <- list(
#' alt1 = SingleCellExperiment(
#' assays = list(
#' counts = assay(es, "subset1"))),
#' alt2 = SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1"))))
#' altExps(es, subsetName = "subset1")
setGeneric(
  name = "altExps",
  def = function(x, withColData = FALSE, subsetName)
  {
    standardGeneric("altExps")
  }
)

#' @rdname altExps
setMethod(
  f = "altExps",
  signature = "ANY",
  definition = function(x, withColData, subsetName)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      SingleCellExperiment::altExps(x@subsets[[subsetName]]@internalAssay, withColData = withColData)
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::altExps(x, withColData = withColData)
      }
      else{
        SingleCellExperiment::altExps(x@root, withColData = withColData)
      }
    }
  }
)

#' @title altExp
#' @description A wrapper to the \code{altExp} from
#'   \link[SingleCellExperiment]{altExps} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{altExp} from \link[SingleCellExperiment]{altExps} method.
#' @param e Same as \code{altExp} from \link[SingleCellExperiment]{altExps}.
#' @param withColData Same as \code{altExp} from
#'   \link[SingleCellExperiment]{altExps}. Default \code{FALSE}.
#' @param subsetName Specify the name of the subset from which the \code{altExp}
#'   should be fetched from. If \code{missing}, \code{altExp} from
#'   \link[SingleCellExperiment]{altExps} method is called on the main object.
#' @return The \code{altExp} from the specified subset or same as \code{altExp}
#'   from \link[SingleCellExperiment]{altExps} when \code{subsetName} is
#'   \code{missing}.
#' @rdname altExp
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExp(es, e = "altExample",
#' subsetName = "subset1") <- SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1")))
#' altExp(es, subsetName = "subset1")
setGeneric(
  name = "altExp",
  def = function(x, e, withColData = FALSE, subsetName)
  {
    standardGeneric("altExp")
  }
)

#' @rdname altExp
setMethod(
  f = "altExp",
  signature = "ANY",
  definition = function(x, e, withColData, subsetName)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      if (missing(e)) {
        SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, withColData = withColData)
      }
      else{
        SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, e, withColData = withColData)
      }
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        if (missing(e)) {
          SingleCellExperiment::altExp(x, withColData = withColData)
        }
        else{
          SingleCellExperiment::altExp(x, e, withColData = withColData)
        }
      } else{
        if (missing(e)) {
          SingleCellExperiment::altExp(x@root, withColData = withColData)
        }
        else{
          SingleCellExperiment::altExp(x@root, e, withColData = withColData)
        }
      }
    }
  }
)

#' @title altExpNames
#' @description A wrapper to the \code{altExpNames} from
#'   \link[SingleCellExperiment]{altExps} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{altExpNames} from \link[SingleCellExperiment]{altExps} method.
#' @param subsetName Specify the name of the subset from which the
#'   \code{altExpNames} should be fetched from. If \code{missing},
#'   \code{altExpNames} from \link[SingleCellExperiment]{altExps} method is
#'   called on the main object.
#' @return The \code{altExpNames} from the specified subset or same as
#'   \code{altExpNames} from \link[SingleCellExperiment]{altExps} when
#'   \code{subsetName} is \code{missing}.
#' @rdname altExpNames
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExp(es, e = "altExample",
#' subsetName = "subset1") <- SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1")))
#' altExpNames(es, subsetName = "subset1")
setGeneric(
  name = "altExpNames",
  def = function(x, subsetName)
  {
    standardGeneric("altExpNames")
  }
)

#' @rdname altExpNames
setMethod(
  f = "altExpNames",
  signature = "ANY",
  definition = function(x, subsetName)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      SingleCellExperiment::altExpNames(x@subsets[[subsetName]]@internalAssay)
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::altExpNames(x)
      }
      else{
        SingleCellExperiment::altExpNames(x@root)
      }
    }
  }
)

#' @title reducedDimNames
#' @description A wrapper to the \code{reducedDimNames} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims}
#'   method.
#' @param subsetName Specify the name of the subset from which the
#'   \code{reducedDimNames} should be fetched from. If \code{missing},
#'   \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method
#'   is called on the main object.
#' @return The \code{reducedDimNames} from the specified subset or same as
#'   \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} when
#'   \code{subsetName} is \code{missing}.
#' @rdname reducedDimNames
#' @export
setGeneric(
  name = "reducedDimNames",
  def = function(x, subsetName)
  {
    standardGeneric("reducedDimNames")
  }
)

#' @rdname reducedDimNames
setMethod(
  f = "reducedDimNames",
  signature = "ANY",
  definition = function(x, subsetName)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      SingleCellExperiment::reducedDimNames(x@subsets[[subsetName]]@internalAssay)
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::reducedDimNames(x)
      }
      else{
        SingleCellExperiment::reducedDimNames(x@root)
      }
    }
  }
)

#' @title altExpNames<-
#' @description A wrapper to the \code{altExpNames<-} from
#'   \link[SingleCellExperiment]{altExps} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{altExpNames<-} from \link[SingleCellExperiment]{altExps} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{altExpNames<-} should be set to. If \code{missing},
#'   \code{altExpNames<-} from \link[SingleCellExperiment]{altExps} method is
#'   called on the main object.
#' @param value Input value same as \code{altExpNames<-} from
#'   \link[SingleCellExperiment]{altExps} method.
#' @return Input object with \code{altExpNames} set.
#' @rdname altExpNames-set
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExp(es, e = "altExample",
#' subsetName = "subset1") <- SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1")))
#' altExpNames(es, subsetName = "subset1") <- c("altExpSubset1")
setGeneric(
  name = "altExpNames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("altExpNames<-")
  }
)

#' @rdname altExpNames-set
setReplaceMethod(
  f = "altExpNames",
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
      SingleCellExperiment::altExpNames(x@subsets[[subsetName]]@internalAssay) <-
        value
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::altExpNames(x) <- value
      }
      else{
        SingleCellExperiment::altExpNames(x@root) <- value
      }
    }
    x
  }
)

#' @title reducedDimNames<-
#' @description A wrapper to the \code{reducedDimNames<-} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims}
#'   method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{reducedDimNames<-} should be set to. If \code{missing},
#'   \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims}
#'   method is called on the main object.
#' @param value Input value same as \code{reducedDimNames<-} from
#'   \link[SingleCellExperiment]{reducedDims} method.
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
      SingleCellExperiment::reducedDimNames(x@subsets[[subsetName]]@internalAssay) <-
        value
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::reducedDimNames(x) <- value
      }
      else{
        SingleCellExperiment::reducedDimNames(x@root) <- value
      }
    }
    x
  }
)

#' @title altExp<-
#' @description A wrapper to the \code{altExp<-} from
#'   \link[SingleCellExperiment]{altExps} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{altExp<-} from \link[SingleCellExperiment]{altExps} method.
#' @param e Same as \code{altExp<-} from \link[SingleCellExperiment]{altExps}
#'   method.
#' @param withColData Same as \code{altExp<-} from
#'   \link[SingleCellExperiment]{altExps} method. Default \code{FALSE}.
#' @param subsetName Specify the name of the subset to which the \code{altExp<-}
#'   should be set to. If \code{missing}, \code{altExp<-} from
#'   \link[SingleCellExperiment]{altExps} method is called on the main object.
#' @param value Input value same as \code{altExp<-} from
#'   \link[SingleCellExperiment]{altExps} method.
#' @return Input object with \code{altExp<-} set.
#' @rdname altExp-set
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExp(es, e = "altExample",
#' subsetName = "subset1") <- SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1")))
setGeneric(
  name = "altExp<-",
  def = function(x, e, withColData = FALSE, subsetName, value)
  {
    standardGeneric("altExp<-")
  }
)

#' @rdname altExp-set
setReplaceMethod(
  f = "altExp",
  signature = "ANY",
  definition = function(x, e, withColData, subsetName, value)
  {
    if (!missing(subsetName)) {
      if (is.null(x@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      if (missing(e)) {
        SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, withColData = withColData) <-
          value
      }
      else{
        SingleCellExperiment::altExp(x@subsets[[subsetName]]@internalAssay, e, withColData = withColData) <-
          value
      }
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        if (missing(e)) {
          SingleCellExperiment::altExp(x, withColData = withColData) <- value
        }
        else{
          SingleCellExperiment::altExp(x, e, withColData = withColData) <-
            value
        }
      }
      else{
        if (missing(e)) {
          SingleCellExperiment::altExp(x@root, withColData = withColData) <-
            value
        }
        else{
          SingleCellExperiment::altExp(x@root, e, withColData = withColData) <-
            value
        }
      }
    }
    x
  }
)

#' @title altExps<-
#' @description A wrapper to the \code{altExps<-} from
#'   \link[SingleCellExperiment]{altExps} method with additional support for
#'   subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{altExps<-} from \link[SingleCellExperiment]{altExps} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{altExps<-} should be set to. If \code{missing}, \code{altExps<-} from
#'   \link[SingleCellExperiment]{altExps} method is called on the main object.
#' @param value Input value same as \code{altExps<-} from
#'   \link[SingleCellExperiment]{altExps} method.
#' @return Input object with \code{altExps<-} set.
#' @rdname altExps-set
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' altExps(es, subsetName = "subset1") <- list(
#' alt1 = SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1"))),
#' alt2 = SingleCellExperiment(
#' assays = list(counts = assay(es, "subset1"))))
#' altExpNames(es, subsetName = "subset1")
setGeneric(
  name = "altExps<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("altExps<-")
  }
)

#' @rdname altExps-set
setReplaceMethod(
  f = "altExps",
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
      SingleCellExperiment::altExps(x@subsets[[subsetName]]@internalAssay) <-
        value
    }
    else{
      if (!inherits(x, "ExperimentSubset")) {
        SingleCellExperiment::altExps(x) <- value
      }
      else{
        SingleCellExperiment::altExps(x@root) <- value
      }
    }
    x
  }
)

#' @title metadata
#' @description Get \code{metadata} from an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to get the \code{metadata} from. If
#'   \code{missing}, \code{metadata} is fetched from the main input object.
#' @return A \code{list} of \code{metadata} elements.
#' @rdname metadata
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' metadata(es, subsetName = "subset1") <- list(
#' meta1 = "This is an example for metadata in subset1")
#' metadata(es, subsetName = "subset1")
setGeneric(
  name = "metadata",
  def = function(object, subsetName)
  {
    standardGeneric("metadata")
  }
)

#' @rdname metadata
setMethod(
  f = "metadata",
  signature = "ANY",
  definition = function(object, subsetName)
  {
    if (!missing(subsetName)) {
      if (is.null(object@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      object@subsets[[subsetName]]@internalAssay@metadata
    }
    else{
      if (!inherits(object, "ExperimentSubset")) {
        object@metadata
      }
      else{
        object@root@metadata
      }
    }
  }
)

#' @title subsetDim
#' @description Retrieves the dimensions of the specified subset in an
#'   \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to retrieve the dimensions from.
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
  def = function(object, subsetName)
  {
    standardGeneric("subsetDim")
  }
)

#' @rdname subsetDim
setMethod(
  f = "subsetDim",
  signature = c("ExperimentSubset", "character"),
  definition = function(object, subsetName)
  {
    dim(object@subsets[[subsetName]]@internalAssay)
  }
)

#' @title metadata<-
#' @description Set \code{metadata} to an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to set the \code{metadata} to. If
#'   \code{missing}, \code{metadata} is set to the main input object.
#' @param value A \code{list} to set to the \code{metadata} slot.
#' @return Input object with \code{metadata} set.
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' metadata(es, subsetName = "subset1") <- list(
#' meta1 = "This is an example for metadata in subset1")
setGeneric(
  name = "metadata<-",
  def = function(object, subsetName, value)
  {
    standardGeneric("metadata<-")
  }
)

#' @title metadata<-
#' @description Set \code{metadata} to an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to set the \code{metadata} to. If \code{missing}, \code{metadata} is set to the main input object.
#' @param value A \code{list} to set to the \code{metadata} slot.
#' @return Input object with \code{metadata} set.
#' @export
setReplaceMethod(
  f = "metadata",
  signature = "ANY",
  definition = function(object, subsetName, value)
  {
    if (!missing(subsetName)) {
      if (is.null(object@subsets[[subsetName]])) {
        stop(paste(
          subsetName,
          "does not exist in the subsets slot of the object."
        ))
      }
      object@subsets[[subsetName]]@internalAssay@metadata <- value
    }
    else{
      if (!inherits(object, "ExperimentSubset")) {
        object@metadata <- value
      }
      else{
        object@root@metadata <- value
      }
    }
    object
  }
)

#' @title subsetCount
#' @description Get the total count of the available subsets in an
#'   \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
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
  def = function(object)
  {
    standardGeneric("subsetCount")
  }
)

#' @rdname subsetCount
setMethod(
  f = "subsetCount",
  signature = "ExperimentSubset",
  definition = function(object)
  {
    return(length(subsetNames(object)))
  }
)

#' @title subsetAssayCount
#' @description Get the count of the total available subsets and the subset
#'   assays inside these subsets in an \code{ExperimentSubset} object.
#' @param object Input \code{ExperimentSubset} object.
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
  def = function(object)
  {
    standardGeneric("subsetAssayCount")
  }
)


#' @rdname subsetAssayCount
setMethod(
  f = "subsetAssayCount",
  signature = "ExperimentSubset",
  definition = function(object)
  {
    return(length(subsetAssayNames(object)))
  }
)

#' @title showSubsetLink
#' @description The function displays the content of an \code{ExperimentSubset}
#'   object including all available main assays, all subsets and the subset
#'   assays inside these subsets. This function also depicts how and in what
#'   order the subsets in the object are linked with their parents. Moreover,
#'   all supplementary data inside the subsets such as \code{reducedDims} and
#'   \code{altExps} are also displayed against each subset entry.
#' @param object Input \code{ExperimentSubset} object.
#' @return Prints all the available subset information against the input
#'   \code{ExperimentSubset} object.
#' @rdname showSubsetLink
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
#' showSubsetLink(es)
setGeneric(
  name = "showSubsetLink",
  def = function(object)
  {
    standardGeneric("showSubsetLink")
  }
)

#' @rdname showSubsetLink
setMethod(
  f = "showSubsetLink",
  signature = "ExperimentSubset",
  definition = function(object)
  {
    cat("Main assay(s):\n",
        SummarizedExperiment::assayNames(object@root),
        "\n\n")
    cat("Subset(s):\n")
    if (!is.null(subsetNames(object))) {
      Name <- list()
      Dimensions <- list()
      Parent <- list()
      Assays <- list()
      Metadata <- list()
      ReducedDims <- list()
      AltExperiments <- list()

      for (i in seq(length(subsetNames(object)))) {
        parent <- subsetParent(object, subsetAssayNames(object)[i])
        Name[[i]] <- subsetNames(object)[i]
        Parent[[i]] <-
          paste(unlist(parent), collapse = ' -> ')
        if (is.null(SummarizedExperiment::assayNames(object@subsets[[i]]@internalAssay))) {
          Assays[[i]] <- ""
        }
        else{
          Assays[[i]] <-
            SummarizedExperiment::assayNames(object@subsets[[i]]@internalAssay)
        }
        Dimensions[[i]] <-
          paste(unlist(subsetDim(object, subsetNames(object)[i])), collapse = ', ')
        ReducedDims[[i]] <-
          paste(unlist(reducedDimNames(object, subsetNames(object)[i])), collapse = ", ")
        AltExperiments[[i]] <-
          paste(unlist(altExpNames(object, subsetName = subsetNames(object)[i])), collapse = ", ")
      }

      df <- data.frame(
        Name = as.character(Name),
        Dim = as.character(Dimensions),
        Parent = as.character(Parent)
      )

      if (length(which(as.character(Assays) == "")) != subsetCount(object)) {
        df <- cbind(df, Assays = as.character(Assays))
      }

      if (length(which(as.character(AltExperiments) == "")) != subsetCount(object)) {
        df <- cbind(df, AltExperiments = as.character(AltExperiments))
      }

      if (length(which(as.character(ReducedDims) == "")) != subsetCount(object)) {
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
#' @description Retrieves a complete subset to parent link from a specified
#'   subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Specify the name of the subset against which the subset to
#'   parent link should be retrieved.
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
  def = function(object, subsetName)
  {
    standardGeneric("subsetParent")
  }
)

#' @rdname subsetParent
setMethod(
  f = "subsetParent",
  signature = "ANY",
  definition = function(object, subsetName)
  {
    parentList <- list()
    if (!subsetName %in% subsetAssayNames(object)) {
      stop(paste(subsetName, "does not exist in the subsets slot of the object."))
    }
    if (!is.null(object@subsets[[subsetName]]) &&
        is.null(object@subsets[[subsetName]]@parentAssay)) {
      return(NULL)
    }
    parent <- subsetName
    while (TRUE) {
      parentList <- c(parentList, parent)
      if (!is.null(object@subsets[[parent]])) {
        parent <- object@subsets[[parent]]@parentAssay
      }
      else{
        for (i in seq(subsetCount(object))) {
          if (parent %in% SummarizedExperiment::assayNames(object@subsets[[i]]@internalAssay)) {
            parent <- object@subsets[[i]]@subsetName
          }
        }
        parentList <- c(parentList, parent)
        parent <- object@subsets[[parent]]@parentAssay
      }
      if (parent %in% SummarizedExperiment::assayNames(object@root)) {
        parentList <- c(parentList, parent)
        break
      }
    }
    parentList[[1]] <- NULL
    return(parentList)
  }
)

#' @title rownames
#' @description Get \code{rownames} from an \code{ExperimentSubset} object or a
#'   subset in the \code{ExperimentSubset} object or any object supported by
#'   \code{rownames} in \code{base} package.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{rownames} in \code{base} package.
#' @param ... Additional parameters and \code{subsetName} parameter to pass the
#'   name of the subset to get \code{rownames} from.
#' @param subsetName Name of the subset to get \code{rownames} from. If
#'   \code{missing}, \code{rownames} from main object are returned.
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
  def = function(object, ...)
  {
    standardGeneric("rownames")
  }
)

#' @rdname rownames
setMethod(
  f = "rownames",
  signature = "ANY",
  definition = function(object, subsetName, ...)
  {
    if (!inherits(object, "ExperimentSubset")) {
      base::rownames(object, ...)
    }
    else if (missing(subsetName)) {
      base::rownames(object@root, ...)
    }
    else{
      if (subsetName %in% subsetNames(object)) {
        esRownames <-
          base::rownames(object@root, ...)[object@subsets[[subsetName]]@rowIndices]
        subsetRownames <-
          base::rownames(object@subsets[[subsetName]]@internalAssay, ...)
        if (is.null(subsetRownames)) {
          subsetRownames <- esRownames
        }
        subsetRownames
      }
      else if (subsetName %in% subsetAssayNames(object)) {
        subsetName <- .getParentAssayName(object, subsetName)
        base::rownames(object@subsets[[subsetName]]@internalAssay, ...)
      }
      else{
        NULL
      }
    }
  }
)

#' @title rownames<-
#' @description Set \code{rownames} to an \code{ExperimentSubset} object or a
#'   subset in the \code{ExperimentSubset} object or any object supported by
#'   \code{rownames} in \code{base} package.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{rownames} in \code{base} package.
#' @param subsetName Name of the subset to get \code{rownames} from. If
#'   \code{missing}, \code{rownames} from main object are returned.
#' @param ... Additional parameters and \code{subsetName} parameter to pass the
#'   name of the subset to get \code{rownames} from.
#' @param value A \code{list} of \code{rownames} to set to the input object.
#' @return Input object with \code{rownames} set.
#' @rdname rownames-set
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' rownames(es, subsetName = "subset1") <-
#' paste0("row", seq(subsetDim(es, subsetName = "subset1")[1]))
setGeneric(
  name = "rownames<-",
  def = function(object, ..., value)
  {
    standardGeneric("rownames<-")
  }
)

#' @rdname rownames-set
setReplaceMethod(
  f = "rownames",
  signature = "ANY",
  definition = function(object, subsetName, ..., value)
  {
    if (!inherits(object, "ExperimentSubset")) {
      base::rownames(object, ...) <- value
    }
    else if (missing(subsetName)) {
      base::rownames(object@root, ...) <- value
    }
    else{
      if (subsetName %in% subsetNames(object)) {
        base::rownames(object@subsets[[subsetName]]@internalAssay, ...) <-
          value
      }
      else if (subsetName %in% subsetAssayNames(object)) {
        subsetName <- .getParentAssayName(object, subsetName)
        base::rownames(object@subsets[[subsetName]]@internalAssay, ...) <-
          value
      }
    }
    object
  }
)

#' @title colnames
#' @description Get \code{colnames} from an \code{ExperimentSubset} object or a
#'   subset in the \code{ExperimentSubset} object or any object supported by
#'   \code{colnames} in \code{base} package.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{colnames} in \code{base} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If
#'   \code{missing}, \code{colnames} from main object are returned.
#' @param ... Additional parameters and \code{subsetName} parameter to pass the
#'   name of the subset to get \code{colnames} from.
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
  def = function(object, ...)
  {
    standardGeneric("colnames")
  }
)

#' @rdname colnames
setMethod(
  f = "colnames",
  signature = "ANY",
  definition = function(object, subsetName, ...)
  {
    if (!inherits(object, "ExperimentSubset")) {
      base::colnames(object, ...)
    }
    else if (missing(subsetName)) {
      base::colnames(object@root, ...)
    }
    else{
      if (subsetName %in% subsetNames(object)) {
        esColnames <-
          base::colnames(object@root, ...)[object@subsets[[subsetName]]@colIndices]
        subsetColnames <-
          base::colnames(object@subsets[[subsetName]]@internalAssay, ...)
        if (is.null(subsetColnames)) {
          subsetColnames <- esColnames
        }
        subsetColnames
      }
      else if (subsetName %in% subsetAssayNames(object)) {
        subsetName <- .getParentAssayName(object, subsetName)
        base::colnames(object@subsets[[subsetName]]@internalAssay, ...)
      }
      else{
        NULL
      }
    }
  }
)

#' @title colnames<-
#' @description Set \code{colnames} to an \code{ExperimentSubset} object or a
#'   subset in the \code{ExperimentSubset} object or any object supported by
#'   \code{colnames} in \code{base} package.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{colnames} in \code{base} package.
#' @param subsetName Name of the subset to get \code{colnames} from. If
#'   \code{missing}, \code{colnames} from main object are returned.
#' @param ... Additional parameters and \code{subsetName} parameter to pass the
#'   name of the subset to get \code{colnames} from.
#' @param value A \code{list} of \code{colnames} to set to the input object.
#' @return Input object with \code{colnames} set.
#' @rdname colnames-set
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' colnames(es, subsetName = "subset1") <-
#' paste0("col", seq(subsetDim(es, subsetName = "subset1")[2]))
setGeneric(
  name = "colnames<-",
  def = function(object, ..., value)
  {
    standardGeneric("colnames<-")
  }
)

#' @rdname colnames-set
setReplaceMethod(
  f = "colnames",
  signature = "ANY",
  definition = function(object, subsetName, ..., value)
  {
    if (!inherits(object, "ExperimentSubset")) {
      base::colnames(object, ...) <- value
    }
    else if (missing(subsetName)) {
      base::colnames(object@root, ...) <- value
    }
    else{
      if (subsetName %in% subsetNames(object)) {
        base::colnames(object@subsets[[subsetName]]@internalAssay, ...) <-
          value
      }
      else if (subsetName %in% subsetAssayNames(object)) {
        subsetName <- .getParentAssayName(object, subsetName)
        base::colnames(object@subsets[[subsetName]]@internalAssay, ...) <-
          value
      }
    }
    object
  }
)

#' @title subsetAssayNames
#' @description Retrieves the names of all the subsets as well as the subset
#'   assays.
#' @param object Input \code{ExperimentSubset} object.
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
  signature = "ExperimentSubset",
  definition = function(object)
  {
    tempNames <- names(object@subsets)
    if (length(object@subsets) > 0) {
      for (i in seq(length(object@subsets))) {
        tempNames <-
          c(
            tempNames,
            SummarizedExperiment::assayNames(object@subsets[[i]]@internalAssay)
          )
      }
    }
    return(tempNames)
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
  signature = "ExperimentSubset",
  definition = function(object)
  {
    cat("class: ExperimentSubset\n")
    cat("root ")
    cat(SingleCellExperiment::show(object@root))
    cat("subsets(", length(subsetNames(object)), "): ",
        sep = "")
    cat(subsetNames(object))
    cat("\nsubsetAssays(", length(subsetAssayNames(object)), "): ",
        sep = "")
    cat(subsetAssayNames(object))
  }
)

#' @title assay<-
#' @description Method to set an \code{assay} to an \code{ExperimentSubset}
#'   object or a subset from an \code{ExperimentSubset} object or any object
#'   supported by \code{assay<-} from \code{SummarizedExperiment}.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{assay} from \code{SummarizedExperiment}.
#' @param i Name of an \code{assay} or name of the subset if storing to an
#'   \code{ExperimentSubset} object.
#' @param withDimnames Set whether dimnames should be applied to \code{assay}.
#'   Default \code{FALSE}.
#' @param subsetAssayName Name of the assay to store if storing to an
#'   \code{ExperimentSubset} object.
#' @param ... Additional parameters.
#' @param value The \code{assay} to store.
#' @return Input object with \code{assay} stored.
#' @export
#' @importMethodsFrom SummarizedExperiment assay<-
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
setReplaceMethod("assay",
                 c("ExperimentSubset", "character"),
                 function(x,
                          i,
                          withDimnames = FALSE,
                          subsetAssayName = NULL,
                          ...,
                          value) {
                   if ((nrow(value) != nrow(x@root))
                       || (ncol(value) != ncol(x@root))) {
                     x <- storeSubset(
                       object = x,
                       subsetName = i,
                       inputMatrix = value,
                       subsetAssayName = subsetAssayName
                     )
                   }
                   else{
                     SummarizedExperiment::assay(
                       x = x@root,
                       i = i,
                       withDimnames = withDimnames,
                       ... = ...
                     ) <- value
                   }
                   x
                 })

#' @title subsetRowData
#' @description Get \code{rowData} from a subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to get \code{rowData} from.
#' @return The \code{rowData} from input object.
#' @rdname subsetRowData
#' @export
setGeneric(
  name = "subsetRowData",
  def = function(object, subsetName)
  {
    standardGeneric("subsetRowData")
  }
)

#' @rdname subsetRowData
setMethod(
  f = "subsetRowData",
  signature = c("ExperimentSubset", "character"),
  definition = function(object, subsetName)
  {
    if (subsetName %in% subsetNames(object)) {
      #is a subset
      out <-
        SummarizedExperiment::rowData(object@root)[object@subsets[[subsetName]]@rowIndices, , drop = F]
      out <-
        cbind(out, rowData(object@subsets[[subsetName]]@internalAssay))
    }
    else if (subsetName %in% subsetAssayNames(object)) {
      #is a subset assay
      subsetName <- .getParentAssayName(object, subsetName)
      out <-
        SummarizedExperiment::rowData(object@root)[object@subsets[[subsetName]]@rowIndices, , drop = F]
      out <-
        cbind(out, rowData(object@subsets[[subsetName]]@internalAssay))
    }
    else{
      #neither a subset nor a subset assay
      stop("Neither a subset nor a subsetAssay.")
    }
    return(out)
  }
)

.getParentAssayName <- function(object, childAssayName) {
  for (i in seq(length(object@subsets))) {
    if (childAssayName %in% SummarizedExperiment::assayNames(object@subsets[[i]]@internalAssay)) {
      return(object@subsets[[i]]@subsetName)
    }
  }
}

#' @title subsetColData
#' @description Get \code{colData} from a subset.
#' @param object Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to get \code{colData} from.
#' @return The \code{colData} from input object.
#' @rdname subsetColData
#' @export
setGeneric(
  name = "subsetColData",
  def = function(object, subsetName)
  {
    standardGeneric("subsetColData")
  }
)

#' @rdname subsetColData
setMethod(
  f = "subsetColData",
  signature = c("ExperimentSubset", "character"),
  definition = function(object, subsetName)
  {
    if (subsetName %in% subsetNames(object)) {
      #is a subset
      out <-
        SummarizedExperiment::colData(object@root)[object@subsets[[subsetName]]@colIndices, , drop = F]
      out <-
        cbind(out, colData(object@subsets[[subsetName]]@internalAssay))
    }
    else if (subsetName %in% subsetAssayNames(object)) {
      #is a subset assay
      subsetName <- .getParentAssayName(object, subsetName)
      out <-
        SummarizedExperiment::colData(object@root)[object@subsets[[subsetName]]@colIndices, , drop = F]
      out <-
        cbind(out, colData(object@subsets[[subsetName]]@internalAssay))
    }
    else{
      #neither a subset nor a subset assay
      stop("Neither a subset nor a subsetAssay.")
    }
    return(out)
  }
)

#' @title reducedDim
#' @description A wrapper to the \code{reducedDim} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDim} from \link[SingleCellExperiment]{reducedDims} method.
#' @param type Same as \code{type} in \code{reducedDim} from
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param withDimnames Same as \code{withDimnames} in \code{reducedDim} from
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset from which the
#'   \code{reducedDim} should be fetched from. If \code{missing},
#'   \code{reducedDim} from \link[SingleCellExperiment]{reducedDims} method is
#'   called on the main object.
#' @return The \code{reducedDim} from the specified subset or same as
#'   \code{reducedDim} from \link[SingleCellExperiment]{reducedDims} when
#'   \code{subsetName} is \code{missing}.
#' @rdname reducedDim
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:1500), cols = c(1:1500),
#' parentAssay = "counts")
#' reducedDim(es, type = "PCA",
#' subsetName = "subset1") <- scater::calculatePCA(
#' assay(es, "subset1"))
#' reducedDim(es, type = "PCA", subsetName = "subset1")
setGeneric(
  name = "reducedDim",
  def = function(object, type, withDimnames, subsetName)
  {
    standardGeneric("reducedDim")
  }
)

#' @rdname reducedDim
setMethod("reducedDim", "ANY", function(object, type, withDimnames, subsetName) {
  if (missing(withDimnames)) {
    withDimnames = FALSE
  }
  if (!missing(subsetName)) {
    out <-
      SingleCellExperiment::reducedDim(object@subsets[[subsetName]]@internalAssay, type, withDimnames)
  }
  else{
    if (!inherits(object, "ExperimentSubset")) {
      out <- SingleCellExperiment::reducedDim(object, type, withDimnames)
    }
    else{
      out <-
        SingleCellExperiment::reducedDim(object@root, type, withDimnames)
    }
  }
  out
})

#' @title reducedDims
#' @description A wrapper to the \link[SingleCellExperiment]{reducedDims} method
#'   with additional support for subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param withDimnames Same as \code{withDimnames} in
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset from which the
#'   \code{reducedDims} should be fetched from. If \code{missing},
#'   \link[SingleCellExperiment]{reducedDims} method is called on the main
#'   object.
#' @return The \code{reducedDims} from the specified subset or same as
#'   link[SingleCellExperiment]{reducedDims} when \code{subsetName} is
#'   \code{missing}.
#' @rdname reducedDims
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:1500), cols = c(1:1500),
#' parentAssay = "counts")
#' reducedDim(es, type = "PCA_1",
#' subsetName = "subset1") <- scater::calculatePCA(
#' assay(es, "subset1"))
#' reducedDim(es, type = "PCA_2",
#' subsetName = "subset1") <- scater::calculatePCA(
#' assay(es, "subset1"))
#' reducedDims(es, subsetName = "subset1")
setGeneric(
  name = "reducedDims",
  def = function(object, withDimnames, subsetName)
  {
    standardGeneric("reducedDims")
  }
)

#' @rdname reducedDims
setMethod("reducedDims", "ANY", function(object, withDimnames, subsetName) {
  if (missing(withDimnames)) {
    withDimnames = FALSE
  }
  if (!missing(subsetName)) {
    out <-
      SingleCellExperiment::reducedDims(object@subsets[[subsetName]]@internalAssay, withDimnames)
  }
  else{
    if (!inherits(object, "ExperimentSubset")) {
      out <- SingleCellExperiment::reducedDims(object, withDimnames)
    }
    else{
      out <- SingleCellExperiment::reducedDims(object@root, withDimnames)
    }
  }
  out
})

#' @title reducedDim<-
#' @description A wrapper to the \code{reducedDim<-} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDim<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param type Same as \code{type} in \code{reducedDim<-} from
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{reducedDim} should be set to. If \code{missing}, \code{reducedDim<-}
#'   from \link[SingleCellExperiment]{reducedDims} method is called on the main
#'   object.
#' @param value Value to set to \code{reducedDim}.
#' @return Updated input object with \code{reducedDim} set.
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:1500), cols = c(1:1500),
#' parentAssay = "counts")
#' reducedDim(es, type = "PCA",
#' subsetName = "subset1") <- scater::calculatePCA(
#' assay(es, "subset1"))
setGeneric(
  name = "reducedDim<-",
  def = function(object, type, subsetName, value)
  {
    standardGeneric("reducedDim<-")
  }
)

#' @title reducedDim<-
#' @description A wrapper to the \code{reducedDim<-} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDim<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param type Same as \code{type} in \code{reducedDim<-} from
#'   \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{reducedDim} should be set to. If \code{missing}, \code{reducedDim<-}
#'   from \link[SingleCellExperiment]{reducedDims} method is called on the main
#'   object.
#' @param value Value to set to \code{reducedDim}.
#' @return Updated input object with \code{reducedDim} set.
#' @export
setReplaceMethod("reducedDim", "ANY", function(object, type, subsetName, value) {
  if (!missing(subsetName)) {
    SingleCellExperiment::reducedDim(object@subsets[[subsetName]]@internalAssay, type) <-
      value
  }
  else{
    if (!inherits(object, "ExperimentSubset")) {
      SingleCellExperiment::reducedDim(object, type) <- value
    }
    else{
      SingleCellExperiment::reducedDim(object@root, type) <- value
    }
  }
  return(object)
})

#' @title reducedDims<-
#' @description A wrapper to the \code{reducedDims<-} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{reducedDims} should be set to. If \code{missing},
#'   \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method
#'   is called on the main object.
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
  def = function(object, subsetName, value)
  {
    standardGeneric("reducedDims<-")
  }
)

#' @title reducedDims<-
#' @description A wrapper to the \code{reducedDims<-} from
#'   \link[SingleCellExperiment]{reducedDims} method with additional support for
#'   subsets.
#' @param object Input \code{ExperimentSubset} object or any object supported by
#'   \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the
#'   \code{reducedDims} should be set to. If \code{missing},
#'   \code{reducedDims<-} from \link[SingleCellExperiment]{reducedDims} method
#'   is called on the main object.
#' @param value A \code{list} of values to set to \code{reducedDims}.
#' @return Updated input object with \code{reducedDims} set.
#' @export
setReplaceMethod("reducedDims", "ANY", function(object, subsetName, value) {
  if (!missing(subsetName)) {
    SingleCellExperiment::reducedDims(object@subsets[[subsetName]]@internalAssay) <-
      value
  }
  else{
    if (!inherits(object, "ExperimentSubset")) {
      SingleCellExperiment::reducedDims(object) <- value
    }
    else{
      SingleCellExperiment::reducedDims(object@root) <- value
    }
  }
  return(object)
})

#' @title rowData
#' @description Get \code{rowData} from a subset of an input object or the
#'   object itself.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{rowData} from \code{SummarizedExperiment}.
#' @param subsetName Name of the subset to get \code{rowData} from. If
#'   \code{NULL} or \code{missing}, \code{rowData} from main input object is
#'   fetched.
#' @param ... Additional parameters.
#' @return The \code{rowData} from input object or subset of an input object.
#' @export
#' @importMethodsFrom SummarizedExperiment rowData
setMethod("rowData", c("ExperimentSubset"), function(x, subsetName = NULL, ...) {
  if (!is.null(subsetName)) {
    out <- subsetRowData(object = x,
                         subsetName = subsetName)
  }
  else{
    out <- SummarizedExperiment::rowData(x = x@root,
                                         ... = ...)
  }
  out
})

#' @title colData
#' @description Get \code{colData} from a subset of an input object or the
#'   object itself.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{colData} from \code{SummarizedExperiment}.
#' @param subsetName Name of the subset to get \code{colData} from. If
#'   \code{NULL} or \code{missing}, \code{colData} from main input object is
#'   fetched.
#' @param ... Additional parameters.
#' @return The \code{colData} from input object or subset of an input object.
#' @export
#' @importMethodsFrom SummarizedExperiment colData
setMethod("colData", c("ExperimentSubset"), function(x, subsetName = NULL, ...) {
  if (!is.null(subsetName)) {
    out <- subsetColData(object = x,
                         subsetName = subsetName)
  }
  else{
    out <- SummarizedExperiment::colData(x = x@root,
                                         ... = ...)
  }
  out
})

#' @title rowData<-
#' @description Set \code{rowData} to a subset of an input object or the object
#'   itself.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{rowData<-} from \code{SummarizedExperiment}.
#' @param ... Additional parameters.
#' @param subsetName Name of the subset to set \code{rowData} to. If \code{NULL}
#'   or \code{missing}, \code{rowData} to main input object is set.
#' @param value The \code{rowData} to store in an object or subset of an object.
#' @return Object with \code{rowData} set.
#' @export
#' @importMethodsFrom SummarizedExperiment rowData<-
setReplaceMethod("rowData", c("ExperimentSubset"), function(x, ..., subsetName, value) {
  #test if this needs DataFrame too
  tempValue <- NULL
  if (!missing(subsetName)) {
    tempValue <- rowData(x@root)
    rowData(x@subsets[[subsetName]]@internalAssay) <-
      value
  }
  else{
    tempValue <- value
  }
  value <- tempValue
  SummarizedExperiment::rowData(x = x@root,
                                ... = ...) <- value
  x
})

#' @title colData<-
#' @description Set \code{colData} to a subset of an input object or the object
#'   itself.
#' @param x Input \code{ExperimentSubset} object or any object supported by
#'   \code{colData<-} from \code{SummarizedExperiment}.
#' @param ... Additional parameters.
#' @param subsetName Name of the subset to set \code{colData} to. If \code{NULL}
#'   or \code{missing}, \code{colData} to main input object is set.
#' @param value The \code{colData} to store in an object or subset of an object.
#' @return Object with \code{colData} set.
#' @export
#' @importMethodsFrom SummarizedExperiment colData<-
setReplaceMethod("colData", c("ExperimentSubset" , "DataFrame"), function(x, ..., subsetName, value) {
  tempValue <- NULL
  if (!missing(subsetName)) {
    tempValue <- colData(x@root)
    colData(x@subsets[[subsetName]]@internalAssay) <-
      value
  }
  else{
    tempValue <- value
  }
  value <- tempValue
  SummarizedExperiment::colData(x = x@root,
                                ... = ...) <- value
  x
})

#' @title dim
#' @description Get dimensions of the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object.
#' @return A \code{list} containing number of rows and number of columns of the
#'   object.
#' @export
#' @importMethodsFrom SingleCellExperiment dim
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' dim(es)
setMethod("dim",
          c("ExperimentSubset"),
          function(x) {
            return(dim(x@root))
          })
