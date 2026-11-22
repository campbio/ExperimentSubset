# The methods the Subset* classes share, each defined once in a loop over
# the class names in .esClasses (all five classes) or .esSCEClasses (the
# three that extend SingleCellExperiment). Each method is a thin wrapper
# around a helper in R/ExperimentSubset.R. Methods for one class only are
# in R/SubsetSpatialExperiment.R and R/SubsetTreeSummarizedExperiment.R.
#
# roxygen can't see methods defined in a loop, so each help page that
# documents them has a block below that generates their aliases and usage
# with @eval .rdMethods(). A new method in a loop needs its generic named
# in the block for its help page.
#
# Method bodies must not use the loop variable cls: the methods don't
# capture it, and it is removed after the loops. A method that needs its
# class can call class(x).

# Methods for all five classes ---------------------------------------------

#' @name subsetColnames
#' @rdname subsetColnames
#' @eval .rdMethods(c("subsetColnames", "subsetColnames<-"))
NULL

#' @name subsetRownames
#' @rdname subsetRownames
#' @eval .rdMethods(c("subsetRownames", "subsetRownames<-"))
NULL

#' @name subsetColData
#' @rdname subsetColData
#' @eval .rdMethods(c("subsetColData", "subsetColData<-"))
NULL

#' @name subsetRowData
#' @rdname subsetRowData
#' @eval .rdMethods(c("subsetRowData", "subsetRowData<-"))
NULL

#' @name getSubsetAssay
#' @rdname getSubsetAssay
#' @eval .rdMethods("getSubsetAssay")
NULL

#' @name getSubsetExperiment
#' @rdname getSubsetExperiment
#' @eval .rdMethods("getSubsetExperiment")
NULL

#' @name setSubsetExperiment
#' @rdname setSubsetExperiment
#' @eval .rdMethods("setSubsetExperiment")
NULL

#' @name ExperimentSubset-metadata
#' @rdname ExperimentSubset-metadata
#' @eval .rdMethods(c("metadata", "metadata<-"))
NULL

#' @name ExperimentSubset-assay
#' @rdname ExperimentSubset-assay
#' @eval .rdMethods(c("assay", "assay<-"))
NULL

#' @name createSubset
#' @rdname createSubset
#' @eval .rdMethods("createSubset")
NULL

#' @name subsetAssayNames
#' @rdname subsetAssayNames
#' @eval .rdMethods("subsetAssayNames")
NULL

#' @name SubsetSummarizedExperiment-class
#' @rdname SubsetSummarizedExperiment-class
#' @param object The object to print, with its subsets and their assays.
#' @return \code{show()} prints the object and returns \code{NULL}
#'   invisibly.
#' @eval .rdMethods("show", "SubsetSummarizedExperiment")
NULL

#' @name SubsetRangedSummarizedExperiment-class
#' @rdname SubsetRangedSummarizedExperiment-class
#' @param object The object to print, with its subsets and their assays.
#' @return \code{show()} prints the object and returns \code{NULL}
#'   invisibly.
#' @eval .rdMethods("show", "SubsetRangedSummarizedExperiment")
NULL

#' @name SubsetSingleCellExperiment-class
#' @rdname SubsetSingleCellExperiment-class
#' @param object The object to print, with its subsets and their assays.
#' @return \code{show()} prints the object and returns \code{NULL}
#'   invisibly.
#' @eval .rdMethods("show", "SubsetSingleCellExperiment")
NULL

#' @name SubsetTreeSummarizedExperiment-class
#' @rdname SubsetTreeSummarizedExperiment-class
#' @param object The object to print, with its subsets and their assays.
#' @return \code{show()} prints the object and returns \code{NULL}
#'   invisibly.
#' @eval .rdMethods("show", "SubsetTreeSummarizedExperiment")
NULL

#' @name SubsetSpatialExperiment-class
#' @rdname SubsetSpatialExperiment-class
#' @param object The object to print, with its subsets and their assays.
#' @return \code{show()} prints the object and returns \code{NULL}
#'   invisibly.
#' @eval .rdMethods("show", "SubsetSpatialExperiment")
NULL

#' @name subsetNames
#' @rdname subsetNames
#' @eval .rdMethods("subsetNames")
NULL

#' @name setSubsetAssay
#' @rdname setSubsetAssay
#' @eval .rdMethods("setSubsetAssay")
NULL

#' @name subsetParent
#' @rdname subsetParent
#' @eval .rdMethods("subsetParent")
NULL

#' @name subsetSummary
#' @rdname subsetSummary
#' @eval .rdMethods("subsetSummary")
NULL

#' @name subsetDim
#' @rdname subsetDim
#' @eval .rdMethods("subsetDim")
NULL

#' @name subsetCount
#' @rdname subsetCount
#' @eval .rdMethods("subsetCount")
NULL

#' @name subsetAssayCount
#' @rdname subsetAssayCount
#' @eval .rdMethods("subsetAssayCount")
NULL

#' @name ExperimentSubset-rowData-colData
#' @rdname ExperimentSubset-rowData-colData
#' @eval .rdMethods(c("rowData", "colData", "rowData<-", "colData<-"))
NULL

#' @name ExperimentSubset-subsetting
#' @rdname ExperimentSubset-subsetting
#' @eval .rdMethods("[")
NULL

#' @name ExperimentSubset-combine
#' @rdname ExperimentSubset-combine
#' @importMethodsFrom SummarizedExperiment cbind rbind
#' @importMethodsFrom S4Vectors combineCols combineRows
#' @eval .rdMethods(c("cbind", "rbind", "combineCols", "combineRows"))
NULL

for (cls in .esClasses) {
  setMethod(
    f = "subsetColnames",
    signature = signature(x = cls),
    definition = function(x, subsetName) {
      .subsetColnames(x, subsetName)
    }
  )

  setReplaceMethod(
    f = "subsetColnames",
    signature = signature(x = cls),
    definition = function(x, subsetName, value) {
      .subsetColnames(x, subsetName) <- value
      return(x)
    }
  )

  setMethod(
    f = "subsetRownames",
    signature = signature(x = cls),
    definition = function(x, subsetName) {
      .subsetRownames(x, subsetName)
    }
  )

  setReplaceMethod(
    f = "subsetRownames",
    signature = signature(x = cls),
    definition = function(x, subsetName, value) {
      .subsetRownames(x, subsetName) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "subsetColData",
    signature = signature(x = cls),
    definition = function(x, subsetName, value) {
      .subsetColData(x, subsetName) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "subsetRowData",
    signature = signature(x = cls),
    definition = function(x, subsetName, value) {
      .subsetRowData(x, subsetName) <- value
      return(x)
    }
  )

  setMethod(
    f = "getSubsetAssay",
    signature = signature(x = cls),
    definition = function(x,
                          subsetName,
                          subsetAssayName = NULL) {
      # so a bad subsetName isn't reported as assay()'s i (#85)
      .checkSubsetName(subsetName)
      if (is.null(subsetAssayName))
        return(assay(x = x, i = subsetName))
      # checked here so the error names subsetName, not assay()'s i (#76)
      .isSubset(x, subsetName)
      assay(x = x, i = subsetName, subsetAssayName = subsetAssayName)
    }
  )

  setMethod(
    f = "getSubsetExperiment",
    signature = signature(x = cls),
    definition = function(x,
                          subsetName,
                          parentData = TRUE) {
      .getSubsetExperiment(x, subsetName, parentData = parentData)
    }
  )

  setMethod(
    f = "setSubsetExperiment",
    signature = signature(x = cls),
    definition = function(x,
                          subsetName,
                          value,
                          unlink = FALSE) {
      .setSubsetExperiment(x, subsetName, value, unlink = unlink)
    }
  )

  setMethod(
    f = "metadata",
    signature = signature(x = cls),
    definition = function(x, ...) {
      .metadata(x, ...)
    }
  )

  setMethod(
    f = "assay",
    signature = signature(x = cls, i = "character"),
    definition = function(x, i, withDimnames = TRUE, ...) {
      .assay(x, i, withDimnames = withDimnames, ...)
    }
  )

  setMethod(
    f = "createSubset",
    signature = c(
      cls
    ),
    definition = function(x,
                          subsetName,
                          rows,
                          cols,
                          parentAssay,
                          altExp) {
      return(.createSubset(x,
                           subsetName,
                           rows,
                           cols,
                           parentAssay,
                           altExp))
    }
  )

  setMethod(
    f = "subsetAssayNames",
    signature = cls,
    definition = function(x, subsetName = NULL) {
      .subsetAssayNames(x, subsetName)
    }
  )

  setMethod(
    f = "show",
    signature = cls,
    definition = function(object) {
      .show(object)
    }
  )

  setMethod(
    f = "subsetNames",
    signature = cls,
    definition = function(x) {
      return(names(.subsets(x)))
    }
  )

  setMethod(
    f = "setSubsetAssay",
    signature = signature(x = cls),
    definition = function(x,
                          subsetName,
                          inputMatrix,
                          subsetAssayName,
                          unlink = FALSE) {
      .setSubsetAssay(x,
                      subsetName,
                      inputMatrix,
                      subsetAssayName,
                      unlink = unlink)
    }
  )

  setReplaceMethod("assay",
                   c(cls, "character"),
                   function(x,
                            i,
                            withDimnames = TRUE,
                            ...,
                            value) {
                     .assay(x, i, withDimnames = withDimnames, ...) <- value
                     return(x)
                   })

  setMethod(
    f = "subsetParent",
    signature = cls,
    definition = function(x, subsetName) {
      .subsetParent(x, subsetName)
    }
  )

  setMethod(
    f = "subsetSummary",
    signature = cls,
    definition = function(x) {
      .subsetSummary(x)
    }
  )

  setMethod(
    f = "subsetDim",
    signature = signature(x = cls),
    definition = function(x, subsetName) {
      .isSubset(x, subsetName)
      dim(.internalAssay(.subsets(x)[[subsetName]]))
    }
  )

  setMethod(
    f = "subsetCount",
    signature = cls,
    definition = function(x) {
      return(length(subsetNames(x)))
    }
  )

  setReplaceMethod(
    f = "metadata",
    signature = signature(x = cls),
    definition = function(x, ..., value) {
      .metadata(x, ...) <- value
      return(x)
    }
  )

  setMethod(
    f = "subsetAssayCount",
    signature = cls,
    definition = function(x) {
      return(length(subsetAssayNames(x)))
    }
  )

  setMethod(
    f = "subsetRowData",
    signature = signature(x = cls),
    definition = function(x, subsetName, parentRowData) {
      .subsetRowData(x, subsetName, parentRowData)
    }
  )

  setMethod(
    f = "subsetColData",
    signature = signature(x = cls),
    definition = function(x, subsetName, parentColData) {
      .subsetColData(x, subsetName, parentColData)
    }
  )

  setMethod(
    f = "rowData",
    signature = signature(x = cls),
    definition = function(x, use.names = TRUE, ...) {
      .rowData(x, use.names = use.names, ...)
    }
  )

  setMethod(
    f = "colData",
    signature = signature(x = cls),
    definition = function(x, ...) {
      .colData(x, ...)
    }
  )

  setReplaceMethod(
    f = "rowData",
    signature = signature(x = cls),
    definition = function(x, ..., value) {
      .rowData(x, ...) <- value
      return(x)
    }
  )

  # DataFrame, NULL and data.frame keep their own methods: with only the
  # "ANY" one, S4 would pick the parent class's methods for DataFrame and
  # NULL, which ignore subsetName, and would do the same for any method
  # for data.frame that a package adds to a parent class (#99)
  setReplaceMethod(
    f = "colData",
    signature = signature(x = cls, value = "DataFrame"),
    definition = function(x, ..., value) {
      .colData(x, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "colData",
    signature = signature(x = cls, value = "NULL"),
    definition = function(x, ..., value) {
      .colData(x, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "colData",
    signature = signature(x = cls, value = "data.frame"),
    definition = function(x, ..., value) {
      .colData(x, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "colData",
    signature = signature(x = cls, value = "ANY"),
    definition = function(x, ..., value) {
      # every value that reaches this method is an error from .colData<-
      `.colData<-`(x, ..., value = value)
    }
  )

  setMethod(
    f = "[",
    signature = signature(x = cls),
    definition = function(x, i, j, ..., drop = TRUE) {
      out <- callNextMethod()
      .remapSubsets(x, out, i, j)
    }
  )

  setMethod(
    f = "cbind",
    signature = cls,
    definition = function(..., deparse.level = 1) {
      objects <- list(...)
      .checkCombine(objects, "cols")
      out <- callNextMethod()
      .combineSubsets(out, objects, "cols")
    }
  )

  setMethod(
    f = "rbind",
    signature = cls,
    definition = function(..., deparse.level = 1) {
      objects <- list(...)
      .checkCombine(objects, "rows")
      out <- callNextMethod()
      .combineSubsets(out, objects, "rows")
    }
  )

  setMethod(
    f = "combineCols",
    signature = cls,
    definition = function(x, ..., delayed = TRUE, fill = NA,
                          use.names = TRUE) {
      objects <- list(x, ...)
      .checkSubsetNames(objects)
      out <- callNextMethod()
      .combineSubsets(out, objects, "cols", use.names)
    }
  )

  setMethod(
    f = "combineRows",
    signature = cls,
    definition = function(x, ..., delayed = TRUE, fill = NA,
                          use.names = TRUE) {
      objects <- list(x, ...)
      .checkSubsetNames(objects)
      out <- callNextMethod()
      .combineSubsets(out, objects, "rows", use.names)
    }
  )
}

# Methods for the classes that extend SingleCellExperiment -----------------

#' @name ExperimentSubset-altExp
#' @rdname ExperimentSubset-altExp
#' @eval .rdMethods(c("altExp", "altExps", "altExpNames", "altExpNames<-",
#'   "altExp<-", "altExps<-"), .esSCEClasses)
NULL

#' @name ExperimentSubset-reducedDim
#' @rdname ExperimentSubset-reducedDim
#' @eval .rdMethods(c("reducedDim", "reducedDim<-", "reducedDims",
#'   "reducedDims<-"), .esSCEClasses)
NULL

for (cls in .esSCEClasses) {
  setMethod(
    f = "altExp",
    signature = signature(x = cls, e = "character"),
    definition = function(x, e,  ...) {
      .altExp(x, e, ...)
    }
  )

  setMethod(
    f = "altExp",
    signature = signature(x = cls, e = "missing"),
    definition = function(x, e,  ...) {
      .altExp(x, e, ...)
    }
  )

  setMethod(
    f = "altExp",
    signature = signature(x = cls, e = "numeric"),
    definition = function(x, e,  ...) {
      .altExp(x, e, ...)
    }
  )

  setMethod(
    f = "altExps",
    signature = signature(x = cls),
    definition = function(x, ...) {
      .altExps(x, ...)
    }
  )

  setMethod(
    f = "altExpNames",
    signature = signature(x = cls),
    definition = function(x, ...) {
      .altExpNames(x, ...)
    }
  )

  setReplaceMethod(
    f = "altExpNames",
    signature = signature(x = cls, value = "character"),
    definition = function(x, ..., value) {
      .altExpNames(x, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "altExpNames",
    signature = signature(x = cls, value = "ANY"),
    definition = function(x, ..., value) {
      # every value that reaches this method is an error from .altExpNames<-
      `.altExpNames<-`(x, ..., value = value)
    }
  )

  setReplaceMethod(
    f = "altExp",
    signature = signature(x = cls, e = "character"),
    definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                          value) {
      .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
              ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "altExp",
    signature = signature(x = cls, e = "missing"),
    definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                          value) {
      .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
              ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "altExp",
    signature = signature(x = cls, e = "numeric"),
    definition = function(x, e, withDimnames = TRUE, withColData = FALSE, ...,
                          value) {
      .altExp(x, e, withDimnames = withDimnames, withColData = withColData,
              ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "altExps",
    signature = signature(x = cls),
    definition = function(x, withDimnames = TRUE, withColData = FALSE, ...,
                          value) {
      .altExps(x, withDimnames = withDimnames, withColData = withColData,
               ...) <- value
      return(x)
    }
  )

  setMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "character"),
    definition = function(x, type, ...) {
      .reducedDim(x, type, ...)
    }
  )

  setMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "missing"),
    definition = function(x, type, ...) {
      .reducedDim(x, type, ...)
    }
  )

  setMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "numeric"),
    definition = function(x, type, ...) {
      .reducedDim(x, type, ...)
    }
  )

  setReplaceMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "character"),
    definition = function(x, type, withDimnames = TRUE, ..., value) {
      .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "missing"),
    definition = function(x, type, withDimnames = TRUE, ..., value) {
      .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
      return(x)
    }
  )

  setReplaceMethod(
    f = "reducedDim",
    signature = signature(x = cls, type = "numeric"),
    definition = function(x, type, withDimnames = TRUE, ..., value) {
      .reducedDim(x, type, withDimnames = withDimnames, ...) <- value
      return(x)
    }
  )

  setMethod(
    f = "reducedDims",
    signature = signature(x = cls),
    definition = function(x, ...) {
      .reducedDims(x, ...)
    }
  )

  setReplaceMethod(
    f = "reducedDims",
    signature = signature(x = cls),
    definition = function(x, withDimnames = TRUE, ..., value) {
      .reducedDims(x, withDimnames = withDimnames, ...) <- value
      x
    }
  )
}

# the loop variable would otherwise stay in the namespace
rm(cls)

# Documentation helpers -----------------------------------------------------

# Returns roxygen tag lines documenting the methods of the generics named in
# `generics` for the classes in `classes`: an @aliases line and an @usage
# line with one \S4method entry per method, as roxygen writes them for a
# method defined outside a loop. roxygen runs @eval after load_all(), so the
# methods are registered by then. Methods are listed by class, then in the
# order of `generics`
.rdMethods <- function(generics, classes = .esClasses) {
  ns <- asNamespace("ExperimentSubset")
  methodLists <- list()
  for (generic in generics) {
    methodList <- list()
    if (methods::isGeneric(generic, where = ns)) {
      methodList <- methods::findMethods(generic, where = ns)
    }
    # names are the signatures, such as "SubsetSingleCellExperiment#ANY"
    missing <- setdiff(classes, sub("#.*", "", names(methodList)))
    if (length(missing) > 0) {
      # a misspelled generic, or a page listing classes its methods aren't
      # defined for, would otherwise drop docs without a message
      stop("'generics' includes '", generic, "', which has no methods ",
           "for ", paste(missing, collapse = ", "), ".")
    }
    methodLists[[generic]] <- methodList
  }
  aliases <- character(0)
  usages <- character(0)
  for (cls in classes) {
    for (generic in generics) {
      methodList <- methodLists[[generic]]
      keys <- names(methodList)
      keys <- keys[sub("#.*", "", keys) == cls]
      for (key in keys) {
        signature <- gsub("#", ",", key, fixed = TRUE)
        aliases <- c(aliases, paste0(generic, ",", signature, "-method"))
        usages <- c(usages,
                    .rdUsage(generic, signature, methodList[[key]]))
      }
    }
  }
  c(paste("@aliases", paste(aliases, collapse = " ")),
    paste("@usage", paste(usages, collapse = "\n\n")))
}

# One \S4method usage entry, such as
# "\S4method{colData}{SubsetSummarizedExperiment,DataFrame}(x, ...) <- value"
.rdUsage <- function(generic, signature, method) {
  # a method with arguments the generic lacks, such as combineCols()'s
  # delayed and fill, is stored wrapped in .local(); use the definition
  args <- formals(methods::unRematchDefinition(method))
  isReplacement <- endsWith(generic, "<-")
  if (isReplacement) {
    args <- args[names(args) != "value"]
  }
  # an argument without a default is the empty symbol
  argText <- vapply(seq_along(args), function(i) {
    if (identical(args[[i]], quote(expr = ))) {
      return(names(args)[i])
    }
    default <- paste(deparse(args[[i]]), collapse = "")
    # \, %, { and } are special in Rd, as in a default such as "%s"
    paste(names(args)[i], "=", gsub("([\\\\%{}])", "\\\\\\1", default))
  }, character(1))
  usage <- sprintf("\\S4method{%s}{%s}(%s)", sub("<-$", "", generic),
                   signature, paste(argText, collapse = ", "))
  if (isReplacement) {
    usage <- paste(usage, "<- value")
  }
  usage
}
