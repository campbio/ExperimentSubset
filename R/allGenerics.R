#' @title Subset creation method for ExperimentSubset objects
#' @description Create a subset from an already available \code{assay} in the
#'   input \code{ExperimentSubset} object by specifying the rows and columns to
#'   include in the subset.
#' @param x \code{ExperimentSubset} Specify the object from which a subset
#'   should be created. Use \code{ExperimentSubset()} first to convert another
#'   experiment object.
#' @param subsetName \code{character(1)} Specify the name of the subset to
#'   create.
#' @param rows \code{numeric} or \code{character} vector. Specify the rows to
#'   include in this subset. If \code{missing} or \code{NULL}, all rows are
#'   included in the subset. Default \code{NULL}.
#' @param cols \code{numeric} or \code{character} vector. Specify the columns
#'   to include in this subset. If \code{missing} or \code{NULL}, all columns
#'   are included in the subset. Default \code{NULL}.
#' @param parentAssay \code{character(1)} Specify the parent \code{assay} of the
#'   subset: the name of a main \code{assay}, of a subset, or of an assay seen
#'   in a subset written as \code{"subset/assay"}. If \code{NULL}, the first
#'   available main \code{assay} will be marked as parent. Default
#'   \code{NULL}.
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

#' @title Name retrieval method for all subset assays in ExperimentSubset
#'   objects
#' @description Retrieves the names of all the subsets as well as the subset
#'   assays.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} If given, list the assays seen in
#'   this subset instead: its own assays and those linked from its parent.
#'   Default \code{NULL}.
#' @return A \code{vector} containing the names of the subsets followed by
#'   each assay stored in a subset as \code{"subset/assay"}, or, with
#'   \code{subsetName}, the names of the assays seen in that subset.
#' @rdname subsetAssayNames
#' @importMethodsFrom SummarizedExperiment assayNames
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
#' subsetAssayNames(es, subsetName = "subset1")
setGeneric(
  name = "subsetAssayNames",
  def = function(x, subsetName = NULL)
  {
    standardGeneric("subsetAssayNames")
  }
)

#' @title Get names of only the subsets in ExperimentSubset objects
#' @description Retrieves the names of the available subsets (not the subset 
#'   assays) in an \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Specify the input ES object.
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

#' @title Method for storing new assays inside subsets in ExperimentSubset 
#'   objects
#' @description Store a new subset \code{assay} inside a specified subset in the
#'   input \code{ExperimentSubset} object.
#' @param x \code{ExperimentSubset} Specify the input object.
#' @param subsetName \code{character(1)} Specify the name of the existing subset
#'   inside which the new subset \code{assay} should be stored.
#' @param inputMatrix A matrix-like object (for example a \code{matrix} or
#'   \code{dgCMatrix}) to store as an \code{assay} in the subset.
#' @param subsetAssayName \code{character(1)} Specify the name of the new
#'   \code{assay} against the \code{inputMatrix} parameter. It can't contain
#'   \code{"/"}. A subset sees its parent's assays as linked assays, so a name
#'   linked from the parent (for example a main \code{assay}) needs
#'   \code{unlink = TRUE}.
#' @param unlink \code{logical(1)} Store the \code{assay} even though the
#'   subset sees an assay of that name linked from its parent. The subset, and
#'   subsets made from it, then use this \code{assay} instead; the parent's is
#'   unchanged. Default \code{FALSE}.
#' @return Updated \code{ExperimentSubset} object with the new \code{assay}
#'   stored inside the specified subset.
#' @rdname setSubsetAssay
#' @importMethodsFrom SummarizedExperiment assay<-
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
#' es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
#' es
setGeneric(
  name = "setSubsetAssay",
  def = function(x,
                 subsetName,
                 inputMatrix,
                 subsetAssayName,
                 unlink = FALSE)
  {
    standardGeneric("setSubsetAssay")
  }
)

#' @title Subset parent hierarchy retrieval method for ExperimentSubset objects
#' @description Retrieves a complete 'subset to parent' link from a specified
#'   subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Specify the name of the subset, or an
#'   assay stored in a subset as \code{"subset/assay"}, against which the
#'   'subset to parent link' should be retrieved.
#' @return A \code{list} containing the 'subset to parent' link.
#' @rdname subsetParent
#' @importMethodsFrom SummarizedExperiment assayNames
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
#' subsetParent(es, "subset1/subset1pAssay")
setGeneric(
  name = "subsetParent",
  def = function(x, subsetName)
  {
    standardGeneric("subsetParent")
  }
)

#' @title Method for displaying 'child-parent' link structure of subsets in
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
#' @importMethodsFrom SingleCellExperiment altExpNames
#' @importMethodsFrom SummarizedExperiment assayNames
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

#' @title reducedDimNames
#' @description A wrapper to the \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} method.
#' @param ... Additional arguments. Use \code{subsetName} to get the
#'   \code{reducedDimNames} of a subset; other arguments are passed to the
#'   SCE method.
#' @return The \code{reducedDimNames} from the specified subset or same as \code{reducedDimNames} from \link[SingleCellExperiment]{reducedDims} when \code{subsetName} is \code{missing}.
#' @rdname reducedDimNames
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:50), cols = c(1:20),
#' parentAssay = "counts")
#' # two small reduced dimensions, one row per cell in the subset
#' set.seed(1)
#' pca1 <- matrix(rnorm(40), nrow = 20)
#' pca2 <- matrix(rnorm(40), nrow = 20)
#' reducedDims(es, subsetName = "subset1") <- list(PCA_1 = pca1, PCA_2 = pca2)
#' reducedDimNames(es, subsetName = "subset1")
setGeneric(
  name = "reducedDimNames",
  def = function(x, ...)
  {
    standardGeneric("reducedDimNames")
  }
)

#' @title Subset count method for ExperimentSubset objects
#' @description Get the total count of the available subsets (excluding subset
#'   assays) in an \code{ExperimentSubset} object.
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

#' @title Accessor method for rowData from subsets in ExperimentSubset objects
#' @description Get \code{rowData} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{rowData} from.
#' @param parentRowData \code{logical(1)} Logical value indicating if parent
#'   rowData should be combined or not. A subset's own column hides the
#'   main object's column of the same name, with a message if their values
#'   differ. Default \code{FALSE}.
#' @return The \code{rowData} from input object.
#' @rdname subsetRowData
#' @importMethodsFrom SingleCellExperiment rowData
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetRowData(es, "subset1")
setGeneric(
  name = "subsetRowData",
  def = function(x, subsetName, parentRowData)
  {
    standardGeneric("subsetRowData")
  }
)

#' @title Accessor method for colData from subsets in ExperimentSubset objects
#' @description Get \code{colData} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{colData} from.
#' @param parentColData \code{logical(1)} Logical value indicating if parent
#'   colData should be combined or not. A subset's own column hides the
#'   main object's column of the same name, with a message if their values
#'   differ. Default \code{FALSE}.
#' @return The \code{colData} from input object.
#' @rdname subsetColData
#' @importMethodsFrom SingleCellExperiment colData
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetColData(es, "subset1")
setGeneric(
  name = "subsetColData",
  def = function(x, subsetName, parentColData)
  {
    standardGeneric("subsetColData")
  }
)

#' @title Setter method for colData to subsets in ExperimentSubset objects
#' @description Set \code{colData} to a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to set
#'   \code{colData} to.
#' @param value Input \code{DataFrame} to store. Its rownames must be
#'   \code{NULL} or the subset's column names, in the same order; other
#'   rownames are an error. Use \code{subsetColnames<-} to rename columns.
#' @return Input object with \code{colData} stored.
#' @rdname subsetColData
#' @importMethodsFrom SingleCellExperiment colData
#' @export
setGeneric(
  name = "subsetColData<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetColData<-")
  }
)

#' @title Setter method for rowData to subsets in ExperimentSubset objects
#' @description Set \code{rowData} to a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to set
#'   \code{rowData} to.
#' @param value Input \code{DataFrame} to store.
#' @return Input object with \code{rowData} stored.
#' @rdname subsetRowData
#' @importMethodsFrom SingleCellExperiment rowData
#' @export
setGeneric(
  name = "subsetRowData<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetRowData<-")
  }
)

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
#' rows = c(1:50), cols = c(1:20),
#' parentAssay = "counts")
#' # two small reduced dimensions, one row per cell in the subset
#' set.seed(1)
#' pca1 <- matrix(rnorm(40), nrow = 20)
#' pca2 <- matrix(rnorm(40), nrow = 20)
#' reducedDims(es, subsetName = "subset1") <- list(PCA_1 = pca1, PCA_2 = pca2)
#' reducedDims(es, subsetName = "subset1")
setGeneric(
  name = "reducedDims<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("reducedDims<-")
  }
)

#' @title subsetColnames
#' @description Get \code{colnames} from a subset in the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object 
#' @param subsetName Name of the subset to get \code{colnames} from.
#' @return A \code{vector} of \code{colnames}.
#' @rdname subsetColnames
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetColnames(es, "subset1")
setGeneric(
  name = "subsetColnames",
  def = function(x, subsetName)
  {
    standardGeneric("subsetColnames")
  }
)

#' @title subsetColnames
#' @description Set \code{colnames} to a subset in the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to set \code{colnames} to.
#' @param value Specify the colname values to replace.
#' @return Input object with colnames set to a subset.
#' @rdname subsetColnames
#' @export
setGeneric(
  name = "subsetColnames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetColnames<-")
  }
)

#' @title subsetRownames
#' @description Get \code{rownames} from a subset in the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to get \code{rownames} from or set
#'   them to.
#' @return A \code{vector} of \code{rownames}.
#' @rdname subsetRownames
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' subsetRownames(es, "subset1")
setGeneric(
  name = "subsetRownames",
  def = function(x, subsetName)
  {
    standardGeneric("subsetRownames")
  }
)

#' @title subsetRownames
#' @description Set \code{rownames} to a subset in the
#'   \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to get \code{rownames} from or set
#'   them to.
#' @param value Specify the rownames values to replace.
#' @return Input object with rownames set to a subset.
#' @rdname subsetRownames
#' @export
setGeneric(
  name = "subsetRownames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetRownames<-")
  }
)

#' Get subset assay from an \code{ExperimentSubset} object.
#'
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName The name of the subset to fetch the assay from, or an
#'   assay seen in a subset written as \code{"subset/assay"}.
#' @param subsetAssayName \code{character(1)} The name of an assay seen in the
#'   subset \code{subsetName}: one stored with \code{setSubsetAssay}, or one
#'   linked from its parent. Default \code{NULL}, which returns the subset's
#'   slice of its \code{parentAssay}.
#' @return Subset assay
#' @importMethodsFrom SummarizedExperiment assay
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es,
#' "subset1",
#' rows = c(10,11,50,56,98,99,102,105,109, 200),
#' cols = c(20,21,40,45,90,99,100,123,166,299),
#' parentAssay = "counts")
#' getSubsetAssay(es, "subset1")
#' # store a changed copy of the subset's counts inside the subset
#' subset1Counts <- getSubsetAssay(es, "subset1")
#' es <- setSubsetAssay(es, "subset1", subset1Counts + 1, "plusOne")
#' getSubsetAssay(es, "subset1", subsetAssayName = "plusOne")
setGeneric(
  name = "getSubsetAssay",
  def = function(x,
                 subsetName,
                 subsetAssayName = NULL)
  {
    standardGeneric("getSubsetAssay")
  }
)

#' @title Take a subset out as a standalone experiment object
#' @description Returns a subset as an object of the main object's class
#'   (for example a \code{SingleCellExperiment}), so that any tool written
#'   for that class can be run on the subset. Store what the tool added with
#'   \code{setSubsetExperiment}.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} The name of the subset.
#' @param parentData \code{logical(1)} Include the main object's rowData and
#'   colData columns for the subset's rows and columns, as with
#'   \code{parentRowData} and \code{parentColData}. They are marked so that
#'   \code{setSubsetExperiment} doesn't store copies of them. A subset's own
#'   column hides the main object's column of the same name, with a message
#'   if their values differ. Default \code{TRUE}.
#' @return An object of the main object's class with the subset's rows and
#'   columns and row and column names. Its assays are every assay the subset
#'   sees, its own and those linked from its parent, under their names, with
#'   the subset's primary assay first. Its rowData and colData are the
#'   subset's, after the main object's columns if \code{parentData = TRUE}.
#'   Its reducedDims, altExps and metadata are the subset's. Everything else,
#'   such as row ranges, row and column trees, spatial coordinates, images,
#'   \code{colPairs} and \code{rowPairs}, comes from the main object.
#' @rdname getSubsetExperiment
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1", rows = 1:100, cols = 1:50)
#' subset1 <- getSubsetExperiment(es, "subset1")
#' subset1
setGeneric(
  name = "getSubsetExperiment",
  def = function(x, subsetName, parentData = TRUE)
  {
    standardGeneric("getSubsetExperiment")
  }
)

#' @title Store a standalone experiment object back into a subset
#' @description Stores what a tool added to an object returned by
#'   \code{getSubsetExperiment}: new assays, rowData and colData columns,
#'   reducedDims, altExps and metadata.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} The name of the subset.
#' @param value An object of the main object's class (not a subclass or a
#'   parent class) with the subset's dimensions and the subset's row and
#'   column names, in the same order. A tool that drops rows or columns is
#'   followed by \code{createSubset} instead.
#' @param unlink \code{logical(1)} Store assays of \code{value} that have
#'   the name of an assay linked from the subset's parent but different
#'   values, as with \code{setSubsetAssay}. Without it, such an assay is an
#'   error. Default \code{FALSE}.
#' @return Updated \code{ExperimentSubset} object. Linked assays that are
#'   in \code{value} hold the same values as what the subset sees (whatever
#'   their dimnames, type or class) are not stored, and neither are the main
#'   object's rowData and colData columns added by \code{getSubsetExperiment}.
#'   Everything else in \code{value} replaces the subset's own data: an own
#'   assay, reducedDim, altExp or column that \code{value} no longer has is
#'   removed from the subset, except an assay another subset was made from,
#'   which is an error. Nothing else is stored: not row ranges, row and
#'   column trees, spatial coordinates, images, \code{colPairs},
#'   \code{rowPairs} or \code{mainExpName}.
#' @rdname setSubsetExperiment
#' @export
#' @examples
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1", rows = 1:100, cols = 1:50)
#' subset1 <- getSubsetExperiment(es, "subset1")
#' # add a reducedDim and a colData column to the standalone object
#' reducedDim(subset1, "random") <- matrix(rnorm(100), nrow = 50)
#' subset1$group <- rep(c("a", "b"), 25)
#' # and store them in the subset
#' es <- setSubsetExperiment(es, "subset1", subset1)
#' reducedDimNames(es, subsetName = "subset1")
setGeneric(
  name = "setSubsetExperiment",
  def = function(x, subsetName, value, unlink = FALSE)
  {
    standardGeneric("setSubsetExperiment")
  }
)

#' @title reducedDimNames<-
#' @description A wrapper to the \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method with additional support for subsets.
#' @param x Input \code{ExperimentSubset} object or any object supported by \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @param subsetName Specify the name of the subset to which the \code{reducedDimNames<-} should be set to. If \code{missing}, \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method is called on the main object.
#' @param value Input value same as \code{reducedDimNames<-} from \link[SingleCellExperiment]{reducedDims} method.
#' @return Input object with \code{reducedDimNames<-} set.
#' @rdname reducedDimNames-set
#' @export
#' @examples 
#' data(sce_chcl, package = "scds")
#' es <- ExperimentSubset(sce_chcl)
#' es <- createSubset(es, "subset1",
#' rows = c(1:50), cols = c(1:20),
#' parentAssay = "counts")
#' # two small reduced dimensions, one row per cell in the subset
#' set.seed(1)
#' pca1 <- matrix(rnorm(40), nrow = 20)
#' pca2 <- matrix(rnorm(40), nrow = 20)
#' reducedDims(es, subsetName = "subset1") <- list(PCA_1 = pca1, PCA_2 = pca2)
#' reducedDimNames(es, subsetName = "subset1") <- c("rDim1", "rDim2")
setGeneric(
  name = "reducedDimNames<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("reducedDimNames<-")
  }
)

#' @title Accessor method for rowLinks from subsets in ExperimentSubset objects
#' @description Get \code{rowLinks} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{rowLinks} from.
#' @param parentRowLinkData \code{logical(1)} Logical value indicating if parent
#'   rowLinks should be combined or not. Default \code{FALSE}.
#' @return The \code{rowLinks} from input object.
#' @rdname subsetRowLinks
#' @importMethodsFrom TreeSummarizedExperiment rowLinks
#' @export
#' @examples
#' # a small TreeSummarizedExperiment with a tree on both rows and columns
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' counts <- matrix(1:100, nrow = 10)
#' # name rows and columns after the tips of the tree
#' rownames(counts) <- tinyTree$tip.label
#' colnames(counts) <- tinyTree$tip.label
#' tse <- TreeSummarizedExperiment(assays = list(counts = counts),
#'                                 rowTree = tinyTree, colTree = tinyTree)
#' es <- ExperimentSubset(tse)
#' es <- createSubset(es, "subset1", rows = 1:5, cols = 1:5)
#' subsetRowLinks(es, "subset1", parentRowLinkData = TRUE)
setGeneric(
  name = "subsetRowLinks",
  def = function(x, subsetName, parentRowLinkData)
  {
    standardGeneric("subsetRowLinks")
  }
)

#' @title Accessor method for colLinks from subsets in ExperimentSubset objects
#' @description Get \code{colLinks} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{colLinks} from.
#' @param parentColLinkData \code{logical(1)} Logical value indicating if parent
#'   colLinks should be combined or not. Default \code{FALSE}.
#' @return The \code{colLinks} from input object.
#' @rdname subsetColLinks
#' @importMethodsFrom TreeSummarizedExperiment colLinks
#' @export
#' @examples
#' # a small TreeSummarizedExperiment with a tree on both rows and columns
#' data("tinyTree", package = "TreeSummarizedExperiment")
#' counts <- matrix(1:100, nrow = 10)
#' # name rows and columns after the tips of the tree
#' rownames(counts) <- tinyTree$tip.label
#' colnames(counts) <- tinyTree$tip.label
#' tse <- TreeSummarizedExperiment(assays = list(counts = counts),
#'                                 rowTree = tinyTree, colTree = tinyTree)
#' es <- ExperimentSubset(tse)
#' es <- createSubset(es, "subset1", rows = 1:5, cols = 1:5)
#' subsetColLinks(es, "subset1", parentColLinkData = TRUE)
setGeneric(
  name = "subsetColLinks",
  def = function(x, subsetName, parentColLinkData)
  {
    standardGeneric("subsetColLinks")
  }
)

#' @title Accessor method for spatialCoords from subsets in ExperimentSubset objects
#' @description Get \code{spatialCoords} from a subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{spatialCoords} from.
#' @return The \code{spatialCoords} from input object.
#' @rdname subsetSpatialCoords
#' @importMethodsFrom SpatialExperiment spatialCoords
#' @export
#' @examples
#' # a small SpatialExperiment with six spots
#' counts <- matrix(1:60, nrow = 10)
#' coords <- cbind(x = 1:6, y = 6:1)
#' spe <- SpatialExperiment(assays = list(counts = counts),
#'                          spatialCoords = coords)
#' es <- ExperimentSubset(spe)
#' es <- createSubset(es, "subset1", cols = 2:4)
#' subsetSpatialCoords(es, "subset1")
setGeneric(
  name = "subsetSpatialCoords",
  def = function(x, subsetName)
  {
    standardGeneric("subsetSpatialCoords")
  }
)

#' @title Deprecated: spatialData of subsets in ExperimentSubset objects
#' @description \code{subsetSpatialData()} and
#'   \code{spatialData(x, subsetName =)} are deprecated, and
#'   \code{subsetSpatialData<-} is defunct, because
#'   \code{spatialData()} and \code{spatialData<-} are deprecated in
#'   SpatialExperiment, which now keeps these columns in \code{colData}.
#'   Use \code{subsetColData(x, subsetName, parentColData = TRUE)} to get
#'   them, and \code{subsetColData(x, subsetName) <- value} to store a
#'   subset's own columns. The getters will be defunct in the next
#'   release.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to get
#'   \code{spatialData} from.
#' @param parentSpatialData \code{logical(1)} Logical value indicating if parent
#'   spatialData should be combined or not. Default \code{FALSE}.
#' @return The \code{spatialData} from input object, with a deprecation
#'   warning. The setter fails with an error.
#' @seealso \code{\link{subsetColData}}
#' @rdname subsetSpatialData
#' @importMethodsFrom SpatialExperiment spatialData
#' @export
setGeneric(
  name = "subsetSpatialData",
  def = function(x, subsetName, parentSpatialData)
  {
    standardGeneric("subsetSpatialData")
  }
)

#' @title Defunct: set spatialData of subsets in ExperimentSubset objects
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Name of the subset to set
#'   \code{spatialData} to.
#' @param value Ignored: the setter is defunct.
#' @rdname subsetSpatialData
#' @importMethodsFrom SpatialExperiment spatialData
#' @export
setGeneric(
  name = "subsetSpatialData<-",
  def = function(x, subsetName, value)
  {
    standardGeneric("subsetSpatialData<-")
  }
)
