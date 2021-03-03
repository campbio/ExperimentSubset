#' @title Subset creation method for ExperimentSubset objects
#' @description Create a subset from an already available \code{assay} in the
#'   input \code{ExperimentSubset} object by specifying the rows and columns to
#'   include in the subset.
#' @param x \code{ExperimentSubset} Specify the object from which a subset
#'   should be created. Input can also be any object inherited from
#'   \code{SummarizedExperiment} for immediate conversion and subset formation.
#'   A list of slots can also be passed to directly construct an ES object from
#'   matrices similar to SE and SCE constructors.
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

#' @title Name retrieval method for all subset assays in ExperimentSubset
#'   objects
#' @description Retrieves the names of all the subsets as well as the subset
#'   assays.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @return A \code{vector} containing the names of the subsets and the subset
#'   assays available in the input \code{ExperimentSubset} object.
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
setGeneric(
  name = "subsetAssayNames",
  def = function(x)
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
#' @param inputMatrix \code{dgCMatrix} The input subset \code{assay}.
#' @param subsetAssayName \code{character(1)} Specify the name of the new
#'   \code{assay} against the \code{inputMatrix} parameter.
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
                 subsetAssayName)
  {
    standardGeneric("setSubsetAssay")
  }
)

#' @title Subset parent hierarchy retrieval method for ExperimentSubset objects
#' @description Retrieves a complete 'subset to parent' link from a specified
#'   subset.
#' @param x \code{ExperimentSubset} Input \code{ExperimentSubset} object.
#' @param subsetName \code{character(1)} Specify the name of the subset against
#'   which the 'subset to parent link' should be retrieved.
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
#' subsetParent(es, "subset1pAssay")
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
#' @param ... Additional arguments to pass to into the SCE method.
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
#'   rowData should be combined or not. Default \code{FALSE}.
#' @return The \code{rowData} from input object.
#' @rdname subsetRowData
#' @importMethodsFrom SingleCellExperiment rowData
#' @export
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
#'   colData should be combined or not. Default \code{FALSE}.
#' @return The \code{colData} from input object.
#' @rdname subsetColData
#' @importMethodsFrom SingleCellExperiment colData
#' @export
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
#' @param value Input \code{DataFrame} to store.
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

#' @title subsetColnames
#' @description Get \code{colnames} from a subset in the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object 
#' @param subsetName Name of the subset to get \code{colnames} from.
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
#' @param subsetName Name of the subset to get \code{colnames} from.
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

#' @title subsetRownames
#' @description Set \code{colnames} to a subset in the \code{ExperimentSubset} object.
#' @param x Input \code{ExperimentSubset} object.
#' @param subsetName Name of the subset to set \code{colnames} to.
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
#' @param subsetName Specify 'subset name' or 'subset assay name' to fetch the
#'   assay from.
#' @return Subset assay
#' @importMethodsFrom SummarizedExperiment assay
#' @export
setGeneric(
  name = "getSubsetAssay",
  def = function(x,
                 subsetName)
  {
    standardGeneric("getSubsetAssay")
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
setGeneric(
  name = "subsetColLinks",
  def = function(x, subsetName, parentColLinkData)
  {
    standardGeneric("subsetColLinks")
  }
)