# ExperimentSubset 1.23.25 (2026-11-01)
* `cbind()` and `rbind()` now keep the subsets of every object, with their own assays, rowData, colData, reducedDims and altExps. Before, only the first object's subsets were kept and the others were dropped without a message. Subsets of later objects move to where that object's columns or rows are in the result. Subset names repeated across the objects are an error, as is `cbind()` of an object with subsets whose rows aren't the first object's, in the same order. For SpatialExperiment objects, each subset's `sample_id` follows the main object's when `cbind()` renames repeated sample IDs. Convert a plain object with `ExperimentSubset()` before combining it with `ExperimentSubset` objects (#67)

# ExperimentSubset 1.23.24 (2026-10-31)
* `subsetSpatialData()` and `spatialData(x, subsetName = )` are deprecated, as `spatialData()` is in SpatialExperiment, which keeps these columns in colData. They warn and still work; use `subsetColData(x, subsetName)` instead, which returns all of the subset's colData columns, with `parentColData = TRUE` to include the main object's. They will be defunct in the next release. `spatialData(x)` on the main object is unchanged (#34)
* `subsetSpatialData<-` is defunct and fails with an error saying to add each column with `colData(x, subsetName = "name")$column <- values`. It skips the usual deprecation release because it already failed for every value: SpatialExperiment's `spatialData<-` takes no `data.frame`, and the method took no `DataFrame` (#34)

# ExperimentSubset 1.23.23 (2026-10-31)
* `createSubset()` from a parent with no rows or no columns, such as a subset none of whose rows remain after `x[i, j]`, now says that the parent is empty ("Parent 's1' has no rows, so a subset can't be made from it."). Before, the error blamed the `rows` or `cols` argument, even when it was left at its default (#68)

# ExperimentSubset 1.23.22 (2026-10-31)
* `x[i, j]` now updates the subsets of `x`: each keeps the rows and columns of the main object that remain, in its own order, with its own assays, rowData, colData, reducedDims and altExps, and subsets of subsets follow what remains of their parent. A subset none of whose rows or columns remain is kept with no rows or no columns. Before, subsets kept their old row and column positions, so reading one failed with "subscript out of bounds" or, when the positions were still in range, silently returned the wrong rows or columns under the subset's names (#4)

# ExperimentSubset 1.23.21 (2026-10-31)
* `rowData(x, subsetName = , parentRowData = TRUE)`, `colData(x, subsetName = , parentColData = TRUE)`, `subsetRowData()`, `subsetColData()` and `getSubsetExperiment(parentData = TRUE)` no longer return two columns with the same name when a subset has its own column with the name of a main-object column, of which `$` and `[[` returned the main object's. The subset's column now hides the main object's, and comes after the main object's columns, with a message if their values differ. A SpatialExperiment subset has one `sample_id` column, the main object's current one (#57)

# ExperimentSubset 1.23.20 (2026-10-31)
* `stats` is no longer listed in Suggests: it is a base package that R attaches by default, including in `R CMD check`

# ExperimentSubset 1.23.19 (2026-10-31)
* Joshua D. Campbell is now the maintainer, taking over from Irzam Sarfraz, who remains an author

# ExperimentSubset 1.23.18 (2026-10-25)
* New article on the website: an end-to-end single-cell workflow on 4,000 PBMCs that keeps empty droplets, cells, quality-controlled cells, variable genes and a T-cell re-analysis as subsets of one object, using both `getSubsetExperiment()`/`setSubsetExperiment()` and direct storage. The vignette links to the website and the article, shows `getSubsetExperiment()`/`setSubsetExperiment()`, and sets subset colData columns directly with `colData(x, subsetName = )$column <-`

# ExperimentSubset 1.23.17 (2026-10-25)
* The package has a pkgdown website at <https://campbio.github.io/ExperimentSubset/>, with the release version at the root and the devel version under `/dev/`

# ExperimentSubset 1.23.16 (2026-10-25)
* New `getSubsetExperiment()` returns a subset as a standalone object of the main object's class (for example a `SingleCellExperiment`), with every assay the subset sees, its rowData and colData (after the main object's columns, with `parentData = TRUE`), reducedDims, altExps and metadata, so tools written for that class can run on the subset. `setSubsetExperiment()` stores the result back: its assays, rowData, colData, reducedDims, altExps and metadata replace the subset's own, so an own assay or column missing from it is removed, while unchanged linked assays and the main object's columns are not stored. A linked assay the tool changed needs `unlink = TRUE`, and removing an assay another subset was made from is an error

# ExperimentSubset 1.23.15 (2026-10-25)
* A subset now sees all of its parent's assays as linked assays: the parent's values sliced to the subset's rows and columns, under the same names, computed when asked for rather than copied. For example a `logcounts` computed on the main object is available in every subset. Read an assay as a subset sees it with `assay(x, "subset", subsetAssayName = "a")`, `getSubsetAssay(x, "subset", subsetAssayName = "a")` or the name `"subset/a"`; `subsetAssayNames(x, subsetName = "subset")` lists them (#2)
* Each subset has its own assay names: two subsets made from the same parent can both store, say, `decontXcounts`. Storing an assay under a name the subset already sees through its parent now needs `unlink = TRUE` in `setSubsetAssay()` or `assay<-`; the subset and subsets made from it then use their own version (including `assay(x, "subset")` when the name is the one its `parentAssay` refers to), and the parent's is unchanged
* A bare name in `assay()`, `getSubsetAssay()`, `createSubset(parentAssay = )`, `subsetParent()`, `subsetRownames()` and `subsetColnames()` now means a main assay or a subset. An assay stored in a subset is named `"subset/assay"`; a bare stored-assay name is an error that shows the qualified form. Subsets made from such an assay record the qualified name, and subsets in objects saved by earlier versions that recorded a bare name are switched to the qualified name the next time a subset is created or an assay stored. `assay(x, "subset/assay") <- value` stores into that subset
* `subsetAssayNames()` lists the assays stored in subsets as `"subset/assay"`, so each listed name can be passed back to `assay()`. Subset names and stored assay names can't contain `/`

# ExperimentSubset 1.23.14 (2026-10-22)
* `withDimnames` (for `assay()`, `assay<-`, `reducedDim<-`, `altExp<-` and `altExps<-`), `withColData` (for `altExp<-` and `altExps<-`) and `use.names` (for `rowData()`) now take effect, on the main object and with `subsetName` (where the `reducedDim()`, `reducedDims()`, `altExp()` and `altExps()` getters also pass them on); they were silently dropped. `reducedDims<-` is unchanged (#26)
* `rowData(x, subsetName = ) <-` and `colData(x, subsetName = ) <-` with `NULL` or a `data.frame` now change the subset; before, they changed the main object's rowData or colData
* `reducedDim(x, subsetName = )` without a type now returns, or with `<-` replaces, the subset's first reducedDim instead of failing
* Writing back rowData or colData read with `parentRowData = TRUE` or `parentColData = TRUE` no longer stores copies of the main object's columns in the subset: those columns are left out, with a message naming them. The getters mark them in the column metadata (`mcols()`) of the result
* `subsetRowLinks()`, `subsetColLinks()` and `subsetSpatialData()` return the subset's row or column names instead of `"1".."n"`, and asking for the parent's links no longer fails when the main object has no tree for that dimension
* `subsetParent()` resolves stored assay names the same way as `assay()` and stops with an error if the chain of parents has a cycle
* Loading the package no longer warns about the `listData` and `elementType` slots in the class prototypes

# ExperimentSubset 1.23.13 (2026-10-22)
* Documentation fixes: the vignette lists all five supported classes, names the `createSubset()` arguments `rows` and `cols`, seeds its random data, drops the pin to Bioconductor 3.11 and installs from GitHub with BiocManager; help pages say which input `ExperimentSubset()` and `createSubset()` accept, fix the `subsetRownames()` page, and add examples for `subsetRowLinks()`, `subsetColLinks()` and `subsetSpatialCoords()`; examples no longer need scater
* DESCRIPTION gains URL and BugReports fields and drops the stale Date and the unused LazyData; the citation uses `bibentry()`

# ExperimentSubset 1.23.12 (2026-10-22)
* `rowData()`, `colData()`, `subsetRowData()` and `subsetColData()` for a subset now return the subset's row or column names as rownames (or `NULL` if it has none), instead of `"1".."n"`, so writing the result back no longer renames the subset's cells
* Setting a subset's colData no longer changes its column names: a value without rownames takes the subset's names, and a value whose rownames differ from them is an error. Rename a subset's columns with `subsetColnames<-`. Subsets whose cells were already renamed to `"1".."n"` this way keep those names; set them back with `subsetColnames<-` or re-create the subset

# ExperimentSubset 1.23.11 (2026-10-22)
* The tests no longer download PBMC data from ExperimentHub: a small simulated data set replaces it, the tests skip, rather than fail, when a Suggests package they need is missing, and new tests cover the TreeSummarizedExperiment and SpatialExperiment methods (except `spatialData`, deprecated upstream)

# ExperimentSubset 1.23.10 (2026-10-21)
* `createSubset()` checks its input more strictly: `rows` and `cols` must select at least one row or column and be positive whole numbers or names, without duplicates (duplicates were accepted in some cases before); names given for an assay without row or column names give a clear error; `subsetName` must be a single non-empty string that is not already the name of an assay or of an assay stored in a subset; `parentAssay` must be a single string
* `ExperimentSubset()` now stops when `subset` contains NA instead of silently skipping the subset, and rejects a list of assays without names
* Unknown subset names now give a clear error in `subsetDim()`, `subsetColnames()`, `subsetRownames()` and their setters, and `reducedDims<-`; `assay(x, i)` and `getSubsetAssay(x, i)` with an unknown name are an error instead of returning `NULL`
* `assay(x, i) <- value` stores the assay inside the subset whenever `i` is a subset name, and then needs `subsetAssayName` (giving `subsetAssayName` with a name that is not a subset is an error); before, the choice depended on the dimensions of `value`, so a subset with the same dimensions as the main object got a new main assay named after it

# ExperimentSubset 1.23.9 (2026-10-21)
* A subset with a single row or column now returns a one-row or one-column matrix from `assay()`, not a vector, so subsets can be made from it
* `setSubsetAssay()` (and `assay(x, i) <- value` for a subset) no longer replaces the subset's row and column names with the input matrix's. A matrix without names is stored under the subset's names; a matrix whose names differ from the subset's, or are in another order, is now an error instead of silently relabelling the subset; a matrix of the wrong size gets an error that says so. `assay(x, subsetName)` now returns the subset's own row and column names, including any set with `subsetRownames<-` or `subsetColnames<-`, so `setSubsetAssay(x, s, assay(x, s), ...)` always works
* A `subsetName` with spaces is now stored under the name with spaces removed, which the warning already announced; before, the subset was stored under the original name and lookups of its stored assays failed. Refer to such a subset by the name without spaces

# ExperimentSubset 1.23.8 (2026-10-21)
* Subsets created from another subset (`parentAssay` set to a subset, or to an assay stored in one) now take their row and column names, `sample_id`, and the parent rowData, colData, rowLinks, colLinks, spatialCoords and spatialData from the right rows and columns of the main object. Before, the indices relative to the parent subset were applied to the main object, so these came from the wrong features and cells while `assay()` returned the right values. The parent metadata is looked up on access, so it is now right for existing objects too; the subset's row and column names (and `sample_id`) are stored when the subset is created, so objects saved with an earlier version keep the wrong names for such subsets. Repair them with `subsetRownames(x, "s2") <- rownames(assay(x, "s2"))` and the same for `subsetColnames()`; a wrong SpatialExperiment `sample_id` needs the subset re-created

# ExperimentSubset 1.23.7 (2026-10-21)
* Accessors and replacement methods called on the main object without `subsetName` now pass on arguments that travel in the base generic's `...`; for example `reducedDim(x, "PCA", withDimnames = FALSE)`, `reducedDims(x, withDimnames = FALSE)` and `altExps(x, withDimnames = FALSE)` failed with "argument x is missing". Arguments that are named formals of the base generic, such as `withDimnames` in `assay<-` or `use.names` in `rowData`, are still dropped (#29)

# ExperimentSubset 1.23.6 (2026-10-21)
* `altExp(x, e, subsetName = )` now returns the altExp named or numbered by `e`; it always returned the first one before
* `altExp(x, <number>) <- value` now replaces the altExp at that position; with `subsetName` it replaced the first altExp, and without it the call failed

# ExperimentSubset 1.23.5 (2026-10-21)
* The vignette's PBMC workflow now uses `scrapper` for QC, normalization, highly variable genes and PCA, replacing `scran` and `scater` functions that are deprecated
* `scrapper` replaces `scran` in Suggests; the tests also use `scrapper` in place of the deprecated functions

# ExperimentSubset 1.23.1 (2026-10-20)
* Fixed a vignette example that filled a 2 x 4 matrix with 100 values, which fails `R CMD check --as-cran`

# ExperimentSubset 1.1.0 (2021-05-09)
* Added support for TreeSummarizedExperiment and SpatialExperiment classes

# ExperimentSubset 0.99.0 (2020-10-02)
PRE-RELEASE
* Pre-release version of the `ExperimentSubset` package to be submitted to Bioconductor
