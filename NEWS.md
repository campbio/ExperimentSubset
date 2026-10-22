# ExperimentSubset 1.23.13 (2026-10-22)
* Documentation fixes: the vignette lists all five supported classes, names the `createSubset()` arguments `rows` and `cols`, seeds its random data, drops the pin to Bioconductor 3.11 and installs from GitHub with BiocManager; help pages say which input `ExperimentSubset()` and `createSubset()` accept, fix the `subsetRownames()` page, and add examples for `subsetRowLinks()`, `subsetColLinks()` and `subsetSpatialCoords()`; examples no longer need scater
* DESCRIPTION gains URL and BugReports fields and drops the stale Date and the unused LazyData; the citation uses `bibentry()`

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
