# ExperimentSubset 1.23.10 (2026-10-21)
* `createSubset()` checks its input more strictly: `rows` and `cols` must select at least one row or column and be positive whole numbers or names, without duplicates (duplicates were accepted in some cases before); names given for an assay without row or column names give a clear error; `subsetName` must be a single non-empty string that is not already the name of an assay or of an assay stored in a subset; `parentAssay` must be a single string
* `ExperimentSubset()` now stops when `subset` contains NA instead of silently skipping the subset, and rejects a list of assays without names
* Unknown subset names now give a clear error in `subsetDim()`, `subsetColnames()`, `subsetRownames()` and their setters, and `reducedDims<-`; `assay(x, i)` and `getSubsetAssay(x, i)` with an unknown name are an error instead of returning `NULL`
* `assay(x, i) <- value` stores the assay inside the subset whenever `i` is a subset name, and then needs `subsetAssayName` (giving `subsetAssayName` with a name that is not a subset is an error); before, the choice depended on the dimensions of `value`, so a subset with the same dimensions as the main object got a new main assay named after it

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
