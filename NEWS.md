# ExperimentSubset 1.23.9 (2026-10-21)
* A subset with a single row or column now returns a one-row or one-column matrix from `assay()`, not a vector, so subsets can be made from it
* `setSubsetAssay()` (and `assay(x, i) <- value` for a subset) no longer replaces the subset's row and column names with the input matrix's. A matrix without names is stored under the subset's names; a matrix whose names differ from the subset's, or are in another order, is now an error instead of silently relabelling the subset; a matrix of the wrong size gets an error that says so. `assay(x, subsetName)` now returns the subset's own row and column names, including any set with `subsetRownames<-` or `subsetColnames<-`, so `setSubsetAssay(x, s, assay(x, s), ...)` always works
* A `subsetName` with spaces is now stored under the name with spaces removed, which the warning already announced; before, the subset was stored under the original name and lookups of its stored assays failed. Refer to such a subset by the name without spaces

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
