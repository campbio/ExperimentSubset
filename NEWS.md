# ExperimentSubset 1.23.7 (2026-10-21)
* Accessors called on the main object without `subsetName` now pass extra arguments on to the base method; for example `reducedDim(x, "PCA", withDimnames = FALSE)`, `reducedDims(x, withDimnames = FALSE)` and `altExps(x, withDimnames = FALSE)` failed with "argument x is missing"

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
