# ExperimentSubset 1.23.8 (2026-10-21)
* Subsets created from another subset (`parentAssay` set to a subset, or to an assay stored in one) now take their row and column names, `sample_id`, and the parent rowData, colData, rowLinks, colLinks, spatialCoords and spatialData from the right rows and columns of the main object. Before, the indices relative to the parent subset were applied to the main object, so these came from the wrong features and cells while `assay()` returned the right values. The parent metadata is looked up on access, so it is now right for existing objects too; the subset's row and column names (and `sample_id`) are stored when the subset is created, so objects saved with an earlier version keep the wrong names for such subsets. Repair them with `subsetRownames(x, "s2") <- rownames(assay(x, "s2"))` and the same for `subsetColnames()`; a wrong SpatialExperiment `sample_id` needs the subset re-created

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
