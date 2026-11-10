# ExperimentSubset
[![Bioc release build](https://bioconductor.org/shields/build/release/bioc/ExperimentSubset.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/ExperimentSubset/)
[![Bioc devel build](https://bioconductor.org/shields/build/devel/bioc/ExperimentSubset.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/ExperimentSubset/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/ExperimentSubset.svg)](https://bioconductor.org/packages/ExperimentSubset)
[![Bioconductor downloads](https://bioconductor.org/shields/downloads/release/ExperimentSubset.svg)](https://bioconductor.org/packages/stats/bioc/ExperimentSubset/)
[![R-CMD-check](https://github.com/campbio/ExperimentSubset/actions/workflows/R-CMD-check.yaml/badge.svg?branch=devel)](https://github.com/campbio/ExperimentSubset/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/campbio/ExperimentSubset/actions/workflows/pkgdown.yaml/badge.svg?branch=devel)](https://github.com/campbio/ExperimentSubset/actions/workflows/pkgdown.yaml)

## Introduction
`ExperimentSubset` package enables users to perform flexible subsetting of Single-Cell data that comes from the same experiment as well as the consequent storage of these subsets back into the same object. It offers the same interface to the users as other `Experiment` classes but with additional subset support.

Documentation, examples and articles are on the [ExperimentSubset website](https://campbio.github.io/ExperimentSubset/), including a [workflow article](https://campbio.github.io/ExperimentSubset/dev/articles/pbmc-workflow.html) that keeps a whole RNA and protein (CITE-seq) analysis of 5,000 PBMCs, from raw droplets to T-cell subtypes, in one object.

## Supported Experiment Classes
1. SummarizedExperiment
2. RangedSummarizedExperiment
3. SingleCellExperiment
4. TreeSummarizedExperiment
5. SpatialExperiment

## Installation
To install the release version from
[Bioconductor](https://bioconductor.org/packages/ExperimentSubset):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ExperimentSubset")
```

To install the development version from the `devel` branch on GitHub
(it may need the Bioconductor devel versions of its dependencies):
```r
BiocManager::install("campbio/ExperimentSubset@devel")
```

## Citation
Irzam Sarfraz, Muhammad Asif, Joshua D Campbell, ExperimentSubset: an R package to manage subsets of Bioconductor Experiment objects, Bioinformatics, 2021;, btab179, https://doi.org/10.1093/bioinformatics/btab179
