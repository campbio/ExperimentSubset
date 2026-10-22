# ExperimentSubset
[![R-CMD-check](https://github.com/campbio/ExperimentSubset/actions/workflows/R-CMD-check.yaml/badge.svg?branch=devel)](https://github.com/campbio/ExperimentSubset/actions/workflows/R-CMD-check.yaml)

## Introduction
`ExperimentSubset` package enables users to perform flexible subsetting of Single-Cell data that comes from the same experiment as well as the consequent storage of these subsets back into the same object. It offers the same interface to the users as other `Experiment` classes but with additional subset support. More details about the package are available in the package vignette.

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

To install the latest development version from GitHub:
```r
BiocManager::install("campbio/ExperimentSubset")
```

## Citation
Irzam Sarfraz, Muhammad Asif, Joshua D Campbell, ExperimentSubset: an R package to manage subsets of Bioconductor Experiment objects, Bioinformatics, 2021;, btab179, https://doi.org/10.1093/bioinformatics/btab179
