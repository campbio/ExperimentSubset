# ExperimentSubset
[![Build Status](https://travis-ci.org/campbio/ExperimentSubset.svg?branch=master)](https://travis-ci.org/campbio/ExperimentSubset)
[![Codecov test coverage](https://codecov.io/gh/campbio/ExperimentSubset/branch/master/graph/badge.svg)](https://codecov.io/gh/campbio/ExperimentSubset?branch=master)

## Introduction
`ExperimentSubset` package enables users to perform flexible subsetting of Single-Cell data that comes from the same experiment as well as the consequent storage of these subsets back into the same object. It offers the same interface to the users as other `Experiment` classes but with additional subset support. More details about the package are available in the package vignette.

<b> Note: Bioconductor 3.12 release has an old version of the `ExperimentSubset` and we suggest to use the current `master` branch of the this github repo to use the latest updated version of `ExperimentSubset` package. </b>

<b> Note: `ExperimentSubset` makes use of development version of `SpatialExperiment` package which should be installed from https://github.com/irzamsarfraz/SpatialExperiment. </b>

## Supported Experiment Classes
1. SummarizedExperiment
2. RangedSummarizedExperiment
3. SingleCellExperiment
4. TreeSummarizedExperiment
5. SpatialExperiment

## Installation
To install the latest version of `ExperimentSubset` from GitHub using `devtools`:
```
library(devtools)
install_github("campbio/ExperimentSubset")
```

## Citation
Irzam Sarfraz, Muhammad Asif, Joshua D Campbell, ExperimentSubset: an R package to manage subsets of Bioconductor Experiment objects, Bioinformatics, 2021;, btab179, https://doi.org/10.1093/bioinformatics/btab179
