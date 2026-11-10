# ExperimentSubset

## Introduction

`ExperimentSubset` package enables users to perform flexible subsetting
of Single-Cell data that comes from the same experiment as well as the
consequent storage of these subsets back into the same object. It offers
the same interface to the users as other `Experiment` classes but with
additional subset support.

Documentation, examples and articles are on the [ExperimentSubset
website](https://campbio.github.io/ExperimentSubset/), including a
[workflow
article](https://campbio.github.io/ExperimentSubset/dev/articles/pbmc-workflow.html)
that keeps a whole RNA and protein (CITE-seq) analysis of 5,000 PBMCs,
from raw droplets to T-cell subtypes, in one object.

## Supported Experiment Classes

1.  SummarizedExperiment
2.  RangedSummarizedExperiment
3.  SingleCellExperiment
4.  TreeSummarizedExperiment
5.  SpatialExperiment

## Installation

To install the release version from
[Bioconductor](https://bioconductor.org/packages/ExperimentSubset):

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"BiocManager"``, quietly ``=`` ``TRUE``)``)`\
`    `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"BiocManager"``)`\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"ExperimentSubset"``)`

To install the development version from the `devel` branch on GitHub (it
may need the Bioconductor devel versions of its dependencies):

\
`BiocManager``::`[`install`](https://bioconductor.github.io/BiocManager/reference/install.html)`(``"campbio/ExperimentSubset@devel"``)`

## Citation

Irzam Sarfraz, Muhammad Asif, Joshua D Campbell, ExperimentSubset: an R
package to manage subsets of Bioconductor Experiment objects,
Bioinformatics, 2021;, btab179,
<https://doi.org/10.1093/bioinformatics/btab179>
