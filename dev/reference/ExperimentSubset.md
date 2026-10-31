# ExperimentSubset constructor

This constructor function sets up an `ExperimentSubset` object from a
`SummarizedExperiment`, `RangedSummarizedExperiment`,
`SingleCellExperiment`, `TreeSummarizedExperiment` or
`SpatialExperiment` object, or from a named `list` of assay matrices. A
subset can also be created directly by passing a named `list` to the
`subset` parameter, with elements named `subsetName`, `rows`, `cols` and
`parentAssay`.

## Usage

``` r
ExperimentSubset(
  x,
  subset = list(subsetName = NA, rows = NA, cols = NA, parentAssay = NA)
)
```

## Arguments

- x:

  An experiment object of one of the five supported classes (subclasses
  of them are not supported), or a named `list` of assay matrices, which
  is passed to the `SingleCellExperiment` constructor as its assays.

- subset:

  A named `list` to create a subset from within the constructor, with
  elements `subsetName` and, optionally, `rows`, `cols` and
  `parentAssay`. Leave out `rows`, `cols` or `parentAssay`, or set them
  to `NULL`, to use all rows, all columns or the first assay. The
  default, a list of `NA` values, creates no subset; an `NA` in a list
  you supply is an error.

## Value

An object of the matching `Subset*` class, for example a
`SubsetSingleCellExperiment` for a `SingleCellExperiment` or a list of
assays.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es
#> class: SubsetSingleCellExperiment 
#> dim: 2000 2000 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3
#> rowData names(1): gene
#> colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT
#>   ATTTCTGGTGATGATA
#> colData names(11): nGene nUMI ... ident doublet_true_labels
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> subsets(0): 
#> subsetAssays(0): 
```
