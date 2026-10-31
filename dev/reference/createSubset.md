# Subset creation method for ExperimentSubset objects

Create a subset from an already available `assay` in the input
`ExperimentSubset` object by specifying the rows and columns to include
in the subset.

## Usage

``` r
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)

# S4 method for class 'SubsetRangedSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)

# S4 method for class 'SubsetSingleCellExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)

# S4 method for class 'SubsetSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)

# S4 method for class 'SubsetSpatialExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)

# S4 method for class 'SubsetTreeSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL)
```

## Arguments

- x:

  `ExperimentSubset` Specify the object from which a subset should be
  created. Use
  [`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
  first to convert another experiment object.

- subsetName:

  `character(1)` Specify the name of the subset to create.

- rows:

  `numeric` or `character` vector. Specify the rows to include in this
  subset. If `missing` or `NULL`, all rows are included in the subset.
  Default `NULL`.

- cols:

  `numeric` or `character` vector. Specify the columns to include in
  this subset. If `missing` or `NULL`, all columns are included in the
  subset. Default `NULL`.

- parentAssay:

  `character(1)` Specify the parent `assay` of the subset: the name of a
  main `assay`, of a subset, or of an assay seen in a subset written as
  `"subset/assay"`. If `NULL`, the first available main `assay` will be
  marked as parent. Default `NULL`.

## Value

An `ExperimentSubset` object that now contains the newly created subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
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
#> subsets(1): subset1
#> subsetAssays(1): subset1
```
