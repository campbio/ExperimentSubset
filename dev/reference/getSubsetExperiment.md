# Take a subset out as a standalone experiment object

Returns a subset as an object of the main object's class (for example a
`SingleCellExperiment`), so that any tool written for that class can be
run on the subset. Store what the tool added with `setSubsetExperiment`.

## Usage

``` r
getSubsetExperiment(x, subsetName, parentData = TRUE)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
getSubsetExperiment(x, subsetName, parentData = TRUE)

# S4 method for class 'SubsetSingleCellExperiment,character'
getSubsetExperiment(x, subsetName, parentData = TRUE)

# S4 method for class 'SubsetSummarizedExperiment,character'
getSubsetExperiment(x, subsetName, parentData = TRUE)

# S4 method for class 'SubsetSpatialExperiment,character'
getSubsetExperiment(x, subsetName, parentData = TRUE)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
getSubsetExperiment(x, subsetName, parentData = TRUE)
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  `character(1)` The name of the subset.

- parentData:

  `logical(1)` Include the main object's rowData and colData columns for
  the subset's rows and columns, as with `parentRowData` and
  `parentColData`. They are marked so that `setSubsetExperiment` doesn't
  store copies of them. A subset's own column hides the main object's
  column of the same name, with a message if their values differ.
  Default `TRUE`.

## Value

An object of the main object's class with the subset's rows and columns
and row and column names. Its assays are every assay the subset sees,
its own and those linked from its parent, under their names, with the
subset's primary assay first. Its rowData and colData are the subset's,
after the main object's columns if `parentData = TRUE`. Its reducedDims,
altExps and metadata are the subset's. Everything else, such as row
ranges, row and column trees, spatial coordinates, images, `colPairs`
and `rowPairs`, comes from the main object.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1", rows = 1:100, cols = 1:50)
subset1 <- getSubsetExperiment(es, "subset1")
subset1
#> class: SingleCellExperiment 
#> dim: 100 50 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(100): C12orf73 RNU6-1256P ... GABRA2 RP11-563K23.1
#> rowData names(1): gene
#> colnames(50): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... ACACTGATCGAATCCA
#>   ACGTCAACAGGTCGTC
#> colData names(11): nGene nUMI ... ident doublet_true_labels
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
