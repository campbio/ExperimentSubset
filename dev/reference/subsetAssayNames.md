# Name retrieval method for all subset assays in ExperimentSubset objects

Retrieves the names of all the subsets as well as the subset assays.

## Usage

``` r
subsetAssayNames(x, subsetName = NULL)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetAssayNames(x, subsetName = NULL)

# S4 method for class 'SubsetSingleCellExperiment'
subsetAssayNames(x, subsetName = NULL)

# S4 method for class 'SubsetSummarizedExperiment'
subsetAssayNames(x, subsetName = NULL)

# S4 method for class 'SubsetSpatialExperiment'
subsetAssayNames(x, subsetName = NULL)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetAssayNames(x, subsetName = NULL)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` If given, list the assays seen in this subset instead:
  its own assays and those linked from its parent. Default `NULL`.

## Value

A `vector` containing the names of the subsets followed by each assay
stored in a subset as `"subset/assay"`, or, with `subsetName`, the names
of the assays seen in that subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
assay(es, "subset1",
subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#> Storing assay inside subset 'subset1'.
subsetAssayNames(es)
#> [1] "subset1"               "subset1/subset1pAssay"
subsetAssayNames(es, subsetName = "subset1")
#> [1] "counts"        "logcounts"     "subset1pAssay"
```
