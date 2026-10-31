# Subset parent hierarchy retrieval method for ExperimentSubset objects

Retrieves a complete 'subset to parent' link from a specified subset.

## Usage

``` r
subsetParent(x, subsetName)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetParent(x, subsetName)

# S4 method for class 'SubsetSingleCellExperiment'
subsetParent(x, subsetName)

# S4 method for class 'SubsetSummarizedExperiment'
subsetParent(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment'
subsetParent(x, subsetName)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetParent(x, subsetName)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Specify the name of the subset, or an assay stored in a
  subset as `"subset/assay"`, against which the 'subset to parent link'
  should be retrieved.

## Value

A `list` containing the 'subset to parent' link.

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
subsetParent(es, "subset1/subset1pAssay")
#> [[1]]
#> [1] "subset1"
#> 
#> [[2]]
#> [1] "counts"
#> 
```
