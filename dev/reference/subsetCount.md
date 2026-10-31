# Subset count method for ExperimentSubset objects

Get the total count of the available subsets (excluding subset assays)
in an `ExperimentSubset` object.

## Usage

``` r
subsetCount(x)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetCount(x)

# S4 method for class 'SubsetSingleCellExperiment'
subsetCount(x)

# S4 method for class 'SubsetSummarizedExperiment'
subsetCount(x)

# S4 method for class 'SubsetSpatialExperiment'
subsetCount(x)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetCount(x)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

## Value

A `numeric` value representing the total count of the subsets.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetCount(es)
#> [1] 1
```
