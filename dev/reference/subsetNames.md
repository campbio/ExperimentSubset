# Get names of only the subsets in ExperimentSubset objects

Retrieves the names of the available subsets (not the subset assays) in
an `ExperimentSubset` object.

## Usage

``` r
subsetNames(x)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetNames(x)

# S4 method for class 'SubsetSingleCellExperiment'
subsetNames(x)

# S4 method for class 'SubsetSummarizedExperiment'
subsetNames(x)

# S4 method for class 'SubsetSpatialExperiment'
subsetNames(x)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetNames(x)
```

## Arguments

- x:

  `ExperimentSubset` Specify the input ES object.

## Value

A `vector` of subset names.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetNames(es)
#> [1] "subset1"
```
