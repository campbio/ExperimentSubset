# Count method for subset assays in ExperimentSubset objects

Get the count of the total available subsets and the subset assays
inside these subsets in an `ExperimentSubset` object.

## Usage

``` r
subsetAssayCount(x)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetAssayCount(x)

# S4 method for class 'SubsetSingleCellExperiment'
subsetAssayCount(x)

# S4 method for class 'SubsetSummarizedExperiment'
subsetAssayCount(x)

# S4 method for class 'SubsetSpatialExperiment'
subsetAssayCount(x)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetAssayCount(x)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

## Value

A `numeric` value representing the sum of the subset count and subset
assay count.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
assay(es, "subset1",
subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#> Storing assay inside subset 'subset1'.
subsetAssayCount(es)
#> [1] 2
```
