# Get dimensions of subsets in ExperimentSubset objects

Retrieves the dimensions of the specified subset in an
`ExperimentSubset` object.

## Usage

``` r
subsetDim(x, subsetName)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetDim(x, subsetName)

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetDim(x, subsetName)

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetDim(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetDim(x, subsetName)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetDim(x, subsetName)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to retrieve the dimensions from.

## Value

A `vector` containing the dimensions of the specified subset i.e. the
number of rows and the number of columns in the subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetDim(es, "subset1")
#> [1] 10 10
```
