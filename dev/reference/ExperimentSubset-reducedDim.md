# Reduced dimensions of subsets in ExperimentSubset objects

Get or set the reduced dimensions, such as a PCA or UMAP, of a subset by
naming it with `subsetName`. Without `subsetName`, these methods act on
the main object, as for `SingleCellExperiment`. They exist for the
classes that extend `SingleCellExperiment`:
`SubsetSingleCellExperiment`, `SubsetSpatialExperiment` and
`SubsetTreeSummarizedExperiment`. `reducedDims<-`,
[`reducedDimNames`](https://campbio.github.io/ExperimentSubset/dev/reference/reducedDimNames.md)
and `reducedDimNames<-` have their own pages.

## Usage

``` r
# S4 method for class 'SubsetSingleCellExperiment,character'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSingleCellExperiment,missing'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSingleCellExperiment,numeric'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSingleCellExperiment,character'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,missing'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,numeric'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment'
reducedDims(x, ...)

# S4 method for class 'SubsetSpatialExperiment,character'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSpatialExperiment,missing'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSpatialExperiment,numeric'
reducedDim(x, type, ...)

# S4 method for class 'SubsetSpatialExperiment,character'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,missing'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,numeric'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment'
reducedDims(x, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
reducedDim(x, type, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,missing'
reducedDim(x, type, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,numeric'
reducedDim(x, type, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,missing'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,numeric'
reducedDim(x, type, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
reducedDims(x, ...)
```

## Arguments

- x:

  An `ExperimentSubset` object.

- type:

  `character(1)` or `integer(1)` The name or position of the reduced
  dimension. Leave it out for the first one.

- ...:

  `subsetName`

  : `character(1)` Name of the subset whose reduced dimensions to get or
    set. Leave it out for the main object's.

  Other arguments are passed on to the `SingleCellExperiment` method,
  for the subset or the main object.

- withDimnames:

  `logical(1)` Passed on to `reducedDim<-`. Default `TRUE`.

- value:

  A matrix with one row per column of the subset (or of the main
  object).

## Value

`reducedDim()` returns a matrix and `reducedDims()` a list of them.
`reducedDim<-` returns `x` with the reduced dimension stored.

## See also

[`reducedDims`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)

## Examples

``` r
# an object with 10 genes and 20 cells
set.seed(1)
counts <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts) <- paste0("gene", 1:10)
colnames(counts) <- paste0("cell", 1:20)
es <- ExperimentSubset(list(counts = counts))
# a subset of the first 5 genes and the first 8 cells
es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
# a 2-dimensional embedding of the subset's 8 cells
embedding <- matrix(rnorm(16), ncol = 2)
# store it in the subset
reducedDim(es, "PCA", subsetName = "first") <- embedding
# read it back
reducedDim(es, "PCA", subsetName = "first")
#>              [,1]       [,2]
#> cell1 -0.62036668  0.3841854
#> cell2  0.04211587  1.6821761
#> cell3 -0.91092165 -0.6357365
#> cell4  0.15802877 -0.4616447
#> cell5 -0.65458464  1.4322822
#> cell6  1.76728727 -0.6506964
#> cell7  0.71670748 -0.2073807
#> cell8  0.91017423 -0.3928079
# all of the subset's reduced dimensions
reducedDims(es, subsetName = "first")
#> List of length 1
#> names(1): PCA
```
