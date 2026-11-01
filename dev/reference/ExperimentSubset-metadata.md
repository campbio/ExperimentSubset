# Metadata of subsets in ExperimentSubset objects

Get or set the `metadata` list of a subset by naming it with
`subsetName`. Without `subsetName`, `metadata()` and `metadata<-` act on
the main object, as for the class `x` extends.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment'
metadata(x, ...)

# S4 method for class 'SubsetRangedSummarizedExperiment'
metadata(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment'
metadata(x, ...)

# S4 method for class 'SubsetSingleCellExperiment'
metadata(x, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment'
metadata(x, ...)

# S4 method for class 'SubsetSummarizedExperiment'
metadata(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment'
metadata(x, ...)

# S4 method for class 'SubsetSpatialExperiment'
metadata(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
metadata(x, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment'
metadata(x, ...) <- value
```

## Arguments

- x:

  An `ExperimentSubset` object.

- ...:

  `subsetName`

  : `character(1)` Name of the subset whose metadata to get or set.
    Leave it out for the main object's.

  Without `subsetName`, all arguments are passed on to the method of the
  class `x` extends.

- value:

  A `list` to store as the metadata.

## Value

`metadata()` returns a `list`. `metadata<-` returns `x` with the
metadata set.

## See also

[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)

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
# store a note in the subset's metadata
metadata(es, subsetName = "first") <- list(note = "the first 8 cells")
# read it back
metadata(es, subsetName = "first")
#> $note
#> [1] "the first 8 cells"
#> 
# the main object's metadata is unchanged
metadata(es)
#> list()
```
