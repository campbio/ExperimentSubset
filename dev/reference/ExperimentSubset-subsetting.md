# Subset an ExperimentSubset object by rows or columns

`x[i, j]` keeps rows `i` and columns `j` of the main object, as for the
class it extends, and also updates every subset. Each subset keeps the
rows and columns that remain, in its own order, together with its own
assays, rowData, colData, reducedDims and altExps. A subset none of
whose rows or columns remain is kept with no rows or no columns. If `i`
or `j` repeats a row or column, subsets keep one copy of it.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SubsetSingleCellExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SubsetSummarizedExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SubsetSpatialExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'SubsetTreeSummarizedExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]
```

## Arguments

- x:

  An `ExperimentSubset` object.

- i, j:

  Rows and columns to keep: positions, names or logical vectors, as for
  the class `x` extends. Leave one out to keep all rows or all columns.

- ...:

  Passed on to the method of the class `x` extends.

- drop:

  Ignored, as for `SummarizedExperiment`. Default `TRUE`.

## Value

An object of the same class as `x`.

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
# keep genes 3 to 10 of the main object
esSmaller <- es[3:10, ]
# the subset keeps genes 3 to 5
subsetRownames(esSmaller, subsetName = "first")
#> [1] "gene3" "gene4" "gene5"
```
