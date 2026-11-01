# Combine ExperimentSubset objects by columns or rows

`cbind(...)` and `rbind(...)` combine objects as for the class they
extend, and keep the subsets of every object, with their own assays,
rowData, colData, reducedDims and altExps. Subsets of the second and
later objects are moved to where that object's columns (for `cbind`) or
rows (for `rbind`) are in the result. Subset names must be unique across
the objects.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetRangedSummarizedExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSingleCellExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSingleCellExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSummarizedExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSpatialExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSpatialExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
rbind(..., deparse.level = 1)
```

## Arguments

- ...:

  `ExperimentSubset` objects of the same class, which must have the same
  rows for `cbind` or the same columns for `rbind`. Convert a plain
  object, such as a `SingleCellExperiment`, with
  [`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
  first, which keeps all of its data: combining it as it is makes R warn
  about "multiple direct matches", and if it comes first, the result is
  a plain object without subsets.

- deparse.level:

  Passed on to the method of the class the objects extend. Default `1`.

## Value

An object of the same class as the objects in `...`.

## Examples

``` r
# two objects with the same 10 genes and 20 cells each
set.seed(1)
counts1 <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts1) <- paste0("gene", 1:10)
colnames(counts1) <- paste0("cellA", 1:20)
counts2 <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts2) <- paste0("gene", 1:10)
colnames(counts2) <- paste0("cellB", 1:20)
es1 <- ExperimentSubset(list(counts = counts1))
es2 <- ExperimentSubset(list(counts = counts2))
# a subset in each object
es1 <- createSubset(es1, subsetName = "firstA", rows = 1:5, cols = 1:8)
es2 <- createSubset(es2, subsetName = "firstB", rows = 1:5, cols = 1:8)
# combine the cells of both objects
esBoth <- cbind(es1, es2)
# both subsets are kept
subsetNames(esBoth)
#> [1] "firstA" "firstB"
subsetColnames(esBoth, subsetName = "firstB")
#> [1] "cellB1" "cellB2" "cellB3" "cellB4" "cellB5" "cellB6" "cellB7" "cellB8"
```
