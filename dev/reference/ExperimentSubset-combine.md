# Combine ExperimentSubset objects by columns or rows

`cbind(...)` and `rbind(...)` combine objects as for the class they
extend, and keep the subsets of every object, with their own assays,
rowData, colData, reducedDims and altExps. Subsets of the second and
later objects are moved to where that object's columns (for `cbind`) or
rows (for `rbind`) are in the result. Subset names must be unique across
the objects, and differ from the assay names of all of them.

`combineCols(x, ...)` and `combineRows(x, ...)` combine objects whose
rows (for `combineCols`) or columns (for `combineRows`) differ, as for
the class they extend: the result has all of the rows (or columns),
matched by name, and fills the gaps. They keep the subsets of every
object too, matching each subset's rows (or columns) by name to the
result's, or by position with `use.names = FALSE`.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetRangedSummarizedExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetRangedSummarizedExperiment'
combineCols(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetRangedSummarizedExperiment'
combineRows(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSingleCellExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSingleCellExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSingleCellExperiment'
combineCols(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSingleCellExperiment'
combineRows(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSummarizedExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSummarizedExperiment'
combineCols(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSummarizedExperiment'
combineRows(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSpatialExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetSpatialExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetSpatialExperiment'
combineCols(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetSpatialExperiment'
combineRows(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetTreeSummarizedExperiment'
cbind(..., deparse.level = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
rbind(..., deparse.level = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
combineCols(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)

# S4 method for class 'SubsetTreeSummarizedExperiment'
combineRows(x, ..., delayed = TRUE, fill = NA, use.names = TRUE)
```

## Arguments

- ...:

  `ExperimentSubset` objects of the same class, which must have the same
  rows for `cbind` or the same columns for `rbind`. Convert a plain
  object, such as a `SingleCellExperiment`, with
  [`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
  first, which keeps all of its data: combining it as it is makes R warn
  about "multiple direct matches", and if it comes first, the result is
  a plain object without subsets. For `combineCols` and `combineRows`,
  the objects to combine with `x`.

- deparse.level:

  Passed on to the method of the class the objects extend. Default `1`.

- x:

  For `combineCols` and `combineRows`, the first `ExperimentSubset`
  object to combine.

- delayed, fill:

  Passed on to the method of the class the objects extend: whether to
  combine assays as `DelayedArray` objects, and the value for the gaps.
  Defaults `TRUE` and `NA`.

- use.names:

  Whether to match rows (for `combineCols`) or columns (for
  `combineRows`) by name, as the method of the class the objects extend
  does; if `FALSE`, they are matched by position and the objects must
  have the same number of them. Default `TRUE`.

## Value

For `cbind` and `rbind`, an object of the same class as the objects in
`...`. For `combineCols` and `combineRows`, the `ExperimentSubset` class
matching the result of the method of the class `x` extends; see Details.

## Details

For some classes, `combineCols()` and `combineRows()` of the class the
objects extend return a different class: `combineCols()` of a
`TreeSummarizedExperiment` or `SpatialExperiment` returns a
`SingleCellExperiment`, and `combineRows()` of a `SingleCellExperiment`,
`TreeSummarizedExperiment` or `SpatialExperiment` returns a
`RangedSummarizedExperiment`. The result is then the matching
`ExperimentSubset` class, such as `SubsetSingleCellExperiment`, and each
subset is converted to that class too, so it loses what the class can't
hold, as the main object does: for example its reducedDims and altExps
in a `SubsetRangedSummarizedExperiment`.

`combineCols()` and `combineRows()` use the method for the class of `x`
only, so if `x` is a plain object, such as a `SingleCellExperiment`, the
result is a plain object without subsets. Convert it with
[`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
first.

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

# a third object, with genes 6 to 10 and 5 genes the others don't have
counts3 <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts3) <- paste0("gene", 6:15)
colnames(counts3) <- paste0("cellC", 1:20)
es3 <- ExperimentSubset(list(counts = counts3))
# a subset of the genes only the third object has
es3 <- createSubset(es3, subsetName = "newGenes", rows = 6:10)
# combine the cells, keeping the genes of both; the gaps are NA
esAll <- combineCols(es1, es3)
dim(esAll)
#> [1] 15 40
# the subsets of both objects are kept
subsetNames(esAll)
#> [1] "firstA"   "newGenes"
# the new subset keeps its genes
subsetRownames(esAll, subsetName = "newGenes")
#> [1] "gene11" "gene12" "gene13" "gene14" "gene15"
```
