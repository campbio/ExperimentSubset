# Assays of subsets in ExperimentSubset objects

Get or set an assay of the main object or of a subset. Unlike the other
methods, `assay()` and `assay<-` name the subset with `i`, not with
`subsetName`.

`assay(x, i)` returns:

- a main assay, when `i` is its name, as for the class `x` extends;

- the subset's rows and columns of the assay it was made from, when `i`
  is a subset name. If the subset, or a subset it was made from, stored
  its own assay of that name with `unlink = TRUE`, that assay is
  returned instead;

- an assay as a subset sees it, when `i` is `"subset/assay"`, or a
  subset name with `subsetAssayName = "assay"`: either one stored in the
  subset, or one linked from its parent, cut to the subset's rows and
  columns.

`assay(x, i) <- value` stores `value` in a subset when `i` is
`"subset/assay"`, or a subset name with `subsetAssayName`, as
[`setSubsetAssay`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetAssay.md)
does. A name with a `"/"` is read as `"subset/assay"` unless it is the
name of a main assay, and a subset name without `subsetAssayName` is an
error. Any other name stores a main assay.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,character'
assay(x, i, withDimnames = TRUE, ...)

# S4 method for class 'SubsetSingleCellExperiment,character'
assay(x, i, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...)

# S4 method for class 'SubsetSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
assay(x, i, withDimnames = TRUE, ...)

# S4 method for class 'SubsetSpatialExperiment,character'
assay(x, i, withDimnames = TRUE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
assay(x, i, withDimnames = TRUE, ...) <- value
```

## Arguments

- x:

  An `ExperimentSubset` object.

- i:

  `character(1)` The name of a main assay, a subset name, or
  `"subset/assay"`.

- withDimnames:

  `logical(1)` As for `SummarizedExperiment`: `TRUE` gives the assay
  returned the row and column names of the main object or the subset,
  and `FALSE` keeps the names stored with the assay, if any. When
  setting a main assay, passed on to the method of the class `x`
  extends. Default `TRUE`.

- ...:

  `subsetAssayName`

  : `character(1)` With a subset name in `i`, the name of the assay to
    get from or store in that subset. It is an error with any other `i`,
    including `"subset/assay"`. Default `NULL`.

  `unlink`

  : `logical(1)` When storing in a subset, store the assay even though
    the subset sees an assay of that name linked from its parent, as for
    [`setSubsetAssay`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetAssay.md).
    It is an error when storing a main assay. Default `FALSE`.

  For a main assay, other arguments are passed on to the method of the
  class `x` extends.

- value:

  A matrix-like assay with the main object's dimensions, or, to store it
  in a subset, the subset's.

## Value

`assay()` returns a matrix-like assay. `assay<-` returns `x` with the
assay stored.

## See also

[`getSubsetAssay`](https://campbio.github.io/ExperimentSubset/dev/reference/getSubsetAssay.md),
[`setSubsetAssay`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetAssay.md),
[`assay`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

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
# the subset's 5 x 8 part of the counts
firstCounts <- assay(es, "first")
# store a new assay in the subset
assay(es, "first", subsetAssayName = "logcounts") <- log1p(firstCounts)
#> Storing assay inside subset 'first'.
# read it back
assay(es, "first/logcounts")
#>          cell1    cell2    cell3    cell4    cell5    cell6    cell7    cell8
#> gene1 1.609438 1.386294 2.302585 1.791759 2.079442 1.791759 2.197225 1.609438
#> gene2 1.609438 1.386294 1.386294 1.791759 1.945910 2.079442 1.609438 2.079442
#> gene3 1.791759 1.945910 1.945910 1.791759 2.079442 1.609438 1.791759 1.609438
#> gene4 2.197225 1.609438 1.386294 1.386294 1.791759 1.386294 1.609438 1.609438
#> gene5 1.386294 2.079442 1.609438 2.079442 1.791759 1.098612 1.945910 1.791759
```
