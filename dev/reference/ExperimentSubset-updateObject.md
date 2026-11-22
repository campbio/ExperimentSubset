# Update an ExperimentSubset object saved by an older version

`updateObject()` brings an `ExperimentSubset` object saved by an older
version of the package up to date: the class it extends updates the main
object, and each subset is updated too. Since version 1.23.55, each
subset records whether it was made from an altExp with
`createSubset(altExp =)`; subsets of objects saved before then record
that they weren't. Such objects also work without `updateObject()`.

## Usage

``` r
# S4 method for class 'SubsetSummarizedExperiment'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'SubsetRangedSummarizedExperiment'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'SubsetSingleCellExperiment'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'SubsetTreeSummarizedExperiment'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'SubsetSpatialExperiment'
updateObject(object, ..., verbose = FALSE)

# S4 method for class 'AssaySubset'
updateObject(object, ..., verbose = FALSE)
```

## Arguments

- object:

  An `ExperimentSubset` object.

- ...:

  Passed on to `updateObject()` of the class the object extends, and of
  each subset's stored object.

- verbose:

  `logical(1)` Whether to say what is updated. Default `FALSE`.

## Value

`object`, updated.

## See also

[`updateObject`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html)

## Examples

``` r
# an object with a subset
counts <- matrix(rpois(60, lambda = 5), nrow = 6)
rownames(counts) <- paste0("gene", 1:6)
colnames(counts) <- paste0("cell", 1:10)
es <- ExperimentSubset(list(counts = counts))
es <- createSubset(es, subsetName = "first", rows = 1:4, cols = 1:8)
# an object saved by an older version is updated the same way
es <- updateObject(es)
subsetNames(es)
#> [1] "first"
```
