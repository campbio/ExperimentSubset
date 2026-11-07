# Store a standalone experiment object back into a subset

Stores what a tool added to an object returned by `getSubsetExperiment`:
new assays, rowData and colData columns, reducedDims, altExps and
metadata.

## Usage

``` r
setSubsetExperiment(x, subsetName, value, unlink = FALSE)

# S4 method for class 'SubsetRangedSummarizedExperiment'
setSubsetExperiment(x, subsetName, value, unlink = FALSE)

# S4 method for class 'SubsetSingleCellExperiment'
setSubsetExperiment(x, subsetName, value, unlink = FALSE)

# S4 method for class 'SubsetSummarizedExperiment'
setSubsetExperiment(x, subsetName, value, unlink = FALSE)

# S4 method for class 'SubsetSpatialExperiment'
setSubsetExperiment(x, subsetName, value, unlink = FALSE)

# S4 method for class 'SubsetTreeSummarizedExperiment'
setSubsetExperiment(x, subsetName, value, unlink = FALSE)
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  `character(1)` The name of the subset.

- value:

  An object of the main object's class (not a subclass or a parent
  class) with the subset's dimensions and the subset's row and column
  names, in the same order. A tool that drops rows or columns is
  followed by `createSubset` instead.

- unlink:

  `TRUE`, `FALSE` or a character vector of `assay` names. Store assays
  of `value` that have the name of an assay linked from the subset's
  parent but different values, as with `setSubsetAssay`: with `TRUE`,
  every such assay; with a character vector, only those it names. An
  assay of an altExp is named `"altExpName/assayName"`, such as
  `"ADT/counts"`. Any other such assay, of the main object or an altExp,
  is an error, and nothing is stored. A message names each assay that is
  unlinked. Names of assays that are not linked, or unchanged, are
  ignored, so `unlink = "decontXcounts"` stores a changed
  `decontXcounts` but stops if `counts` was changed too. Default
  `FALSE`.

## Value

Updated `ExperimentSubset` object. Linked assays that are in `value`
hold the same values as what the subset sees (whatever their dimnames,
type or class) are not stored, and neither are the main object's rowData
and colData columns added by `getSubsetExperiment`. altExps are stored
as `altExp<-` stores them (see
[`ExperimentSubset-altExp`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-altExp.md)):
only the assays and data that differ from a linked altExp, and a linked
altExp missing from `value` stays linked. Everything else in `value`
replaces the subset's own data: an own assay, reducedDim, altExp or
column that `value` no longer has is removed from the subset, except an
assay another subset was made from, which is an error. Nothing else is
stored: not row ranges, row and column trees, spatial coordinates,
images, `colPairs`, `rowPairs` or `mainExpName`.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1", rows = 1:100, cols = 1:50)
subset1 <- getSubsetExperiment(es, "subset1")
# add a reducedDim and a colData column to the standalone object
reducedDim(subset1, "random") <- matrix(rnorm(100), nrow = 50)
subset1$group <- rep(c("a", "b"), 25)
# and store them in the subset
es <- setSubsetExperiment(es, "subset1", subset1)
reducedDimNames(es, subsetName = "subset1")
#> [1] "random"
```
