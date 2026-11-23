# Remove subsets, and what happens to the subsets made from them

`removeSubset` removes subsets from an `ExperimentSubset` object. Its
`dependents` argument says what happens to the subsets made from them,
and the same argument does so wherever what a subset was made from can
go away:

- a main assay, removed with `assay(x, i) <- NULL` or
  `assays(x) <- value`;

- an assay of a subset, removed with `setSubsetExperiment`;

- an assay of an altExp, removed with `altExp<-`, `altExps<-` or
  `setSubsetExperiment`, or lost when the main altExp gets other rows
  and what a subset stored in it no longer fits;

- an altExp that subsets were made from with `createSubset(altExp =)`,
  removed with `altExp<-` or `altExps<-`;

- a subset, removed with `removeSubset`.

## Usage

``` r
removeSubset(x, subsetName, dependents = "error")

# S4 method for class 'SubsetSummarizedExperiment'
removeSubset(x, subsetName, dependents = "error")

# S4 method for class 'SubsetRangedSummarizedExperiment'
removeSubset(x, subsetName, dependents = "error")

# S4 method for class 'SubsetSingleCellExperiment'
removeSubset(x, subsetName, dependents = "error")

# S4 method for class 'SubsetTreeSummarizedExperiment'
removeSubset(x, subsetName, dependents = "error")

# S4 method for class 'SubsetSpatialExperiment'
removeSubset(x, subsetName, dependents = "error")
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  `character` The names of one or more subsets to remove.

- dependents:

  `"error"`, `"keep"` or `"remove"`: what happens to the subsets made
  from what is removed.

  - `"error"` stops, naming them, and changes nothing.

  - `"keep"` keeps them intact: each stores its own copy of what it saw
    of what is removed, with a message such as "Subset 's' now stores
    its own assay 'counts'.", and sees the same data as before. That
    copy isn't linked: an assay of the same name added later doesn't
    replace it. For `removeSubset`, a subset made from a removed subset
    is then made from that subset's parent instead, with its own copy of
    the assays, altExps, spatial coordinates and trees it saw from the
    removed subset; one made with `createSubset(altExp =)` from the
    removed subset's altExp can only be kept if its parent's altExp has
    the same rows. When a whole altExp is removed, `"keep"` is not
    supported yet.

  - `"remove"` removes them, and every subset made from them, with a
    warning naming them.

  A subset that only saw a removed assay linked from its parent, but was
  made from another one, isn't affected: it no longer sees it. Default
  `"error"`.

## Value

Updated `ExperimentSubset` object, of the same class.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, subsetName = "subset1", rows = 1:100, cols = 1:50)
es <- createSubset(es, subsetName = "subset2", rows = 1:10, cols = 1:5,
                   parentAssay = "subset1")
# subset2 was made from subset1, so removing subset1 stops unless
# dependents says what happens to subset2
es2 <- removeSubset(es, subsetName = "subset1", dependents = "keep")
subsetNames(es2)
#> [1] "subset2"
# the same for a main assay that subsets were made from
assay(es, "counts", dependents = "keep") <- NULL
#> Subset 'subset1' now stores its own assay 'counts'.
subsetAssayNames(es, subsetName = "subset1")
#> [1] "logcounts" "counts"   
```
