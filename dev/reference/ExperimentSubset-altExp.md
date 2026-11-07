# Alternative experiments of subsets in ExperimentSubset objects

Get or set the alternative experiments (`altExps`) of a subset by naming
it with `subsetName`. Without `subsetName`, these methods act on the
main object, as for `SingleCellExperiment`. They exist for the classes
that extend `SingleCellExperiment`: `SubsetSingleCellExperiment`,
`SubsetSpatialExperiment` and `SubsetTreeSummarizedExperiment`.

## Usage

``` r
# S4 method for class 'SubsetSingleCellExperiment,character'
altExp(x, e, ...)

# S4 method for class 'SubsetSingleCellExperiment,missing'
altExp(x, e, ...)

# S4 method for class 'SubsetSingleCellExperiment,numeric'
altExp(x, e, ...)

# S4 method for class 'SubsetSingleCellExperiment'
altExps(x, ...)

# S4 method for class 'SubsetSingleCellExperiment'
altExpNames(x, ...)

# S4 method for class 'SubsetSingleCellExperiment,character'
altExpNames(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,character'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,missing'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,numeric'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment'
altExps(x, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
altExp(x, e, ...)

# S4 method for class 'SubsetSpatialExperiment,missing'
altExp(x, e, ...)

# S4 method for class 'SubsetSpatialExperiment,numeric'
altExp(x, e, ...)

# S4 method for class 'SubsetSpatialExperiment'
altExps(x, ...)

# S4 method for class 'SubsetSpatialExperiment'
altExpNames(x, ...)

# S4 method for class 'SubsetSpatialExperiment,character'
altExpNames(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,missing'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,numeric'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetSpatialExperiment'
altExps(x, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
altExp(x, e, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,missing'
altExp(x, e, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,numeric'
altExp(x, e, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment'
altExps(x, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment'
altExpNames(x, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
altExpNames(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,missing'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,numeric'
altExp(x, e, withDimnames = TRUE, withColData = FALSE, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
altExps(x, withDimnames = TRUE, withColData = FALSE, ...) <- value
```

## Arguments

- x:

  An `ExperimentSubset` object.

- e:

  `character(1)` or `integer(1)` The name or position of the alternative
  experiment. Leave it out for the first one. With `subsetName`,
  positions follow `altExpNames()`.

- ...:

  `subsetName`

  : `character(1)` Name of the subset whose alternative experiments to
    get or set. Leave it out for the main object's.

  `unlink`

  : For `altExp<-` and `altExps<-` with `subsetName`: `TRUE`, `FALSE` or
    a character vector of `"altExpName/assayName"` names. Store assays
    of `value` that have the name of a linked altExp assay but different
    values: with `TRUE`, every such assay; with names, only those. Any
    other such assay is an error, and nothing is stored. A message names
    each assay that is unlinked. Default `FALSE`. Any other value is an
    error without `subsetName`.

  Other arguments are passed on to the `SingleCellExperiment` method.
  With `subsetName`, `altExpNames()` and `altExpNames<-` ignore other
  arguments.

- value:

  For `altExp<-`, a `SummarizedExperiment` with one column per column of
  the subset (or of the main object). For `altExps<-`, a list of them.
  For `altExpNames<-`, a `character` vector of new names, one for each
  name `altExpNames()` returns; with `subsetName`, they must all differ,
  and only the subset's own altExps can get new names.

- withDimnames:

  `logical(1)` Passed on to `altExp<-` or `altExps<-`. Default `TRUE`.

- withColData:

  `logical(1)` Passed on to `altExp<-` or `altExps<-`. Default `FALSE`.

## Value

`altExp()` returns a `SummarizedExperiment`, `altExps()` a list of them
and `altExpNames()` a `character` vector. With `subsetName`, these are
every altExp the subset sees: the main object's first, then those stored
only on subsets. The setters return `x` with the change made.

## Details

A subset sees the main object's altExps, and those stored on the subsets
it was made from, cut to its own columns (cells). Their rows (features,
such as the proteins of an ADT altExp) are not cut: a subset's `rows`
pick features of the main object only. These altExps are linked, as
assays are: the subset sees their current values, and stores only what
it changes.

Storing an altExp in a subset with `altExp<-` or `altExps<-` compares it
with the linked one of the same name, assay by assay:

- A new assay is stored as the subset's own, and the linked assays stay
  linked.

- An assay with the values of the linked one is not stored.

- An assay with other values than the linked one needs `unlink`, named
  as `"altExpName/assayName"`, such as `"ADT/counts"`.

- rowData, colData, metadata and reducedDims are stored when they differ
  from the linked altExp's, and then replace them whole. An empty one
  shows the linked one, so a subset can't clear them. Other parts, such
  as the altExp's own altExps, `colPairs` and `rowPairs`, always come
  from the linked altExp.

- An assay the subset already stores itself, after an unlink, is not
  linked, so storing it again needs no `unlink`.

- An altExp with other rows, rownames or class than the linked one
  replaces it entirely, so every linked assay counts as changed. It
  keeps replacing it for as long as the subset stores it, even if the
  linked altExp later has the same rows, and so does an altExp a subset
  stores while nothing of that name is linked. To add to the linked
  altExp instead, remove the subset's with
  `altExp(x, e, subsetName = s) <- NULL` and store it again.

A subset can't remove or rename a linked altExp, or give one of its own
the name of another altExp it sees: an altExp missing from `value` in
`altExps<-` stays visible, and `altExp(x, e, subsetName = s) <- NULL`
removes only what the subset stores, so the linked altExp shows again.

Changes to a linked altExp carry along to the subsets below it, made
from the main object or from the subset that stores it. Its new values
show in them. Renaming it renames what they store of it. Removing it, or
storing it with other rows or another class, drops what they added to
it, with a warning; an altExp a subset stores whole is kept.

## See also

[`altExps`](https://rdrr.io/pkg/SingleCellExperiment/man/altExps.html)

## Examples

``` r
# an object with 10 genes and 20 cells
set.seed(1)
counts <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts) <- paste0("gene", 1:10)
colnames(counts) <- paste0("cell", 1:20)
es <- ExperimentSubset(list(counts = counts))
# counts of 2 proteins in the same 20 cells, as an altExp
adtCounts <- matrix(rpois(40, lambda = 20), nrow = 2)
rownames(adtCounts) <- c("CD3", "CD4")
altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts))
#> Warning: 'colnames(value)' are not the same as 'colnames(x)' for 'altExp<-'.
#> This will be an error in the next release of Bioconductor.
# a subset of the first 5 genes and the first 8 cells
es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
# the subset sees ADT for its 8 cells, with both proteins
altExpNames(es, subsetName = "first")
#> [1] "ADT"
adt <- altExp(es, "ADT", subsetName = "first")
dim(adt)
#> [1] 2 8
# add an assay to ADT in the subset only; ADT's counts stay linked
assay(adt, "logcounts") <- log1p(assay(adt, "counts"))
altExp(es, "ADT", subsetName = "first") <- adt
assayNames(altExp(es, "ADT", subsetName = "first"))
#> [1] "counts"    "logcounts"
assayNames(altExp(es, "ADT"))
#> [1] "counts"
# a changed linked assay needs unlink
assay(adt, "counts") <- assay(adt, "counts") + 1L
altExp(es, "ADT", subsetName = "first", unlink = "ADT/counts") <- adt
#> Subset 'first' now stores its own 'ADT/counts', unlinked from the main altExp 'ADT'.
```
