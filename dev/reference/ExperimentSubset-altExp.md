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
  experiment. Leave it out for the first one.

- ...:

  `subsetName`

  : `character(1)` Name of the subset whose alternative experiments to
    get or set. Leave it out for the main object's.

  Other arguments are passed on to the `SingleCellExperiment` method.
  With `subsetName`, `altExpNames()` and `altExpNames<-` ignore other
  arguments.

- value:

  For `altExp<-`, a `SummarizedExperiment` with one column per column of
  the subset (or of the main object). For `altExps<-`, a list of them.
  For `altExpNames<-`, a `character` vector of new names.

- withDimnames:

  `logical(1)` Passed on to `altExp<-` or `altExps<-`. Default `TRUE`.

- withColData:

  `logical(1)` Passed on to `altExp<-` or `altExps<-`. Default `FALSE`.

## Value

`altExp()` returns a `SummarizedExperiment`, `altExps()` a list of them
and `altExpNames()` a `character` vector. The setters return `x` with
the change made.

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
# a subset of the first 5 genes and the first 8 cells
es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
# counts of 2 spike-in genes in the subset's 8 cells
spikeCounts <- matrix(rpois(16, lambda = 5), nrow = 2)
spikes <- SingleCellExperiment(list(counts = spikeCounts))
# store them as an alternative experiment of the subset
altExp(es, "spikes", subsetName = "first") <- spikes
#> Warning: 'colnames(value)' are not the same as 'colnames(x)' for 'altExp<-'.
#> This will be an error in the next release of Bioconductor.
# list the subset's alternative experiments
altExpNames(es, subsetName = "first")
#> [1] "spikes"
# read it back
altExp(es, "spikes", subsetName = "first")
#> class: SingleCellExperiment 
#> dim: 2 8 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames(8): cell1 cell2 ... cell7 cell8
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
```
