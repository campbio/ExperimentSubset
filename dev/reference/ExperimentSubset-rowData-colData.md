# rowData and colData of subsets in ExperimentSubset objects

Get or set the `rowData` or `colData` of a subset by naming it with
`subsetName`. Without `subsetName`, these methods act on the main
object, as for the class `x` extends.

A subset stores only its own columns. With `parentRowData = TRUE` or
`parentColData = TRUE`, the main object's columns come first, for the
subset's rows or columns, and a subset's own column hides the main
object's column of the same name, with a message if their values differ.
When setting a subset's data, columns that came from the main object, in
data read with `parentRowData = TRUE` or `parentColData = TRUE`, are
left out with a message, so only the subset's own are stored. To add one
column, use `colData(x, subsetName = "name")$column <- values`.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'SubsetRangedSummarizedExperiment'
colData(x, ...)

# S4 method for class 'SubsetRangedSummarizedExperiment'
rowData(x, ...) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,data.frame'
colData(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'SubsetSingleCellExperiment'
colData(x, ...)

# S4 method for class 'SubsetSingleCellExperiment'
rowData(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SubsetSingleCellExperiment,data.frame'
colData(x, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'SubsetSummarizedExperiment'
colData(x, ...)

# S4 method for class 'SubsetSummarizedExperiment'
rowData(x, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SubsetSummarizedExperiment,data.frame'
colData(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'SubsetSpatialExperiment'
colData(x, ...)

# S4 method for class 'SubsetSpatialExperiment'
rowData(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SubsetSpatialExperiment,data.frame'
colData(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
rowData(x, use.names = TRUE, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment'
colData(x, ...)

# S4 method for class 'SubsetTreeSummarizedExperiment'
rowData(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,DataFrame'
colData(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,NULL'
colData(x, ...) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,data.frame'
colData(x, ...) <- value
```

## Arguments

- x:

  An `ExperimentSubset` object.

- use.names:

  `logical(1)` Whether the `rowData` returned has row names. Default
  `TRUE`.

- ...:

  `subsetName`

  : `character(1)` Name of the subset whose data to get or set. Leave it
    out for the main object's.

  `parentRowData`

  : `logical(1)` For `rowData()` of a subset, whether to include the
    main object's columns. Default `FALSE`.

  `parentColData`

  : `logical(1)` For `colData()` of a subset, whether to include the
    main object's columns. Default `FALSE`.

  Without `subsetName`, other arguments are passed on to the method of
  the class `x` extends.

- value:

  A `DataFrame` with one row per row (`rowData<-`) or column
  (`colData<-`) of the subset or the main object, or `NULL` to remove
  all columns (a `SpatialExperiment` keeps its `sample_id` column). For
  a subset, a `data.frame` is also accepted and converted to a
  `DataFrame`; for the main object, `colData<-` with a `data.frame` is
  an error.

## Value

`rowData()` and `colData()` return a `DataFrame`. The setters return `x`
with the data set.

## See also

[`subsetRowData`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRowData.md),
[`subsetColData`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColData.md),
[`rowData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)

## Examples

``` r
# an object with 10 genes and 20 cells
set.seed(1)
counts <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts) <- paste0("gene", 1:10)
colnames(counts) <- paste0("cell", 1:20)
es <- ExperimentSubset(list(counts = counts))
# a sample label for every cell of the main object
es$sample <- rep(c("A", "B"), each = 10)
# a subset of the first 5 genes and the first 8 cells
es <- createSubset(es, subsetName = "first", rows = 1:5, cols = 1:8)
# a cluster label for each of the subset's 8 cells
colData(es, subsetName = "first")$cluster <- rep(1:2, times = 4)
# the subset's own columns
colData(es, subsetName = "first")
#> DataFrame with 8 rows and 1 column
#>         cluster
#>       <integer>
#> cell1         1
#> cell2         2
#> cell3         1
#> cell4         2
#> cell5         1
#> cell6         2
#> cell7         1
#> cell8         2
# with the main object's columns first
colData(es, subsetName = "first", parentColData = TRUE)
#> DataFrame with 8 rows and 2 columns
#>            sample   cluster
#>       <character> <integer>
#> cell1           A         1
#> cell2           A         2
#> cell3           A         1
#> cell4           A         2
#> cell5           A         1
#> cell6           A         2
#> cell7           A         1
#> cell8           A         2
# the subset's rowData, which has no columns yet
rowData(es, subsetName = "first")
#> DataFrame with 5 rows and 0 columns
```
