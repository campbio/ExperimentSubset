# Accessor method for rowData from subsets in ExperimentSubset objects

Get `rowData` from a subset.

Set `rowData` to a subset.

## Usage

``` r
subsetRowData(x, subsetName, parentRowData)

subsetRowData(x, subsetName) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetRowData(x, subsetName, parentRowData)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetRowData(x, subsetName) <- value

# S4 method for class 'SubsetSingleCellExperiment'
subsetRowData(x, subsetName, parentRowData)

# S4 method for class 'SubsetSingleCellExperiment'
subsetRowData(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment'
subsetRowData(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment'
subsetRowData(x, subsetName, parentRowData)

# S4 method for class 'SubsetSpatialExperiment'
subsetRowData(x, subsetName, parentRowData)

# S4 method for class 'SubsetSpatialExperiment'
subsetRowData(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetRowData(x, subsetName, parentRowData)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to set `rowData` to.

- parentRowData:

  `logical(1)` Logical value indicating if parent rowData should be
  combined or not. A subset's own column hides the main object's column
  of the same name, with a message if their values differ. Default
  `FALSE`.

- value:

  Input `DataFrame` to store; a `data.frame` is converted to one, and
  `NULL` removes the subset's own `rowData` columns. Its rownames must
  be `NULL` or the subset's row names, in the same order; other rownames
  are an error, since rows are matched by position. To store a table
  with other rownames, such as gene IDs for rows named by symbol, put
  its rows in the subset's order and set `rownames(value) <- NULL`
  first.

## Value

The `rowData` from input object.

Input object with `rowData` stored.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetRowData(es, "subset1")
#> DataFrame with 10 rows and 0 columns
```
