# Accessor method for colData from subsets in ExperimentSubset objects

Get `colData` from a subset.

Set `colData` to a subset.

## Usage

``` r
subsetColData(x, subsetName, parentColData)

subsetColData(x, subsetName) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetColData(x, subsetName, parentColData)

# S4 method for class 'SubsetRangedSummarizedExperiment,character,DataFrame'
subsetColData(x, subsetName) <- value

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetColData(x, subsetName, parentColData)

# S4 method for class 'SubsetSingleCellExperiment,character,DataFrame'
subsetColData(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment,character,DataFrame'
subsetColData(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetColData(x, subsetName, parentColData)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetColData(x, subsetName, parentColData)

# S4 method for class 'SubsetSpatialExperiment,character,DataFrame'
subsetColData(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetColData(x, subsetName, parentColData)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to set `colData` to.

- parentColData:

  `logical(1)` Logical value indicating if parent colData should be
  combined or not. A subset's own column hides the main object's column
  of the same name, with a message if their values differ. Default
  `FALSE`.

- value:

  Input `DataFrame` to store. Its rownames must be `NULL` or the
  subset's column names, in the same order; other rownames are an error.
  Use `subsetColnames<-` to rename columns.

## Value

The `colData` from input object.

Input object with `colData` stored.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetColData(es, "subset1")
#> DataFrame with 10 rows and 0 columns
```
