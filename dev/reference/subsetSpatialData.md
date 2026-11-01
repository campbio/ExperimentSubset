# Deprecated: spatialData of subsets in ExperimentSubset objects

`subsetSpatialData()` and `spatialData(x, subsetName =)` are deprecated,
and `subsetSpatialData<-` is defunct, because `spatialData()` and
`spatialData<-` are deprecated in SpatialExperiment, which now keeps
these columns in `colData`. Use `subsetColData(x, subsetName)` to get
them, which returns all of the subset's colData columns, with
`parentColData = TRUE` to include the main object's, and
`colData(x, subsetName = "name")$column <- values` to add a column to a
subset. The getters will be defunct in the next release.

## Usage

``` r
subsetSpatialData(x, subsetName, parentSpatialData)

subsetSpatialData(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
subsetSpatialData(x, subsetName, parentSpatialData)

# S4 method for class 'SubsetSpatialExperiment'
subsetSpatialData(x, subsetName) <- value
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to set `spatialData` to.

- parentSpatialData:

  `logical(1)` Logical value indicating if parent spatialData should be
  combined or not. Default `FALSE`.

- value:

  Ignored: the setter is defunct.

## Value

The `spatialData` from input object, with a deprecation warning. The
setter fails with an error.

## See also

[`subsetColData`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColData.md)
