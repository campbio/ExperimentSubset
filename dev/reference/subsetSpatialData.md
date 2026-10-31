# Accessor method for spatialData from subsets in ExperimentSubset objects

Get `spatialData` from a subset.

Set `spatialData` to a subset.

## Usage

``` r
subsetSpatialData(x, subsetName, parentSpatialData)

subsetSpatialData(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
subsetSpatialData(x, subsetName, parentSpatialData)

# S4 method for class 'SubsetSpatialExperiment,character,data.frame'
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

  Input `data.frame` to store.

## Value

The `spatialData` from input object.

Input object with `spatialData` stored.
