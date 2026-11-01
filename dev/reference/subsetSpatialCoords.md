# Accessor method for spatialCoords from subsets in ExperimentSubset objects

Get `spatialCoords` from a subset, with
`subsetSpatialCoords(x, subsetName)` or
`spatialCoords(x, subsetName = )`. Without `subsetName`,
`spatialCoords()` returns the main object's, as for `SpatialExperiment`.

## Usage

``` r
subsetSpatialCoords(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetSpatialCoords(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment'
spatialCoords(x, ...)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `spatialCoords` from.

- ...:

  For `spatialCoords()`:

  `subsetName`

  : `character(1)` Name of the subset to get `spatialCoords` from. Leave
    it out for the main object's.

  Without `subsetName`, all arguments are passed on to the
  `SpatialExperiment` method.

## Value

The `spatialCoords` from input object.

## Examples

``` r
# a small SpatialExperiment with six spots
counts <- matrix(1:60, nrow = 10)
coords <- cbind(x = 1:6, y = 6:1)
spe <- SpatialExperiment(assays = list(counts = counts),
                         spatialCoords = coords)
es <- ExperimentSubset(spe)
es <- createSubset(es, "subset1", cols = 2:4)
subsetSpatialCoords(es, "subset1")
#>      x y
#> [1,] 2 5
#> [2,] 3 4
#> [3,] 4 3
# the same, through spatialCoords()
spatialCoords(es, subsetName = "subset1")
#>      x y
#> [1,] 2 5
#> [2,] 3 4
#> [3,] 4 3
```
