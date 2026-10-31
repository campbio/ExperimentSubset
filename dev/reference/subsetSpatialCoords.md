# Accessor method for spatialCoords from subsets in ExperimentSubset objects

Get `spatialCoords` from a subset.

## Usage

``` r
subsetSpatialCoords(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetSpatialCoords(x, subsetName)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `spatialCoords` from.

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
```
