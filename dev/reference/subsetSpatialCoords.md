# Get or set spatialCoords of subsets in ExperimentSubset objects

Get `spatialCoords` from a subset, with
`subsetSpatialCoords(x, subsetName)` or
`spatialCoords(x, subsetName = )`. Without `subsetName`,
`spatialCoords()` returns the main object's, as for `SpatialExperiment`.

`subsetSpatialCoords(x, subsetName) <- value` stores a subset's own
`spatialCoords`, for example after a transformation, without changing
the main object's. A subset with its own coordinates returns them from
`subsetSpatialCoords()`, `spatialCoords(x, subsetName = )` and
[`getSubsetExperiment()`](https://campbio.github.io/ExperimentSubset/dev/reference/getSubsetExperiment.md),
and subsets made from it see them, cut to their spots. A subset without
its own coordinates sees those of the nearest subset it was made from
that has them, or else the main object's.

## Usage

``` r
subsetSpatialCoords(x, subsetName)

subsetSpatialCoords(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment'
subsetSpatialCoords(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment'
subsetSpatialCoords(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment'
spatialCoords(x, ...)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `spatialCoords` from, or to
  store them in.

- value:

  A numeric `matrix` with one row for each spot in the subset and one or
  more columns. Its rownames must be `NULL` or the subset's column
  names, in the same order. `NULL` removes the subset's own coordinates,
  so it sees those of the subset it was made from, or the main object's,
  again.

- ...:

  For `spatialCoords()`:

  `subsetName`

  : `character(1)` Name of the subset to get `spatialCoords` from. Leave
    it out for the main object's.

  Without `subsetName`, all arguments are passed on to the
  `SpatialExperiment` method.

## Value

For `subsetSpatialCoords()` and `spatialCoords()`, the `spatialCoords`
the subset sees.

For `subsetSpatialCoords<-`, the input object with the subset's
`spatialCoords` stored.

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
# a small SpatialExperiment with six spots
counts <- matrix(1:60, nrow = 10)
coords <- cbind(x = 1:6, y = 6:1)
spe <- SpatialExperiment(assays = list(counts = counts),
                         spatialCoords = coords)
es <- ExperimentSubset(spe)
es <- createSubset(es, "subset1", cols = 2:4)
# store coordinates moved 100 to the right
movedCoords <- subsetSpatialCoords(es, subsetName = "subset1")
movedCoords[, "x"] <- movedCoords[, "x"] + 100
subsetSpatialCoords(es, subsetName = "subset1") <- movedCoords
subsetSpatialCoords(es, subsetName = "subset1")
#>        x y
#> [1,] 102 5
#> [2,] 103 4
#> [3,] 104 3
# the main object's coordinates are unchanged
spatialCoords(es)
#>      x y
#> [1,] 1 6
#> [2,] 2 5
#> [3,] 3 4
#> [4,] 4 3
#> [5,] 5 2
#> [6,] 6 1
```
