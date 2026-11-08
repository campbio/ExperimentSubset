# Accessor method for imgData of subsets in ExperimentSubset objects

Get the images of a subset: the rows of the main object's `imgData`
whose `sample_id` is the sample of one or more of the subset's spots.
They are looked up from the main object each time, so they are always
current; subsets don't store images of their own.

## Usage

``` r
subsetImgData(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment'
subsetImgData(x, subsetName)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `imgData` from.

## Value

A `DataFrame` of the subset's images, with the columns of `imgData()`.

## Examples

``` r
# a small SpatialExperiment with six spots in two samples, and an image
# of each sample
counts <- matrix(1:60, nrow = 10)
coords <- cbind(x = 1:6, y = 6:1)
sampleIds <- rep(c("A", "B"), each = 3)
image <- SpatialImage(as.raster(matrix(c(0, 1, 1, 0), nrow = 2)))
images <- DataFrame(sample_id = c("A", "B"), image_id = "lowres",
                    data = I(list(image, image)), scaleFactor = 1)
spe <- SpatialExperiment(assays = list(counts = counts),
                         sample_id = sampleIds,
                         spatialCoords = coords,
                         imgData = images)
es <- ExperimentSubset(spe)
# a subset of spots from sample B only
es <- createSubset(es, "sampleB", cols = 4:6)
subsetImgData(es, subsetName = "sampleB")
#> DataFrame with 1 row and 4 columns
#>     sample_id    image_id   data scaleFactor
#>   <character> <character> <list>   <numeric>
#> 1           B      lowres   ####           1
```
