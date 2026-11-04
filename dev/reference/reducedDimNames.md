# reducedDimNames

A wrapper to the `reducedDimNames` from
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
method with additional support for subsets.

## Usage

``` r
reducedDimNames(x, ...)

# S4 method for class 'ANY'
reducedDimNames(x, ...)
```

## Arguments

- x:

  Input `ExperimentSubset` object or any object supported by
  `reducedDimNames` from
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
  method.

- ...:

  Only `subsetName`, the name of the subset whose `reducedDimNames` to
  get, given by name. Leave it out for the main object's. A subset name
  by position, or any other argument, is an error.

## Value

The `reducedDimNames` from the specified subset or same as
`reducedDimNames` from
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
when `subsetName` is `missing`.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1",
rows = c(1:50), cols = c(1:20),
parentAssay = "counts")
# two small reduced dimensions, one row per cell in the subset
set.seed(1)
pca1 <- matrix(rnorm(40), nrow = 20)
pca2 <- matrix(rnorm(40), nrow = 20)
reducedDims(es, subsetName = "subset1") <- list(PCA_1 = pca1, PCA_2 = pca2)
reducedDimNames(es, subsetName = "subset1")
#> [1] "PCA_1" "PCA_2"
```
