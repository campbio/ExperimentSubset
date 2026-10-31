# reducedDims\<-

A wrapper to the `reducedDims<-` from
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
method with additional support for subsets.

## Usage

``` r
reducedDims(x, subsetName) <- value
```

## Arguments

- x:

  Input `ExperimentSubset` object or any object supported by
  `reducedDims<-` from
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
  method.

- subsetName:

  Specify the name of the subset to which the `reducedDims` should be
  set to. If `missing`, `reducedDims<-` from
  [reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
  method is called on the main object.

- value:

  A `list` of values to set to `reducedDims`.

## Value

Updated input object with `reducedDims` set.

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
reducedDims(es, subsetName = "subset1")
#> List of length 2
#> names(2): PCA_1 PCA_2
```
