# reducedDims\<-

A wrapper to the `reducedDims<-` from
[reducedDims](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html)
method with additional support for subsets.

## Usage

``` r
# S4 method for class 'ANY'
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
