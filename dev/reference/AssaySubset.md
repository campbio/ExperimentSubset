# AssaySubset constructor

Constructor for creating a experiment object internally by the
`ExperimentSubset` object. Should not be used directly by the user.

## Usage

``` r
AssaySubset(
  subsetName = "subset",
  rowIndices = NULL,
  colIndices = NULL,
  parentAssay = "counts",
  internalAssay = NULL
)
```

## Arguments

- subsetName:

  `character(1)` Name of the subset.

- rowIndices:

  `vector("numeric")` Indices of the rows to include in the subset.

- colIndices:

  `vector("numeric")` Indices of the columns to include in the subset.

- parentAssay:

  `character(1)` Name of the parent of this subset.

- internalAssay:

  An internal Experiment object to store additional subset data.

## Value

A `AssaySubset` object.
