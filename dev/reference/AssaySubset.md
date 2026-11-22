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
  internalAssay = NULL,
  altExpName = NULL
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

- altExpName:

  `character(1)` or `NULL`. The altExp whose features the subset's rows
  index, for a subset made with `createSubset(altExp =)`.

## Value

A `AssaySubset` object.
