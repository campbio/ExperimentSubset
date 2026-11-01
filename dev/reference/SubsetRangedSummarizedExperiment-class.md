# An S4 class for `RangedSummarizedExperiment` objects with added support for subsets.

An S4 class for `RangedSummarizedExperiment` objects with added support
for subsets.

## Usage

``` r
# S4 method for class 'SubsetRangedSummarizedExperiment'
show(object)
```

## Arguments

- object:

  The object to print, with its subsets and their assays.

## Value

`show()` prints the object and returns `NULL` invisibly.

## Slots

- `subsets`:

  A list of `AssaySubset` objects.
