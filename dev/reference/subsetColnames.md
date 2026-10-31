# subsetColnames

Get `colnames` from a subset in the `ExperimentSubset` object.

Set `colnames` to a subset in the `ExperimentSubset` object.

## Usage

``` r
subsetColnames(x, subsetName)

subsetColnames(x, subsetName) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetColnames(x, subsetName)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetColnames(x, subsetName) <- value

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetColnames(x, subsetName)

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetColnames(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetColnames(x, subsetName)

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetColnames(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
subsetColnames(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetColnames(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetColnames(x, subsetName)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetColnames(x, subsetName) <- value
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  Name of the subset to set `colnames` to.

- value:

  Specify the colname values to replace.

## Value

A `vector` of `colnames`.

Input object with colnames set to a subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetColnames(es, "subset1")
#>  [1] "TTAGTTCTCGTAGGAG" "GGCTGGTGTCTCGTTC" "TCGCGTTGTTAAGTAG" "GGATGTTGTAGGCATG"
#>  [5] "AATCGGTTCTGATACG" "GTTCTCGGTATGAAAC" "TCCACACTCTTTAGTC" "GCGCAGTAGATGCGAC"
#>  [9] "CGGACACTCAAGAAGT" "TGACTTTGTCGCGAAA"
```
