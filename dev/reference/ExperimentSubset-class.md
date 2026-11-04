# The class of every ExperimentSubset object

A virtual class that the five `Subset*` classes extend:
[`SubsetSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSummarizedExperiment-class.md),
[`SubsetRangedSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetRangedSummarizedExperiment-class.md),
[`SubsetSingleCellExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSingleCellExperiment-class.md),
[`SubsetTreeSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetTreeSummarizedExperiment-class.md)
and
[`SubsetSpatialExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSpatialExperiment-class.md).
Use `is(x, "ExperimentSubset")` to test whether `x` is an
`ExperimentSubset` object, whichever class it extends. A
`SubsetSpatialExperiment`, for example, is not a
`SubsetSingleCellExperiment`, though it extends `SingleCellExperiment`.

## Details

The class is virtual, so there are no objects of the class itself.
Create an `ExperimentSubset` object with
[`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md).

## Examples

``` r
# simulate counts for 10 genes in 20 spots
set.seed(1)
values <- rpois(200, lambda = 5)
counts <- matrix(values, nrow = 10)
rownames(counts) <- paste0("gene", 1:10)
colnames(counts) <- paste0("spot", 1:20)

# make a SpatialExperiment from the counts
spe <- SpatialExperiment(assays = list(counts = counts))

# convert it to an ExperimentSubset object
es <- ExperimentSubset(spe)

# es is an ExperimentSubset: TRUE
is(es, "ExperimentSubset")
#> [1] TRUE

# but not a SubsetSingleCellExperiment, though it extends
# SingleCellExperiment: FALSE
is(es, "SubsetSingleCellExperiment")
#> [1] FALSE

# the SpatialExperiment it was made from is not an ExperimentSubset: FALSE
is(spe, "ExperimentSubset")
#> [1] FALSE
```
