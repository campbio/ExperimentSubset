# Method for displaying 'child-parent' link structure of subsets in ExperimentSubset objects

The function displays the content of an `ExperimentSubset` object
including all available main assays, all subsets and the subset assays
inside these subsets. This function also depicts how and in what order
the subsets in the object are linked with their parents. Moreover, all
supplementary data inside the subsets such as `reducedDims` and
`altExps` are also displayed against each subset entry.

## Usage

``` r
subsetSummary(x)

# S4 method for class 'SubsetRangedSummarizedExperiment'
subsetSummary(x)

# S4 method for class 'SubsetSingleCellExperiment'
subsetSummary(x)

# S4 method for class 'SubsetSummarizedExperiment'
subsetSummary(x)

# S4 method for class 'SubsetSpatialExperiment'
subsetSummary(x)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetSummary(x)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

## Value

Prints all the available subset information against the input
`ExperimentSubset` object.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
assay(es, "subset1",
subsetAssayName = "subset1pAssay") <- assay(es, "subset1")[,] + 1
#> Storing assay inside subset 'subset1'.
subsetSummary(es)
#> Main assay(s):
#>  counts logcounts 
#> 
#> Subset(s):
#>      Name    Dim Parent        Assays
#> 1 subset1 10, 10 counts subset1pAssay
```
