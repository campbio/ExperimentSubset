# Accessor method for rowLinks from subsets in ExperimentSubset objects

Get `rowLinks` from a subset.

## Usage

``` r
subsetRowLinks(x, subsetName, parentRowLinkData)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetRowLinks(x, subsetName, parentRowLinkData)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `rowLinks` from.

- parentRowLinkData:

  `logical(1)` Logical value indicating if parent rowLinks should be
  combined or not. Default `FALSE`.

## Value

The `rowLinks` from input object.

## Examples

``` r
# a small TreeSummarizedExperiment with a tree on both rows and columns
data("tinyTree", package = "TreeSummarizedExperiment")
counts <- matrix(1:100, nrow = 10)
# name rows and columns after the tips of the tree
rownames(counts) <- tinyTree$tip.label
colnames(counts) <- tinyTree$tip.label
tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                rowTree = tinyTree, colTree = tinyTree)
es <- ExperimentSubset(tse)
es <- createSubset(es, "subset1", rows = 1:5, cols = 1:5)
subsetRowLinks(es, "subset1", parentRowLinkData = TRUE)
#> LinkDataFrame with 5 rows and 5 columns
#>        nodeLab nodeLab_alias   nodeNum    isLeaf   whichTree
#>    <character>   <character> <integer> <logical> <character>
#> t2          t2       alias_1         1      TRUE       phylo
#> t7          t7       alias_2         2      TRUE       phylo
#> t6          t6       alias_3         3      TRUE       phylo
#> t9          t9       alias_4         4      TRUE       phylo
#> t4          t4       alias_5         5      TRUE       phylo
```
