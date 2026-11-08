# Accessor method for colLinks from subsets in ExperimentSubset objects

Get `colLinks` from a subset: the links of its columns to the nodes of
the tree that
[`subsetColTree`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColTree.md)
returns.

## Usage

``` r
subsetColLinks(x, subsetName, parentColLinkData = TRUE)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetColLinks(x, subsetName, parentColLinkData = TRUE)
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get `colLinks` from.

- parentColLinkData:

  `logical(1)` If `TRUE` (the default), a subset without a column tree
  of its own, or one it was made from, gets the links to the main
  object's column tree. If `FALSE`, it gets only links to a column tree
  stored by the subset or one it was made from, and `NULL` if there is
  none.

## Value

The `colLinks` of the subset's columns, with the subset's column names,
or `NULL`.

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
subsetColLinks(es, subsetName = "subset1")
#> LinkDataFrame with 5 rows and 5 columns
#>        nodeLab nodeLab_alias   nodeNum    isLeaf   whichTree
#>    <character>   <character> <integer> <logical> <character>
#> t2          t2       alias_1         1      TRUE       phylo
#> t7          t7       alias_2         2      TRUE       phylo
#> t6          t6       alias_3         3      TRUE       phylo
#> t9          t9       alias_4         4      TRUE       phylo
#> t4          t4       alias_5         5      TRUE       phylo
```
