# Get or set the rowTree of subsets in ExperimentSubset objects

Get the row tree of a subset: its own, if it stores one; otherwise that
of the nearest subset it was made from that stores one; otherwise the
main object's. TreeSummarizedExperiment's `rowTree` takes no
`subsetName`, so subsets use `subsetRowTree`.

`subsetRowTree<-` stores a tree of the subset's own, linking each row to
the node whose label is the row's name in the subset, without changing
the main object or other subsets. Subsets made from it use the tree too,
unless they store their own. `NULL` removes the subset's own tree.
`setSubsetExperiment` also stores a row tree that a tool changed, but
not the removal of one, as a subset always uses its own tree, a parent's
or the main object's. A subset's own tree is its only row tree, so get
it with `whichTree = 1`. `rbind` and `cbind` keep subsets' own trees;
`combineCols` and `combineRows` return classes without trees (see
[ExperimentSubset-combine](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)).

## Usage

``` r
subsetRowTree(x, subsetName, whichTree = 1)

subsetRowTree(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetRowTree(x, subsetName, whichTree = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetRowTree(x, subsetName) <- value
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get the tree of, or to store it
  in.

- whichTree:

  The name or position of the tree to get, as in `rowTree`. Default `1`.

- value:

  A `phylo` tree whose node labels include every rowname of the subset,
  or `NULL`.

## Value

For `subsetRowTree`, a `phylo` tree, or `NULL` if there is none. For
`subsetRowTree<-`, the input object with the subset's tree stored.

## Examples

``` r
# a small TreeSummarizedExperiment with a tree on its rows
data("tinyTree", package = "TreeSummarizedExperiment")
counts <- matrix(1:100, nrow = 10)
rownames(counts) <- tinyTree$tip.label
tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                rowTree = tinyTree)
es <- ExperimentSubset(tse)
es <- createSubset(es, "subset1", rows = 1:5)
# the subset uses the main object's tree
subsetRowTree(es, subsetName = "subset1")
#> 
#> Phylogenetic tree with 10 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   t2, t7, t6, t9, t4, t8, ...
#> Node labels:
#>   Node_11, Node_12, Node_13, Node_14, Node_15, Node_16, ...
#> 
#> Rooted; includes branch length(s).
# store a tree with the same labels on other tips
ownTree <- tinyTree
ownTree$tip.label <- rev(tinyTree$tip.label)
subsetRowTree(es, subsetName = "subset1") <- ownTree
subsetRowLinks(es, subsetName = "subset1")
#> LinkDataFrame with 5 rows and 5 columns
#>        nodeLab   nodeNum nodeLab_alias    isLeaf   whichTree
#>    <character> <integer>   <character> <logical> <character>
#> t2          t2        10      alias_10      TRUE       phylo
#> t7          t7         9       alias_9      TRUE       phylo
#> t6          t6         8       alias_8      TRUE       phylo
#> t9          t9         7       alias_7      TRUE       phylo
#> t4          t4         6       alias_6      TRUE       phylo
# the main object's links are unchanged
rowLinks(es)[1:5, ]
#> LinkDataFrame with 5 rows and 5 columns
#>        nodeLab nodeLab_alias   nodeNum    isLeaf   whichTree
#>    <character>   <character> <integer> <logical> <character>
#> t2          t2       alias_1         1      TRUE       phylo
#> t7          t7       alias_2         2      TRUE       phylo
#> t6          t6       alias_3         3      TRUE       phylo
#> t9          t9       alias_4         4      TRUE       phylo
#> t4          t4       alias_5         5      TRUE       phylo
```
