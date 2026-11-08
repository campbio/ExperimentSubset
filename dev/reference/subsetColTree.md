# Get or set the colTree of subsets in ExperimentSubset objects

Get the column tree of a subset: its own, if it stores one; otherwise
that of the nearest subset it was made from that stores one; otherwise
the main object's. TreeSummarizedExperiment's `colTree` takes no
`subsetName`, so subsets use `subsetColTree`.

`subsetColTree<-` stores a tree of the subset's own, linking each column
to the node whose label is the column's name in the subset, without
changing the main object or other subsets. Subsets made from it use the
tree too, unless they store their own. `NULL` removes the subset's own
tree. `setSubsetExperiment` also stores a column tree that a tool
changed, but not the removal of one, as a subset always uses its own
tree, a parent's or the main object's. A subset's own tree is its only
column tree, so get it with `whichTree = 1`. `rbind` and `cbind` keep
subsets' own trees; `combineCols` and `combineRows` return classes
without trees (see
[ExperimentSubset-combine](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)).

## Usage

``` r
subsetColTree(x, subsetName, whichTree = 1)

subsetColTree(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetColTree(x, subsetName, whichTree = 1)

# S4 method for class 'SubsetTreeSummarizedExperiment'
subsetColTree(x, subsetName) <- value
```

## Arguments

- x:

  `ExperimentSubset` Input `ExperimentSubset` object.

- subsetName:

  `character(1)` Name of the subset to get the tree of, or to store it
  in.

- whichTree:

  The name or position of the tree to get, as in `colTree`. Default `1`.

- value:

  A `phylo` tree whose node labels include every column name of the
  subset, or `NULL`.

## Value

For `subsetColTree`, a `phylo` tree, or `NULL` if there is none. For
`subsetColTree<-`, the input object with the subset's tree stored.

## Examples

``` r
# a small TreeSummarizedExperiment with a tree on its columns
data("tinyTree", package = "TreeSummarizedExperiment")
counts <- matrix(1:100, nrow = 10)
colnames(counts) <- tinyTree$tip.label
tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                colTree = tinyTree)
es <- ExperimentSubset(tse)
es <- createSubset(es, "subset1", cols = 1:5)
# store a tree with the same labels on other tips
ownTree <- tinyTree
ownTree$tip.label <- rev(tinyTree$tip.label)
subsetColTree(es, subsetName = "subset1") <- ownTree
subsetColTree(es, subsetName = "subset1")
#> 
#> Phylogenetic tree with 10 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   t3, t5, t1, t10, t8, t4, ...
#> Node labels:
#>   Node_11, Node_12, Node_13, Node_14, Node_15, Node_16, ...
#> 
#> Rooted; includes branch length(s).
subsetColLinks(es, subsetName = "subset1")
#> LinkDataFrame with 5 rows and 5 columns
#>        nodeLab   nodeNum nodeLab_alias    isLeaf   whichTree
#>    <character> <integer>   <character> <logical> <character>
#> t2          t2        10      alias_10      TRUE       phylo
#> t7          t7         9       alias_9      TRUE       phylo
#> t6          t6         8       alias_8      TRUE       phylo
#> t9          t9         7       alias_7      TRUE       phylo
#> t4          t4         6       alias_6      TRUE       phylo
```
