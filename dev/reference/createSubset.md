# Subset creation method for ExperimentSubset objects

Create a subset from an already available `assay` in the input
`ExperimentSubset` object by specifying the rows and columns to include
in the subset.

## Usage

``` r
createSubset(
  x,
  subsetName,
  rows = NULL,
  cols = NULL,
  parentAssay = NULL,
  altExp = NULL
)

# S4 method for class 'SubsetSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL, altExp = NULL)

# S4 method for class 'SubsetRangedSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL, altExp = NULL)

# S4 method for class 'SubsetSingleCellExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL, altExp = NULL)

# S4 method for class 'SubsetTreeSummarizedExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL, altExp = NULL)

# S4 method for class 'SubsetSpatialExperiment'
createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay = NULL, altExp = NULL)
```

## Arguments

- x:

  `ExperimentSubset` Specify the object from which a subset should be
  created. Use
  [`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
  first to convert another experiment object.

- subsetName:

  `character(1)` Specify the name of the subset to create.

- rows:

  `numeric` or `character` vector. Specify the rows to include in this
  subset. If `missing` or `NULL`, all rows are included in the subset.
  Default `NULL`.

- cols:

  `numeric` or `character` vector. Specify the columns to include in
  this subset. If `missing` or `NULL`, all columns are included in the
  subset. Default `NULL`.

- parentAssay:

  `character(1)` Specify the parent `assay` of the subset: the name of a
  main `assay`, of a subset, or of an assay seen in a subset written as
  `"subset/assay"`. If `NULL`, the first available main `assay` will be
  marked as parent. Default `NULL`. With `altExp`, it names an assay of
  that altExp instead; see the section below.

- altExp:

  `character(1)` or `NULL`. The name of an alternative experiment
  (altExp) of the main object, such as `"ADT"`, to make the subset from:
  `rows` then pick its features. Only for the classes that extend
  `SingleCellExperiment`. If `NULL`, the subset is made from the main
  object's features. Default `NULL`.

## Value

An `ExperimentSubset` object that now contains the newly created subset.

## Subsets made from an altExp

With `altExp = "ADT"`, the subset's rows are features of the altExp
`"ADT"`, such as some of its proteins, and its columns are cells, as for
any subset. `parentAssay` names an assay of the altExp, as the main
object or a subset sees it:

- `NULL` or `"counts"`: the first assay, or the assay `counts`, of the
  main object's altExp. `cols` pick the main object's cells.

- `"clean"` or `"clean/clr"`: the first assay, or the assay `clr`, of
  the altExp as subset `clean` sees it, with the assays `clean` stores
  in it. `cols` pick `clean`'s cells.

`rows` pick the altExp's features, by position or name. The subset sees
every assay of the altExp as its parent does, linked and cut to its rows
and cells, and stores its own with
[`setSubsetAssay()`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetAssay.md),
with the usual `unlink` rules. Its rowData with `parentRowData = TRUE`
is the altExp's, and
[`getSubsetExperiment()`](https://campbio.github.io/ExperimentSubset/dev/reference/getSubsetExperiment.md)
returns an object of the altExp's class. Subsets made from it, without
`altExp`, pick its rows and cells. It doesn't see the main object's
altExps. A row tree doesn't apply to it; column trees and spatial data
are read at its cells.

Changes to the altExp carry along: renaming it, on the main object or in
the subset it was made from, renames what the subset records; if it gets
other rows, the subset keeps the features whose names remain; if it is
removed, the subset and those made from it are dropped, with a warning.
If the assay the subset was made from goes, it behaves as for a main
assay: storing a subset's view of the altExp without that assay is an
error, and a change to the main object's altExp that drops it leaves the
subset without its parent assay.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
es
#> class: SubsetSingleCellExperiment 
#> dim: 2000 2000 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3
#> rowData names(1): gene
#> colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT
#>   ATTTCTGGTGATGATA
#> colData names(11): nGene nUMI ... ident doublet_true_labels
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> subsets(1): subset1
#> subsetAssays(1): subset1

# an object with 10 genes and 20 cells
set.seed(1)
counts <- matrix(rpois(200, lambda = 5), nrow = 10)
rownames(counts) <- paste0("gene", 1:10)
colnames(counts) <- paste0("cell", 1:20)
es <- ExperimentSubset(list(counts = counts))
# counts of 4 proteins in the same cells, as the altExp "ADT"
adtCounts <- matrix(rpois(80, lambda = 20), nrow = 4)
rownames(adtCounts) <- c("CD3", "CD4", "CD8", "Isotype")
altExp(es, "ADT") <- SingleCellExperiment(list(counts = adtCounts))
#> Warning: 'colnames(value)' are not the same as 'colnames(x)' for 'altExp<-'.
#> This will be an error in the next release of Bioconductor.
# a subset of the first 12 cells
es <- createSubset(es, subsetName = "clean", cols = 1:12)
# a subset of clean's cells without the isotype control
es <- createSubset(es, subsetName = "clean.antibodies",
                   rows = c("CD3", "CD4", "CD8"),
                   parentAssay = "clean", altExp = "ADT")
dim(assay(es, "clean.antibodies"))
#> [1]  3 12
subsetParent(es, "clean.antibodies")
#> [[1]]
#> [1] "clean/ADT/counts"
#> 
#> [[2]]
#> [1] "clean"
#> 
#> [[3]]
#> [1] "counts"
#> 
```
