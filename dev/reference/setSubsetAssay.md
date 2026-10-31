# Method for storing new assays inside subsets in ExperimentSubset objects

Store a new subset `assay` inside a specified subset in the input
`ExperimentSubset` object.

## Usage

``` r
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)

# S4 method for class 'SubsetRangedSummarizedExperiment,character,ANY,character'
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)

# S4 method for class 'SubsetSingleCellExperiment,character,ANY,character'
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)

# S4 method for class 'SubsetSummarizedExperiment,character,ANY,character'
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)

# S4 method for class 'SubsetSpatialExperiment,character,ANY,character'
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)

# S4 method for class 'SubsetTreeSummarizedExperiment,character,ANY,character'
setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName, unlink = FALSE)
```

## Arguments

- x:

  `ExperimentSubset` Specify the input object.

- subsetName:

  `character(1)` Specify the name of the existing subset inside which
  the new subset `assay` should be stored.

- inputMatrix:

  A matrix-like object (for example a `matrix` or `dgCMatrix`) to store
  as an `assay` in the subset.

- subsetAssayName:

  `character(1)` Specify the name of the new `assay` against the
  `inputMatrix` parameter. It can't contain `"/"`. A subset sees its
  parent's assays as linked assays, so a name linked from the parent
  (for example a main `assay`) needs `unlink = TRUE`.

- unlink:

  `logical(1)` Store the `assay` even though the subset sees an assay of
  that name linked from its parent. The subset, and subsets made from
  it, then use this `assay` instead; the parent's is unchanged. Default
  `FALSE`.

## Value

Updated `ExperimentSubset` object with the new `assay` stored inside the
specified subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es, "subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
counts1p <- assay(es, "subset1")
counts1p[,] <- counts1p[,] + 1
es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
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
#> subsetAssays(2): subset1 subset1/scaledSubset1
```
