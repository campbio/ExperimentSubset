# Get subset assay from an `ExperimentSubset` object.

Get subset assay from an `ExperimentSubset` object.

## Usage

``` r
getSubsetAssay(x, subsetName, subsetAssayName = NULL)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
getSubsetAssay(x, subsetName, subsetAssayName = NULL)

# S4 method for class 'SubsetSingleCellExperiment,character'
getSubsetAssay(x, subsetName, subsetAssayName = NULL)

# S4 method for class 'SubsetSummarizedExperiment,character'
getSubsetAssay(x, subsetName, subsetAssayName = NULL)

# S4 method for class 'SubsetSpatialExperiment,character'
getSubsetAssay(x, subsetName, subsetAssayName = NULL)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
getSubsetAssay(x, subsetName, subsetAssayName = NULL)
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  The name of the subset to fetch the assay from, or an assay seen in a
  subset written as `"subset/assay"`.

- subsetAssayName:

  `character(1)` The name of an assay seen in the subset `subsetName`:
  one stored with `setSubsetAssay`, or one linked from its parent.
  Default `NULL`, which returns the subset's slice of its `parentAssay`.

## Value

Subset assay

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
getSubsetAssay(es, "subset1")
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 10 column names ‘TTAGTTCTCGTAGGAG’, ‘GGCTGGTGTCTCGTTC’, ‘TCGCGTTGTTAAGTAG’ ... ]]
#>                                  
#> OPCML         . . . . . . . . . .
#> RP11-338K17.8 . . . . . . . . . .
#> RP11-770J1.3  . . . . . . . . . .
#> NPSR1         . . . . . . . . . .
#> AC060835.1    . . . . . . . . . .
#> GABRA2        . . . . . . . . . .
#> PRDM16        . . . . . . . . . .
#> CAB39         . . . . . . . . . .
#> LYRM1         . . 1 2 . . . . . 1
#> RP1-178F10.3  . . . . . . . . . .
# store a changed copy of the subset's counts inside the subset
subset1Counts <- getSubsetAssay(es, "subset1")
es <- setSubsetAssay(es, "subset1", subset1Counts + 1, "plusOne")
getSubsetAssay(es, "subset1", subsetAssayName = "plusOne")
#> 10 x 10 Matrix of class "dgeMatrix"
#>               TTAGTTCTCGTAGGAG GGCTGGTGTCTCGTTC TCGCGTTGTTAAGTAG
#> OPCML                        1                1                1
#> RP11-338K17.8                1                1                1
#> RP11-770J1.3                 1                1                1
#> NPSR1                        1                1                1
#> AC060835.1                   1                1                1
#> GABRA2                       1                1                1
#> PRDM16                       1                1                1
#> CAB39                        1                1                1
#> LYRM1                        1                1                2
#> RP1-178F10.3                 1                1                1
#>               GGATGTTGTAGGCATG AATCGGTTCTGATACG GTTCTCGGTATGAAAC
#> OPCML                        1                1                1
#> RP11-338K17.8                1                1                1
#> RP11-770J1.3                 1                1                1
#> NPSR1                        1                1                1
#> AC060835.1                   1                1                1
#> GABRA2                       1                1                1
#> PRDM16                       1                1                1
#> CAB39                        1                1                1
#> LYRM1                        3                1                1
#> RP1-178F10.3                 1                1                1
#>               TCCACACTCTTTAGTC GCGCAGTAGATGCGAC CGGACACTCAAGAAGT
#> OPCML                        1                1                1
#> RP11-338K17.8                1                1                1
#> RP11-770J1.3                 1                1                1
#> NPSR1                        1                1                1
#> AC060835.1                   1                1                1
#> GABRA2                       1                1                1
#> PRDM16                       1                1                1
#> CAB39                        1                1                1
#> LYRM1                        1                1                1
#> RP1-178F10.3                 1                1                1
#>               TGACTTTGTCGCGAAA
#> OPCML                        1
#> RP11-338K17.8                1
#> RP11-770J1.3                 1
#> NPSR1                        1
#> AC060835.1                   1
#> GABRA2                       1
#> PRDM16                       1
#> CAB39                        1
#> LYRM1                        2
#> RP1-178F10.3                 1
```
