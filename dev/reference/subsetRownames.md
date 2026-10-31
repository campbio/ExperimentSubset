# subsetRownames

Get `rownames` from a subset in the `ExperimentSubset` object.

Set `rownames` to a subset in the `ExperimentSubset` object.

## Usage

``` r
subsetRownames(x, subsetName)

subsetRownames(x, subsetName) <- value

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetRownames(x, subsetName)

# S4 method for class 'SubsetRangedSummarizedExperiment,character'
subsetRownames(x, subsetName) <- value

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetRownames(x, subsetName)

# S4 method for class 'SubsetSingleCellExperiment,character'
subsetRownames(x, subsetName) <- value

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetRownames(x, subsetName)

# S4 method for class 'SubsetSummarizedExperiment,character'
subsetRownames(x, subsetName) <- value

# S4 method for class 'SubsetSpatialExperiment,character'
subsetRownames(x, subsetName)

# S4 method for class 'SubsetSpatialExperiment,character'
subsetRownames(x, subsetName) <- value

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetRownames(x, subsetName)

# S4 method for class 'SubsetTreeSummarizedExperiment,character'
subsetRownames(x, subsetName) <- value
```

## Arguments

- x:

  Input `ExperimentSubset` object.

- subsetName:

  Name of the subset to get `rownames` from or set them to.

- value:

  Specify the rownames values to replace.

## Value

A `vector` of `rownames`.

Input object with rownames set to a subset.

## Examples

``` r
data(sce_chcl, package = "scds")
es <- ExperimentSubset(sce_chcl)
es <- createSubset(es,
"subset1",
rows = c(10,11,50,56,98,99,102,105,109, 200),
cols = c(20,21,40,45,90,99,100,123,166,299),
parentAssay = "counts")
subsetRownames(es, "subset1")
#>  [1] "OPCML"         "RP11-338K17.8" "RP11-770J1.3"  "NPSR1"        
#>  [5] "AC060835.1"    "GABRA2"        "PRDM16"        "CAB39"        
#>  [9] "LYRM1"         "RP1-178F10.3" 
```
