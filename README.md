# ExperimentSubset
Manages subsets of data with Bioconductor Experiment objects

# Sample Workflow
Load sample data:
```r
data(sce_chcl, package = "scds")
```
Load ExperimentSubset and create ES object:
```r
library(ExperimentSubset)
es <- ExperimentSubset::ExperimentSubset(assays = list(counts = assay(sce_chcl, "counts")))
```
Create multiple new subset assays:
```r
es <- subsetAssay(es, "subset1", rowIndices = c(1:5), colIndices = c(1:3))
es <- subsetAssay(es, "subset2", rowIndices = c(5:10), colIndices = c(10:12))
es <- subsetAssay(es, "subset3", rowIndices = c(20:23), colIndices = c(14:19))
```
See available subset assays:
```r
subsetNames(es)
```
Extract "subset1" subsetted assay:
```r
subsetAssay <- assay(es, "subset1")
```
# Features

# Installation
