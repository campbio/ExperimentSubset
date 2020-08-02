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
#subset by continuous indices
es <- subsetAssay(es, "subset1", subsetRows = c(1:5), subsetCols = c(1:3))

#subset by intermittent indices
es <- subsetAssay(es, "subset2", subsetRows = c(2,3,4,5), subsetCols = c(4,5,6))

#subset by feature and cell names
es <- subsetAssay(es, "subset3", subsetRows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), subsetCols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"))

#subset by both indices and names
es <- subsetAssay(es, "subset3", subsetRows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), subsetCols = c(1:10))
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
