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
es
```
<blockquote>
class: ExperimentSubset </br>
dim: 2000 2000 </br>
metadata(0):</br>
assays(1): counts</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(0):</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(0):</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):</br>
subsets(0): 
</blockquote>

Create multiple new subset assays:
```r
#subset by continuous indices
es <- subsetAssay(es, "subset1", subsetRows = c(1:5), subsetCols = c(1:3))

#subset by intermittent indices
es <- subsetAssay(es, "subset2", subsetRows = c(2,3,4,5), subsetCols = c(4,5,6))

#subset by feature and cell names
es <- subsetAssay(es, "subset3", subsetRows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), subsetCols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"))

#subset by both indices and names
es <- subsetAssay(es, "subset4", subsetRows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), subsetCols = c(1:10))

#view es object
es
```
<blockquote>
class: ExperimentSubset</br>
dim: 2000 2000 </br>
metadata(0):</br>
assays(1): counts</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(0):</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(0):</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):</br>
subsets(4): subset1 subset2 subset3 subset4
</blockquote>

See available subset assays:
```r
subsetNames(es)
```
>[1] "subset1" "subset2" "subset3" "subset4"

Extract subsetted assays:
```r
subsetAssay1 <- assay(es, "subset1")
subsetAssay2 <- assay(es, "subset2")
subsetAssay3 <- assay(es, "subset3")
subsetAssay4 <- assay(es, "subset4")
```

Store an external subset assay with different data:
```r
#imitate an external assay
logCounts <- log1p(assay(sce_chcl, "counts"))

#subset the assay
logCounts <- logCounts[1:10, 1:2]

#store subset through saveSubset function
es <- saveSubset(es, "logSubset1", logCounts)

#or store directly through assay<- function
assay(es, "logSubset2") <- logCounts

#view es object
es
```
<blockquote>
class: ExperimentSubset</br>
dim: 2000 2000</br>
metadata(0):</br>
assays(3): counts logSubset1_internal logSubset2_internal</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(0):</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(0):</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):</br>
subsets(6): subset1 subset2 subset3 subset4 logSubset1 logSubset2
</blockquote>

View internal storage mechanism of subsets:
```r
es@subsets$logSubset2
```
<blockquote>
An object of class "SingleCellSubset"</br></br>
Slot "subsetName":</br>
[1] "logSubset2"</br></br>
Slot "rowIndices":</br>
 [1]  1  2  3  4  5  6  7  8  9 10</br></br>
Slot "colIndices":</br>
[1] 1 2</br></br>
Slot "useAssay":</br>
[1] "logSubset2_internal"
</blockquote>

# Features

# Installation
