# ExperimentSubset
Manages subsets of data with Bioconductor Experiment objects

## Sample Workflow
**Load sample data:**
```r
data(sce_chcl, package = "scds")
```

**Load ExperimentSubset and create ES object:**
```r
library(ExperimentSubset)
es <- ExperimentSubset::ExperimentSubset(assays = list(counts = assay(sce_chcl, "counts"), logcounts = assay(sce_chcl, "logcounts")), colData=colData(sce_chcl), rowData= rowData(sce_chcl))
es
```
<blockquote>
class: ExperimentSubset</br>
dim: 2000 2000</br> 
metadata(0):</br>
assays(2): counts logcounts</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(1): gene</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(11): nGene nUMI ... ident doublet_true_labels</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):</br>
subsets(0):
</blockquote>

**Create multiple new subset assays:**
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
dim: 2000 2000</br> 
metadata(0):</br>
assays(2): counts logcounts</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(1): gene</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(11): nGene nUMI ... ident doublet_true_labels</br>
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

**Extract subsetted assays:**
```r
subsetAssay1 <- assay(es, "subset1")
subsetAssay2 <- assay(es, "subset2")
subsetAssay3 <- assay(es, "subset3")
subsetAssay4 <- assay(es, "subset4")
```

**Store an external subset assay with different data:**
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
assays(4): counts logcounts logSubset1_internal logSubset2_internal</br>
rownames(2000): C12orf73 RNU6-1256P ... FOXF1 PRR3</br>
rowData names(1): gene</br>
colnames(2000): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC ... GATGCTATCAGCAACT ATTTCTGGTGATGATA</br>
colData names(11): nGene nUMI ... ident doublet_true_labels</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):</br>
subsets(6): subset1 subset2 subset3 subset4 logSubset1 logSubset2
</blockquote>

**View internal storage mechanism of subsets:**
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

**View colData from a subset:**
```r
subsetColData(es, "logSubset1")
```
<blockquote>
 DataFrame with 2 rows and 11 columns</br>
                     nGene      nUMI    orig.ident hash_maxID hash_secondID      hash_margin hto_classification</br>
                 <integer> <numeric>      <factor>   <factor>      <factor>        <numeric>           <factor></br>
CTGCTGTCAGGGTATG      1583      3193 SeuratProject     THP1_A         KG1_A 2.20423805202945             THP1_A</br>
CAGTCCTTCGGTTAAC       775      1207 SeuratProject      KG1_A        THP1_A 1.75647763012464              KG1_A</br>
                 hto_classification_global  hash_ID    ident doublet_true_labels</br>
                                  <factor> <factor> <factor>         <character></br>
CTGCTGTCAGGGTATG                   Singlet   THP1_A  Singlet             Singlet</br>
CAGTCCTTCGGTTAAC                   Singlet    KG1_A  Singlet             Singlet
</blockquote>
                                   
**View rowData from a subset:**
```r
subsetRowData(es, "logSubset1")
```
<blockquote>
 DataFrame with 10 rows and 1 column</br>
                       gene</br>
                <character></br>
C12orf73           C12orf73</br>
RNU6-1256P       RNU6-1256P</br>
RN7SL749P         RN7SL749P</br>
RNU6-157P         RNU6-157P</br>
SP1                     SP1</br>
RP11-100E13.1 RP11-100E13.1</br>
AC005498.3       AC005498.3</br>
DCLRE1C             DCLRE1C</br>
RP11-1H15.2     RP11-1H15.2</br>
OPCML                 OPCML
</blockquote>

## Features

## Installation
