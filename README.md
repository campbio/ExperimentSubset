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
es <- createSubset(es, "subset1", rows = c(1:5), cols = c(1:3), useAssay = "counts")

#subset by intermittent indices
es <- createSubset(es, "subset2", rows = c(2,3,4,5), cols = c(4,5,6), useAssay = "counts")

#subset by feature and cell names
es <- createSubset(es, "subset3", rows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), useAssay = "counts")

#subset by both indices and names
es <- createSubset(es, "subset4", rows = c("C12orf73", "RNU6-1256P", "RN7SL749P", "RNU6-157P"), cols = c(1:10), useAssay = "counts")

#or use assay<- function to input a subset assay using "useAssay" parameter as long as the assay already resides inside the es object
c <- assay(es, "counts")
c <- c[1:10, 1:5]
assay(es, "subset5", useAssay = "counts") <- c

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
subsets(4): subset1 subset2 subset3 subset4 subset5
</blockquote>


**See available subset assays:**
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
subsetAssay5 <- assay(es, "subset5")
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
subsets(6): subset1 subset2 subset3 subset4 subset5 logSubset1 logSubset2
</blockquote>


**View internal storage mechanism of subsets:**
```r
es@subsets$subset1
```
<blockquote>
 An object of class "SingleCellSubset"</br>
Slot "subsetName":</br>
[1] "subset1"</br></br>
Slot "rowIndices":</br>
[1] 1 2 3 4 5</br></br>
Slot "colIndices":</br>
[1] 1 2 3</br></br>
Slot "useAssay":</br>
[1] "counts"</br></br>
Slot "internalAssay":</br>
class: SingleCellExperiment</br> 
dim: 5 3 </br>
metadata(0):</br>
assays(1): counts</br>
rownames: NULL</br>
rowData names(0):</br>
colnames: NULL</br>
colData names(0):</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):
</blockquote>

```r
es@subsets$logSubset2
```
<blockquote>
An object of class "SingleCellSubset"</br>
Slot "subsetName":</br>
[1] "logSubset2"</br></br>
Slot "rowIndices":</br>
 [1]  1  2  3  4  5  6  7  8  9 10</br></br>
Slot "colIndices":</br>
[1] 1 2</br></br>
Slot "useAssay":</br>
NULL</br></br>
Slot "internalAssay":</br>
class: SingleCellExperiment</br> 
dim: 10 2 </br>
metadata(0):</br>
assays(1): counts</br>
rownames(10): C12orf73 RNU6-1256P ... RP11-1H15.2 OPCML</br>
rowData names(0):</br>
colnames(2): CTGCTGTCAGGGTATG CAGTCCTTCGGTTAAC</br>
colData names(0):</br>
reducedDimNames(0):</br>
spikeNames(0):</br>
altExpNames(0):
</blockquote>


**View colData from a subset:**
```r
#either use subsetcolData
subsetColData(es, "logSubset1")

#or use colData with subsetName parameter
colData(es, subsetName = "logSubset1")
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
#either use subsetRowData
subsetRowData(es, "logSubset1")

#or use rowData with subsetName parameter
rowData(es, subsetName = "logSubset1")
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

**Add new colData or rowData only to the existing subset**
```r
#add a "DummyGeneID" column to rowData only against the "subset1"
rowData(es, subsetName = "subset1") <- DataFrame(DummyGeneID = c("a", "b", "c", "d", "e"))

#add a "SampleID" column to colData only against the "subset1"
colData(es, subsetName = "subset1") <- DataFrame(SampleID = c("1", "2", "3"))
```

**View subset colData or rowData**
```r
#displays the previously added "DummyGeneID" column as well which is only present against the "subset1"
rowData(es, subsetName = "subset1")

#similarly for colData
colData(es, subsetName = "subset1")
```
<blockquote>
 DataFrame with 5 rows and 2 columns</br>
                  gene DummyGeneID</br>
           <character> <character></br>
C12orf73      C12orf73           a</br>
RNU6-1256P  RNU6-1256P           b</br>
RN7SL749P    RN7SL749P           c</br>
RNU6-157P    RNU6-157P           d</br>
SP1                SP1           e
</blockquote>

**Sample Nested Assay Script for Testing**
```r
es <- ExperimentSubset::ExperimentSubset(assays = list(counts = assay(sce_chcl, "counts"), logcounts = assay(sce_chcl, "logcounts")), colData=colData(sce_chcl), rowData= rowData(sce_chcl))
es

#subset1
es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "counts")

#scaledSubset1
scaledCounts <- assay(es, "subset1")
scaledCounts[,] <- scaledCounts[,] + 1
es <- storeSubset(es, "subset1", scaledCounts, "scaledSubset1")


#topHVGSubset1
es <- createSubset(es, "topHVGSubset1", rows = c("NPSR1", "GABRA2", "CAB39", "LYRM1", "RP1-178F10.3"), cols = c(1:10), parentAssay = "subset1")

#topHVGSubset2
scaledCounts <- assay(es, "topHVGSubset1")
scaledCounts[,] <- scaledCounts[,] + 1.77
es <- storeSubset(es, "topHVGSubset1", scaledCounts, "topHVGSubset2")

#subsetHVG
es <- createSubset(es, "subsetHVG", rows = c(1:5), cols = c(1:2), parentAssay = "topHVGSubset1")


#scaledTopHVGSubset1
scaledCounts <- assay(es, "subsetHVG")
scaledCounts[,] <- scaledCounts[,] + 56
es <- storeSubset(es, "subsetHVG", scaledCounts, "scaledTopHVGSubset1")


#smallScaledSubset1
es <- createSubset(es, "smallScaledSubset1", c("NPSR1", "LYRM1"), c(3:4), "scaledSubset1") 


#logSmallScaledSubset1
scaledCounts <- assay(es, "smallScaledSubset1")
scaledCounts[,] <- scaledCounts[,] + 1000
es <- storeSubset(es, "smallScaledSubset1", scaledCounts, "logSmallScaledSubset1")

```

## Features

## Installation
