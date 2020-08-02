# ExperimentSubset
Efficiently subset single cell RNA seqeuncing data

# Sample Workflow
```r
library(ExperimentSubset)
scs <- ExperimentSubset::ExperimentSubset(assays = list(counts = counts))
scs <- subsetAssay(scs, "redDim", rowIndices = c(1:5), colIndices = c(1:3))
subsetNames(scs)
subsetAssay <- assay(scs, "redDim")
```

# Features

# Installation
