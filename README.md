# ExperimentSubset
This branch contains code for possible support of ExperimentSubset with SpatialExperiment & VisiumExperiment classes. It should be noted that this branch is not a replacement for the devel branch and contains quite possibly many errors and documentation issues.

Warning: SpatialExperiment class is not exported by default. Required changes can be seen in the repo below:
https://github.com/irzamsarfraz/SpatialExperiment

## Setting up SpatialExperiment & VisiumExperiment objects
```
#loading libraries
library("SpatialExperiment")
library("Matrix")
library("rjson")


#Setting up SpatialExperiment object as described in the official vignette

fishCoordFile <- system.file(file.path("extdata", "seqFISH",
                            "fcortex.coordinates.txt"), 
                            package="SpatialExperiment")
fishCoordinates <- read.table(fishCoordFile, header=FALSE, sep=" ")
colnames(fishCoordinates) <- c("Cell_ID", "Irrelevant", "x", "y")

fishCellLabsFile <- system.file(file.path("extdata", "seqFISH", 
                            "seqfish_cell_labels.tsv"),
                            package="SpatialExperiment")
fishCellLabels <- read.table(file=fishCellLabsFile, header=FALSE, sep="\t")
colnames(fishCellLabels) <- c("Cell_ID", "cluster", "class", "classID", 
                                "Irrelevant", "Prob")
fishFeatCountsFile <- system.file(file.path("extdata", "seqFISH",
                            "seqfish_normalized_cortex_b2_testing.txt"), 
                            package="SpatialExperiment")
fishFeaturesCounts <- read.table(file=fishFeatCountsFile, 
                                header=FALSE, sep="\t", row.names=1)
                                
se <- SpatialExperiment(rowData=rownames(fishFeaturesCounts),
                        colData=fishCellLabels,
                        assays=SimpleList(counts=as.matrix(fishFeaturesCounts)),
                        spatialCoords=fishCoordinates)
                        
#Setting up VisiumExperiment object as described in the official vignette

barcodesFile <- system.file(file.path("extdata", "10x_visium",
                            "barcodes.tsv"), package="SpatialExperiment")
barcodesEx <- read.csv(barcodesFile, sep="\t", 
                     header=FALSE, col.names=c("Cell_ID"))
featuresFile <- system.file(file.path("extdata", "10x_visium",
                            "features.tsv"), package="SpatialExperiment")
featuresEx <- read.csv(featuresFile, sep="\t", 
                     header=FALSE, col.names=c("Feature_ID", "Feature_name", 
                                               "Feature_type"))

countsFile <- system.file(file.path("extdata", "10x_visium",
                            "matrix.mtx"), package="SpatialExperiment")
countsEx <- readMM(file=countsFile)
posFile <- system.file(file.path("extdata", "10x_visium",
                        "tissue_positions_list.tsv"), 
                        package="SpatialExperiment")
tissPosEx <- read.csv(posFile, 
                        sep="\t", header=FALSE, 
                        col.names=c("Cell_ID", "in_tissue", 
                                    "array_row", "array_col",
                                    "pxl_col_in_fullres", "pxl_row_in_fullres"))
imageFilePath <- list.files(system.file(file.path("extdata", "10x_visium", 
                              "images"), package="SpatialExperiment"), 
                            full.names=TRUE)
                            
scaleFile <- system.file(file.path("extdata", "10x_visium",
                        "scalefactors_json.json"), 
                        package="SpatialExperiment")

scalefactors <- fromJSON(file=scaleFile)

ve <- VisiumExperiment(rowData=featuresEx, colData=barcodesEx, 
                            assays=c(counts=countsEx), 
                            spatialCoords=tissPosEx,
                            scaleFactors=scalefactors)
```

## Creating ExperimentSubset object from SpatialExperiment
```
es <- ExperimentSubset(se)
es
```

Output:
```
class: SubsetSpatialExperiment 
dim: 113 1597 
metadata(0):
assays(1): counts
rownames(113): abca15 abca9 ... zfp715 zfp90
rowData names(1): X
colnames(1597): V2 V3 ... V1597 V1598
colData names(6): Cell_ID cluster ... Irrelevant Prob
reducedDimNames(0):
altExpNames(0):
spatialCoords(4): Cell_ID Irrelevant x y
subsets(0): 
subsetAssays(0): 
```

## Creating ExperimentSubset object from VisiumExperiment
```
es <- ExperimentSubset(ve)
es
```

Output:
```
class: SubsetVisiumExperiment 
dim: 50 50 
metadata(0):
assays(1): counts
rownames: NULL
rowData names(3): Feature_ID Feature_name Feature_type
colnames: NULL
colData names(1): Cell_ID
reducedDimNames(0):
altExpNames(0):
spatialCoords(6): Cell_ID in_tissue ... pxl_col_in_fullres pxl_row_in_fullres
inTissue(1): 22
imagePaths(0):
subsets(0): 
subsetAssays(0):
```

## Creating a subset from es object (from ve)
```
es <- createSubset(es, "subset1", cols = c(1:30))
es
```

Output:
```
class: SubsetVisiumExperiment 
dim: 50 50 
metadata(0):
assays(1): counts
rownames: NULL
rowData names(3): Feature_ID Feature_name Feature_type
colnames: NULL
colData names(1): Cell_ID
reducedDimNames(0):
altExpNames(0):
spatialCoords(6): Cell_ID in_tissue ... pxl_col_in_fullres pxl_row_in_fullres
inTissue(1): 22
imagePaths(0):
subsets(1): subset1
subsetAssays(1): subset1
```

## Showing internal structure of 'subset1' for demonstration purposes
```
es@subsets[["subset1"]]
```

Output:
```
An object of class "AssaySubset"
Slot "subsetName":
[1] "subset1"

Slot "rowIndices":
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42
[43] 43 44 45 46 47 48 49 50

Slot "colIndices":
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30

Slot "parentAssay":
[1] "counts"

Slot "internalAssay":
class: VisiumExperiment 
dim: 50 30 
metadata(0):
assays(0):
rownames: NULL
rowData names(0):
colnames: NULL
colData names(0):
reducedDimNames(0):
altExpNames(0):
spatialCoords(6): Cell_ID in_tissue ... pxl_col_in_fullres pxl_row_in_fullres
inTissue(1): 11
imagePaths(0):
```
