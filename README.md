# ExperimentSubset
This branch contains code for possible support of ExperimentSubset with SpatialExperiment & VisiumExperiment classes. It should be noted that this branch is not a replacement for the devel branch and contains quite possibly many errors and documentation issues.

Warning: SpatialExperiment class is not exported by default. Required changes can be seen in the repo below:
https://github.com/irzamsarfraz/SpatialExperiment

## Setting up SpatialExperiment & VisiumExperiment objects
```
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
```

Output:
```

```

## Creating ExperimentSubset object from VisiumExperiment
```
es <- ExperimentSubset(ve)
```

Output:
```
```
