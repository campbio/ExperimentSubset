# Package index

## Constructor and classes

- [`ExperimentSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset.md)
  : ExperimentSubset constructor

- [`SubsetRangedSummarizedExperiment-class`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetRangedSummarizedExperiment-class.md)
  [`.SubsetRangedSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetRangedSummarizedExperiment-class.md)
  :

  An S4 class for `RangedSummarizedExperiment` objects with added
  support for subsets.

- [`SubsetSingleCellExperiment-class`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSingleCellExperiment-class.md)
  [`.SubsetSingleCellExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSingleCellExperiment-class.md)
  :

  An S4 class for `SingleCellExperiment` objects with added support for
  subsets.

- [`SubsetSpatialExperiment-class`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSpatialExperiment-class.md)
  [`.SubsetSpatialExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSpatialExperiment-class.md)
  :

  An S4 class for `SpatialExperiment` objects with added support for
  subsets.

- [`SubsetSummarizedExperiment-class`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSummarizedExperiment-class.md)
  [`.SubsetSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetSummarizedExperiment-class.md)
  :

  An S4 class for `SummarizedExperiment` objects with added support for
  subsets.

- [`SubsetTreeSummarizedExperiment-class`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetTreeSummarizedExperiment-class.md)
  [`.SubsetTreeSummarizedExperiment`](https://campbio.github.io/ExperimentSubset/dev/reference/SubsetTreeSummarizedExperiment-class.md)
  :

  An S4 class for `TreeSummarizedExperiment` objects with added support
  for subsets.

## Creating and managing subsets

- [`createSubset()`](https://campbio.github.io/ExperimentSubset/dev/reference/createSubset.md)
  : Subset creation method for ExperimentSubset objects
- [`subsetNames()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetNames.md)
  : Get names of only the subsets in ExperimentSubset objects
- [`subsetCount()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetCount.md)
  : Subset count method for ExperimentSubset objects
- [`subsetDim()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetDim.md)
  : Get dimensions of subsets in ExperimentSubset objects
- [`subsetParent()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetParent.md)
  : Subset parent hierarchy retrieval method for ExperimentSubset
  objects
- [`subsetSummary()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetSummary.md)
  : Method for displaying 'child-parent' link structure of subsets in
  ExperimentSubset objects
- [`subsetAssayNames()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetAssayNames.md)
  : Name retrieval method for all subset assays in ExperimentSubset
  objects
- [`subsetAssayCount()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetAssayCount.md)
  : Count method for subset assays in ExperimentSubset objects
- [`` `[`( ``*`<SubsetRangedSummarizedExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-subsetting.md)
  [`` `[`( ``*`<SubsetSingleCellExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-subsetting.md)
  [`` `[`( ``*`<SubsetSummarizedExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-subsetting.md)
  [`` `[`( ``*`<SubsetSpatialExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-subsetting.md)
  [`` `[`( ``*`<SubsetTreeSummarizedExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-subsetting.md)
  : Subset an ExperimentSubset object by rows or columns
- [`cbind(`*`<SubsetRangedSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`rbind(`*`<SubsetRangedSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`cbind(`*`<SubsetSingleCellExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`rbind(`*`<SubsetSingleCellExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`cbind(`*`<SubsetSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`rbind(`*`<SubsetSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`cbind(`*`<SubsetSpatialExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`rbind(`*`<SubsetSpatialExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`cbind(`*`<SubsetTreeSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  [`rbind(`*`<SubsetTreeSummarizedExperiment>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/ExperimentSubset-combine.md)
  : Combine ExperimentSubset objects by columns or rows

## Getting and setting subset data

- [`getSubsetAssay()`](https://campbio.github.io/ExperimentSubset/dev/reference/getSubsetAssay.md)
  :

  Get subset assay from an `ExperimentSubset` object.

- [`setSubsetAssay()`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetAssay.md)
  : Method for storing new assays inside subsets in ExperimentSubset
  objects

- [`getSubsetExperiment()`](https://campbio.github.io/ExperimentSubset/dev/reference/getSubsetExperiment.md)
  : Take a subset out as a standalone experiment object

- [`setSubsetExperiment()`](https://campbio.github.io/ExperimentSubset/dev/reference/setSubsetExperiment.md)
  : Store a standalone experiment object back into a subset

- [`subsetRowData()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRowData.md)
  [`` `subsetRowData<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRowData.md)
  : Accessor method for rowData from subsets in ExperimentSubset objects

- [`subsetColData()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColData.md)
  [`` `subsetColData<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColData.md)
  : Accessor method for colData from subsets in ExperimentSubset objects

- [`subsetRownames()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRownames.md)
  [`` `subsetRownames<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRownames.md)
  : subsetRownames

- [`subsetColnames()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColnames.md)
  [`` `subsetColnames<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColnames.md)
  : subsetColnames

## Class-specific accessors

Tree links of a TreeSummarizedExperiment and spatial data of a
SpatialExperiment, for a subset.

- [`subsetRowLinks()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetRowLinks.md)
  : Accessor method for rowLinks from subsets in ExperimentSubset
  objects
- [`subsetColLinks()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetColLinks.md)
  : Accessor method for colLinks from subsets in ExperimentSubset
  objects
- [`subsetSpatialCoords()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetSpatialCoords.md)
  : Accessor method for spatialCoords from subsets in ExperimentSubset
  objects
- [`subsetSpatialData()`](https://campbio.github.io/ExperimentSubset/dev/reference/subsetSpatialData.md)
  [`` `subsetSpatialData<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/subsetSpatialData.md)
  : Deprecated: spatialData of subsets in ExperimentSubset objects

## Methods that take subsetName

Methods of Bioconductor generics that read or set a subset’s data when
given subsetName, and otherwise act on the main object. Besides those
below, these include assay(), assay\<-, rowData(), colData(),
reducedDim(), reducedDims(), altExp(), altExps(), altExpNames(),
metadata() and their setters.

- [`` `reducedDimNames<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/reducedDimNames-set.md)
  : reducedDimNames\<-
- [`reducedDimNames()`](https://campbio.github.io/ExperimentSubset/dev/reference/reducedDimNames.md)
  : reducedDimNames
- [`` `reducedDims<-`( ``*`<ANY>`*`)`](https://campbio.github.io/ExperimentSubset/dev/reference/reducedDims-set-ANY-method.md)
  : reducedDims\<-
- [`` `reducedDims<-`() ``](https://campbio.github.io/ExperimentSubset/dev/reference/reducedDims-set.md)
  : reducedDims\<-
