library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
library(scds)
library(TENxPBMCData)
library(scater)
library(airway)

#load data once
data(sce_chcl, package = "scds")
tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")

context("Testing ExperimentSubset functions")

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SingleCellExperiment object",{
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "SubsetSingleCellExperiment")
})

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SummarizedExperiment object",{
  se <- SummarizedExperiment::SummarizedExperiment(sce_chcl)
  testthat::expect_error(es <- ExperimentSubset(se), "Experiment objects with unnamed assays are not supported. Use assayNames<- setter method to set assay names before creating ES object.")
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = assay(sce_chcl, "counts"), logcounts = assay(sce_chcl, "logcounts")))
  testthat::expect_s4_class(es <- ExperimentSubset(se), "SubsetSummarizedExperiment")
  testthat::expect_true(validObject(es))
  testthat::expect_equal(length(metadata(es)), 0)
  testthat::expect_s4_class(assay(es, "counts"), "dgCMatrix")
  es <- createSubset(es, "s1", rows = c(1:10))
  testthat::expect_true("s1" %in% subsetNames(es))
  testthat::expect_true("s1" %in% subsetAssayNames(es))
  show(es)
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "newS1")
  testthat::expect_true("newS1" %in% subsetAssayNames(es))
  testthat::expect_false("newS1" %in% subsetNames(es))
  assay(es, "assay2") <- assay(es, "counts")
  testthat::expect_true("assay2" %in% assayNames(es))
  testthat::expect_equal(length(subsetParent(es, "s1")), 1)
  subsetSummary(es)
  testthat::expect_equal(subsetDim(es, "s1")[1], 10)
  testthat::expect_equal(subsetDim(es, "s1")[2], 2000)
  testthat::expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  testthat::expect_equal(length(metadata(es)), 1)
  testthat::expect_s4_class(getSubsetAssay(es, "s1"), "dgCMatrix")
  subsetRownames(es, "s1") <- subsetRownames(es, "s1")
  testthat::expect_equal(length(subsetRownames(es, "s1")), 10)
  subsetColnames(es, "s1") <- subsetColnames(es, "s1")
  testthat::expect_equal(length(subsetColnames(es, "s1")), 2000)
  subsetRownames(es, "newS1") <- subsetRownames(es, "newS1")
  testthat::expect_equal(length(subsetRownames(es, "newS1")), 10)
  subsetColnames(es, "newS1") <- subsetColnames(es, "newS1")
  testthat::expect_equal(length(subsetColnames(es, "newS1")), 2000)
  colData(es, subsetName = "s1") <- cbind(colData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[2]))))
  testthat::expect_equal(ncol(colData(es, subsetName = "s1")), 1)
  rowData(es, subsetName = "s1") <- cbind(rowData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[1]))))
  testthat::expect_equal(ncol(rowData(es, subsetName = "s1")), 1)
  
  
  data(airway, package="airway")
  rse <- airway
  testthat::expect_s4_class(es <- ExperimentSubset(rse), "SubsetRangedSummarizedExperiment")
  subsetSummary(es)
  testthat::expect_equal(length(metadata(es)), 1)

  es <- createSubset(es, "s2", rows = c(1:15))
  testthat::expect_true("s2" %in% subsetNames(es))
  testthat::expect_true("s2" %in% subsetAssayNames(es))
  show(es)
  es <- setSubsetAssay(es, "s2", assay(es, "s2"), "newS2")
  testthat::expect_true("newS2" %in% subsetAssayNames(es))
  testthat::expect_false("newS2" %in% subsetNames(es))
  assay(es, "assay2") <- assay(es, "counts")
  testthat::expect_true("assay2" %in% assayNames(es))
  testthat::expect_equal(length(subsetParent(es, "s2")), 1)
  subsetSummary(es)
  testthat::expect_equal(subsetDim(es, "s2")[1], 15)
  testthat::expect_equal(subsetDim(es, "s2")[2], 8)
  testthat::expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  testthat::expect_equal(length(metadata(es)), 1)

  subsetRownames(es, "s2") <- subsetRownames(es, "s2")
  testthat::expect_equal(length(subsetRownames(es, "s2")), 15)
  subsetColnames(es, "s2") <- subsetColnames(es, "s2")
  testthat::expect_equal(length(subsetColnames(es, "s2")), 8)
  subsetRownames(es, "newS2") <- subsetRownames(es, "newS2")
  testthat::expect_equal(length(subsetRownames(es, "newS2")), 15)
  subsetColnames(es, "newS2") <- subsetColnames(es, "newS2")
  testthat::expect_equal(length(subsetColnames(es, "newS2")), 8)
  colData(es, subsetName = "s2") <- cbind(colData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[2]))))
  testthat::expect_equal(ncol(colData(es, subsetName = "s2", parentColData = TRUE)), 10)
  rowData(es, subsetName = "s2") <- cbind(rowData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[1]))))
  testthat::expect_equal(ncol(rowData(es, subsetName = "s2")), 1)
  
  testthat::expect_error(es <- createSubset(es, "s2", rows = c(1:15)),
                         "A subset with the specified 'subsetName' parameter already exists
         in the object. Subset names must be unique.")
  testthat::expect_error(es <- createSubset(es, subsetName = 12, rows = c(1:15)),
                         "'subsetName' parameter must be a unique a character value.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", rows = list()),
                         "'rows' parameter must be either a 'NULL' to include all rows,
         or a numeric vector or a character vector that specify the rows to 
         include in the subset.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", cols = list()),
                         "'cols' parameter must be either a 'NULL' to include all rows,
         or a numeric vector or a character vector that specify the columns to 
         include in the subset.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", parentAssay = 1),
                         "'parentAssay' parameter can either be 'NULL' to use the default
         assay in the input object or a character value that specifies the
         parentAssay to use from parent object.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", rows = c(64103:64103)),
                         "Selected rows not available in the specified assay.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", rows = c(seq(1:64102), 64102)),
                         "Selected rows not available in the specified assay.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", rows = c(1,3, NA)),
                         "Selected rows not available in the specified assay.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", cols = c(10:12)),
                         "Selected columns not available in the specified assay.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", cols = c(seq(1:8), 8)),
                         "Selected columns not available in the specified assay.")
  testthat::expect_error(es <- createSubset(es, subsetName = "s3", cols = c(1,3, NA)),
                         "Selected columns not available in the specified assay.")
  
})

testthat::test_that("Testing ExperimentSubset constructor with subset & createSubset without parentAssay",{
  es <- ExperimentSubset(sce_chcl, subset = list(subsetName = "subset10", rows = c(1:10), cols = c(1:2), parentAssay = "counts"))
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset10")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(1:10))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(1:2))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109, 200))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
})

testthat::test_that("Testing createSubset and setSubsetAssay",{
  es <- ExperimentSubset(sce_chcl)
  #createSubset
  es <- createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109,200))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  subsetRownames(es, "subset1") <- subsetRownames(es, "subset1")
  subsetColnames(es, "subset1") <- subsetColnames(es, "subset1")
  subsetRownames(es, "subset1") <- subsetRownames(es, "subset1")
  subsetColnames(es, "subset1") <- subsetColnames(es, "subset1")
  
  es <- createSubset(es, "subset2",
                                       rows = c("OPCML", "RP11-338K17.8", "RP11-770J1.3", "NPSR1", "AC060835.1", "GABRA2", "PRDM16", "CAB39", "LYRM1", "RP1-178F10.3"),
                                       cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
                                       parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset2@subsetName, "subset2")
  testthat::expect_equal(es@subsets$subset2@rowIndices, c(10,11,50,56,98,99,102,105,109,200))
  testthat::expect_equal(es@subsets$subset2@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset2@parentAssay, "counts")
  
  # expect_error(es <- createSubset(es, "subset3",
  #                                                   rows = c("OPCML", "RP11-338K17.8", 50,56,98,99,102,105,109,200),
  #                                                   cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
  #                                                   parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  
  #setSubsetAssay
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  expect_error(es <- setSubsetAssay(es, "subset4", counts1p, "scaledSubset1"),
               "subset4 does not exist in the subsets slot of the object.")
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  expect_true("scaledSubset1" %in% assayNames(es@subsets$subset1@internalAssay))
  getSubsetAssay(es, "scaledSubset1")
  # expect_error(es <- setSubsetAssay(es, "subset1", counts1p, subsetAssayName = NULL),
  #              "subset1 already exists. Please choose a different subsetName parameter.")
  # es <- setSubsetAssay(es, "subset4", counts1p, subsetAssayName = NULL)
  # expect_equal(rownames(counts1p), rownames(es, "subset4"))
  # expect_equal(colnames(counts1p), colnames(es, "subset4"))
  #cant have null now
})

testthat::test_that("Testing rownames, colnames, rowData and colData",{
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "subset1",
                     rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  
  # expect_equal(rownames(es, subsetName = "scaledSubset1"), rownames(counts1p))
  # expect_equal(colnames(es, subsetName = "scaledSubset1"), colnames(counts1p))
  # expect_equal(rownames(es, subsetName = "subset1"), rownames(es)[c(10,11,50,56,98,99,102,105,109,200)])
  # expect_equal(colnames(es, subsetName = "subset1"), colnames(es)[c(20,21,40,45,90,99,100,123,166,299)])
  testthat::expect_error(colData(es, subsetName = "subset1Internal"),
                         "subset1Internal subset does not exist.")
  testthat::expect_error(rowData(es, subsetName = "subset1Internal"),
                         "subset1Internal subset does not exist.")
  rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[1])))
  colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[2])))
  expect_equal(ncol(rowData(es, subsetName = "subset1", parentRowData = TRUE)), 2)
  expect_equal(nrow(rowData(es, subsetName = "subset1", parentRowData = TRUE)), subsetDim(es, subsetName = "subset1")[1])
  expect_equal(ncol(colData(es, subsetName = "subset1", parentColData = TRUE)), 12)
  expect_equal(nrow(rowData(es, subsetName = "subset1", parentColData = TRUE)), subsetDim(es, subsetName = "subset1")[2])
})

testthat::test_that("Testing subset helper/supplementary functions",{
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  expect_equal(subsetCount(es), 1)
  expect_equal(subsetAssayCount(es), 2)
  expect_equal(subsetAssayNames(es), c("subset1", "scaledSubset1"))
  expect_equal(subsetNames(es), "subset1")
  subsetSummary(es)
  expect_equal(subsetParent(es, "scaledSubset1"), list("subset1", "counts"))
  show(es)
})

testthat::test_that("Testing assay function",{
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  es@subsets$subset1@parentAssay <- "counts"
  subset1Assay <- assay(es, "subset1")
  expect_true(all.equal(subset1Assay, assay(es, "counts")[c(10,11,50,56,98,99,102,105,109,200), c(20,21,40,45,90,99,100,123,166,299)]))
  assay(es, "counts2") <- assay(es, "counts")
  expect_true(all.equal(assay(es, "counts"), assay(es, "counts2")))
})

#create es with pbmc4k data, compute colsums and rowsums, create a hvg1000 subset for further testing
es <- ExperimentSubset(tenx_pbmc4k)
colData(es) <- cbind(colData(es), colSums = colSums(assay(es, "counts")))
rowData(es) <- cbind(rowData(es), rowSums = rowSums(assay(es, "counts")))
es <- createSubset(es, "filteredAssay", rows = rownames(es), cols = which(colData(es)$colSums > mean(colData(es)$colSums)), parentAssay = "counts")
#normalize and choose HVGs with scrapper, as in the vignette
filteredCounts <- assay(es, "filteredAssay")
filteredSizeFactors <- scrapper::centerSizeFactors(colSums(filteredCounts))
assay(
  es, "filteredAssay", withDimnames = FALSE,
  subsetAssayName = "filteredAssayNormalized"
) <- scrapper::normalizeCounts(
  filteredCounts, size.factors = filteredSizeFactors
)
geneVariances <- scrapper::modelGeneVariances(
  assay(es, "filteredAssayNormalized")
)
hvgIndices <- scrapper::chooseHighlyVariableGenes(
  geneVariances$statistics$residuals, top = 1000
)
#pass gene names to keep character rows with a subset parent covered
hvgNames <- rownames(assay(es, "filteredAssayNormalized"))[hvgIndices]
stopifnot(length(hvgNames) > 0, !anyNA(hvgNames))
es <- createSubset(
  es, "hvg1000", rows = hvgNames, cols = seq_len(1622),
  parentAssay = "filteredAssayNormalized"
)

testthat::test_that("Testing reducedDims",{
  set.seed(20)
  colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
  expect_true("cluster" %in% colnames(colData(es, subsetName = "hvg1000")))

  reducedDim(es, type = "PCA", subsetName = "hvg1000") <- scater::calculatePCA(assay(es, "hvg1000"))
  expect_true("PCA" %in% reducedDimNames(es, subsetName = "hvg1000"))
  reducedDims(es, subsetName = "hvg1000") <- list(a = scater::calculatePCA(assay(es, "hvg1000")), b = scater::calculatePCA(assay(es, "hvg1000")))
  expect_true(all.equal(reducedDimNames(es, subsetName = "hvg1000"), c("a","b")))
  
  reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1", "PCA_2")
  expect_true(all.equal(reducedDimNames(es, subsetName = "hvg1000"), c("PCA_1", "PCA_2")))

  testthat::expect_error(reducedDimNames(es, subsetName = "hvg123"),
                         "hvg123 subset does not exist.")
  
  testthat::expect_error(reducedDimNames(es, subsetName = "hvg123") <- c("PCA_1"),
                         "hvg123 does not exist in the subsets slot of the object.")
  
  reducedDim(es, type = "PCA") <- scater::calculatePCA(assay(es, "counts"))
  expect_true(!is.null(reducedDim(es, type = "PCA")))

  reducedDims(es) <- list(a = scater::calculatePCA(assay(es, "counts")), b = scater::calculatePCA(assay(es, "counts")))
  expect_true(!is.null(reducedDims(es)))
  expect_true(all.equal(reducedDimNames(es), c("a", "b")))

  reducedDimNames(es) <- c("PCA_Counts_1", "PCA_Counts_2")
  expect_true(all.equal(reducedDimNames(es), c("PCA_Counts_1", "PCA_Counts_2")))

})

testthat::test_that("Testing metadata",{
  metadata(es, subsetName = "hvg1000") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(metadata(es, subsetName = "hvg1000")$meta1, "This is testing meta1")
  expect_equal(metadata(es, subsetName = "hvg1000")$meta2, "This is testing meta2")
  
  testthat::expect_error(metadata(es, subsetName = "hvg") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2"),
                         "hvg subset does not exist.")
  
  testthat::expect_error(metadata(es, subsetName = "hvg"),
                         "hvg subset does not exist.")
  
  metadata(es) <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(metadata(es)$meta1, "This is testing meta1")
  expect_equal(metadata(es)$meta2, "This is testing meta2")
})

testthat::test_that("Testing altExperiments",{
  altExp(es, "alt1", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  altExp(es, "alt2", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  expect_true(all.equal(altExpNames(es, subsetName = "hvg1000"), c("alt1", "alt2")))

  testthat::expect_error(altExp(es, "alt1", subsetName = "hvg") <- es@subsets$hvg1000@internalAssay,
                         "hvg subset does not exist.")
  
  altExp(es, subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  expect_equal(length(altExpNames(es, subsetName = "hvg1000")), 2)
  
  expect_true(!is.null(altExp(es, subsetName = "hvg1000")))
  
  altExp(es) <- tenx_pbmc4k
  expect_true(!is.null(altExp(es)))

  altExp(es, "a1") <- tenx_pbmc4k
  expect_true("a1" %in% altExpNames(es))
  expect_true(!is.null(altExp(es, "a1")))
  
  testthat::expect_error(altExpNames(es, subsetName = "na"),
                         "na subset does not exist.")
  
  testthat::expect_error(altExpNames(es, subsetName = "na") <- c("na2"),
                         "na subset does not exist.")
  
  altExpNames(es) <- c("a1", "a2")
  expect_true(all.equal(altExpNames(es), c("a1", "a2")))
  
  testthat::expect_error(altExp(es, "alt1", subsetName = "hvg10002"),
                         "hvg10002 subset does not exist.")
  
  altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2")
  expect_true(all.equal(altExpNames(es, subsetName = "hvg1000"), c("a1", "a2")))

  testthat::expect_error(altExpNames(es, subsetName = "hvg1000x"),
                         "hvg1000x subset does not exist.")
  
  altExps(es, subsetName = "hvg1000") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay)
  expect_true(all.equal(altExpNames(es, subsetName = "hvg1000"), c("a","b")))

  testthat::expect_error(altExps(es, subsetName = "hvg") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay),
                         "hvg subset does not exist.")
  
  testthat::expect_error(altExps(es, subsetName = "hvg10002"),
                         "hvg10002 subset does not exist.")
  
  altExps(es) <- list(a = tenx_pbmc4k, b = tenx_pbmc4k)
  expect_true(all.equal(altExpNames(es), c("a", "b")))
} 
)

testthat::test_that("altExp gets and sets the altExp chosen by name or index", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "small", rows = 1:20, cols = 1:10)
  # altExps told apart by their number of rows
  makeAlt <- function(value, nrow, cells) {
    counts <- matrix(value, nrow, length(cells), dimnames = list(NULL, cells))
    SingleCellExperiment(list(counts = counts))
  }
  subsetCells <- subsetColnames(es, "small")
  alt1 <- makeAlt(1, 2, subsetCells)
  alt2 <- makeAlt(2, 3, subsetCells)
  alt3 <- makeAlt(3, 4, subsetCells)
  altExp(es, "alt1", subsetName = "small") <- alt1
  altExp(es, "alt2", subsetName = "small") <- alt2

  testthat::expect_equal(nrow(altExp(es, "alt2", subsetName = "small")), 3L)
  testthat::expect_equal(nrow(altExp(es, 2, subsetName = "small")), 3L)
  testthat::expect_equal(nrow(altExp(es, subsetName = "small")), 2L)

  altExp(es, 2, subsetName = "small") <- alt3
  testthat::expect_equal(nrow(altExp(es, "alt1", subsetName = "small")), 2L)
  testthat::expect_equal(nrow(altExp(es, "alt2", subsetName = "small")), 4L)

  # same calls on the main object, without subsetName
  altExp(es, "a1") <- makeAlt(1, 2, colnames(es))
  altExp(es, "a2") <- makeAlt(2, 3, colnames(es))
  altExp(es, 2) <- makeAlt(3, 4, colnames(es))
  testthat::expect_equal(nrow(altExp(es, "a1")), 2L)
  testthat::expect_equal(nrow(altExp(es, 2)), 4L)
})

testthat::test_that("Accessors on the main object pass extra arguments on", {
  # withDimnames is passed through the generic's ... for these getters
  sce <- sce_chcl
  reducedDim(sce, "PCA") <- matrix(seq_len(ncol(sce) * 2), ncol = 2)
  # stored altExp colnames differ, so withDimnames changes the result
  altCells <- paste0("alt_", seq_len(ncol(sce)))
  counts <- matrix(1, 2, ncol(sce), dimnames = list(NULL, altCells))
  alt <- SingleCellExperiment(list(counts = counts))
  altExps(sce, withDimnames = FALSE) <- list(alt1 = alt)
  es <- ExperimentSubset(sce)

  testthat::expect_equal(reducedDim(es, "PCA", withDimnames = FALSE),
                         reducedDim(sce, "PCA", withDimnames = FALSE))
  testthat::expect_equal(reducedDims(es, withDimnames = FALSE),
                         reducedDims(sce, withDimnames = FALSE))
  testthat::expect_equal(altExps(es, withDimnames = FALSE),
                         altExps(sce, withDimnames = FALSE))
  testthat::expect_equal(colnames(altExps(es, withDimnames = FALSE)[[1]]),
                         altCells)
})

testthat::test_that("Subsets of subsets use the right rows and columns", {
  sce <- sce_chcl
  rowData(sce)$geneIndex <- seq_len(nrow(sce))
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 101:200, cols = 501:600)
  # s2 is made from the subset s1
  es <- createSubset(es, "s2", rows = 1:10, cols = 1:10, parentAssay = "s1")
  # s3 is made from s2, three levels below the main assay
  es <- createSubset(es, "s3", rows = 2:3, cols = 4:5, parentAssay = "s2")
  # s4 is made from an assay stored inside s1
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "s1copy")
  es <- createSubset(es, "s4", rows = 5:6, cols = 2:3, parentAssay = "s1copy")

  testthat::expect_equal(rownames(assay(es, "s2")), rownames(sce)[101:110])
  testthat::expect_equal(subsetRownames(es, "s2"), rownames(sce)[101:110])
  testthat::expect_equal(subsetColnames(es, "s2"), colnames(sce)[501:510])
  testthat::expect_equal(subsetRownames(es, "s3"), rownames(sce)[102:103])
  testthat::expect_equal(subsetColnames(es, "s3"), colnames(sce)[504:505])
  testthat::expect_equal(subsetRownames(es, "s4"), rownames(sce)[105:106])
  testthat::expect_equal(subsetColnames(es, "s4"), colnames(sce)[502:503])

  s2RowData <- rowData(es, subsetName = "s2", parentRowData = TRUE)
  s2ColData <- colData(es, subsetName = "s2", parentColData = TRUE)
  testthat::expect_equal(s2RowData$geneIndex, 101:110)
  testthat::expect_equal(s2ColData$cellIndex, 501:510)
  s3ColData <- colData(es, subsetName = "s3", parentColData = TRUE)
  testthat::expect_equal(s3ColData$cellIndex, 504:505)
})

testthat::test_that("Subsets of subsets use the right spatial coordinates", {
  spots <- paste0("spot", 1:30)
  counts <- matrix(seq_len(600), nrow = 20,
                   dimnames = list(paste0("gene", 1:20), spots))
  coords <- cbind(x = as.numeric(1:30), y = as.numeric(1:30) * 10)
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 15),
                           spatialCoords = coords)
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", rows = 1:10, cols = 11:20)
  es <- createSubset(es, "s2", rows = 1:5, cols = 4:8, parentAssay = "s1")

  testthat::expect_equal(subsetColnames(es, "s2"), spots[14:18])
  testthat::expect_equal(unname(subsetSpatialCoords(es, "s2")[, "x"]),
                         as.numeric(14:18))
  testthat::expect_equal(colData(es, subsetName = "s2")$sample_id,
                         c("A", "A", "B", "B", "B"))
})

testthat::test_that("Subsets of subsets use the right rowLinks", {
  data("tinyTree", package = "TreeSummarizedExperiment")
  counts <- matrix(seq_len(100), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:10)))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  es <- createSubset(es, "s1", rows = 3:8, cols = 1:10)
  es <- createSubset(es, "s2", rows = 2:4, cols = 1:5, parentAssay = "s1")

  links <- subsetRowLinks(es, "s2", parentRowLinkData = TRUE)
  testthat::expect_equal(links$nodeLab, rowLinks(tse)$nodeLab[4:6])
})

testthat::test_that("Parent lookups follow .assay() for ambiguous names", {
  sce <- sce_chcl
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 101:200, cols = 501:600)
  es <- createSubset(es, "t1", rows = 1:10, cols = 1:10)
  # two subsets store an assay named "dup"; .assay() uses the last one
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "dup")
  es <- setSubsetAssay(es, "t1", assay(es, "t1"), "dup")
  es <- createSubset(es, "s5", rows = 2:3, cols = 2:3, parentAssay = "dup")
  testthat::expect_equal(subsetRownames(es, "s5"), rownames(assay(es, "s5")))
  testthat::expect_equal(subsetRownames(es, "s5"), rownames(sce)[2:3])

  # renaming the main assay must not break subsets made from it
  assayNames(es)[1] <- "raw"
  t1ColData <- colData(es, subsetName = "t1", parentColData = TRUE)
  testthat::expect_equal(t1ColData$cellIndex, 1:10)
})

testthat::test_that("Subsets with one row or column keep their dimensions", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "oneRow", rows = 5, cols = 1:10)
  es <- createSubset(es, "oneCol", rows = 1:10, cols = 5)
  testthat::expect_equal(dim(assay(es, "oneRow")), c(1L, 10L))
  testthat::expect_equal(dim(assay(es, "oneCol")), c(10L, 1L))
  # a subset made from a one-row subset
  es <- createSubset(es, "child", rows = 1, cols = 1:3, parentAssay = "oneRow")
  testthat::expect_equal(dim(assay(es, "child")), c(1L, 3L))
})

testthat::test_that("setSubsetAssay keeps the subset's row and column names", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  genes <- subsetRownames(es, "s1")
  cells <- subsetColnames(es, "s1")
  counts <- as.matrix(assay(es, "s1"))

  # input without names: the subset keeps its names, and so does the assay
  es <- setSubsetAssay(es, "s1", unname(counts), "unnamed")
  testthat::expect_equal(subsetRownames(es, "s1"), genes)
  testthat::expect_equal(subsetColnames(es, "s1"), cells)
  testthat::expect_equal(dimnames(assay(es, "unnamed")), list(genes, cells))

  # input with the same names is accepted
  es <- setSubsetAssay(es, "s1", counts, "named")
  testthat::expect_equal(subsetColnames(es, "s1"), cells)

  # input with reordered or different names is rejected
  reordered <- counts[, c(2, 1, 3, 4)]
  testthat::expect_error(setSubsetAssay(es, "s1", reordered, "reordered"),
                         regexp = "inputMatrix")
  renamed <- counts
  rownames(renamed) <- paste0("other", 1:5)
  testthat::expect_error(setSubsetAssay(es, "s1", renamed, "renamed"),
                         regexp = "inputMatrix")
  testthat::expect_error(assay(es, "s1", subsetAssayName = "x") <- reordered,
                         regexp = "inputMatrix")
  testthat::expect_equal(subsetColnames(es, "s1"), cells)

  # a matrix of the wrong size gets a dimension error, not a names error
  testthat::expect_error(setSubsetAssay(es, "s1", counts[1:3, ], "small"),
                         regexp = "same dimensions")

  # after renaming the subset's cells, its assay round-trips
  subsetColnames(es, "s1") <- paste0("cell", 1:4)
  s1Counts <- assay(es, "s1")
  testthat::expect_equal(colnames(s1Counts), paste0("cell", 1:4))
  es <- setSubsetAssay(es, "s1", s1Counts, "roundTrip")
  testthat::expect_equal(colnames(assay(es, "roundTrip")), paste0("cell", 1:4))
})

testthat::test_that("A subset of a subset round-trips through setSubsetAssay", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "p", rows = 5:8, cols = 3:6)
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2, parentAssay = "p")
  childCounts <- assay(es, "child")
  testthat::expect_equal(rownames(childCounts), rownames(sce_chcl)[5:6])
  testthat::expect_equal(colnames(childCounts), colnames(sce_chcl)[3:4])
  es <- setSubsetAssay(es, "child", childCounts, "copy")
  testthat::expect_equal(dimnames(assay(es, "copy")), dimnames(childCounts))
})

testthat::test_that("A subset of an object without names takes the input's", {
  sce <- SingleCellExperiment(list(counts = matrix(1, 6, 5)))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:3, cols = 1:2)
  counts <- matrix(2, 3, 2, dimnames = list(paste0("g", 1:3), c("a", "b")))
  es <- setSubsetAssay(es, "s1", counts, "named")
  testthat::expect_equal(subsetRownames(es, "s1"), paste0("g", 1:3))
  testthat::expect_equal(subsetColnames(es, "s1"), c("a", "b"))
})

testthat::test_that("Spaces are removed from a new subset's name", {
  es <- ExperimentSubset(sce_chcl)
  spaced <- "my subset"
  testthat::expect_warning(es <- createSubset(es, spaced, 1:5, 1:5),
                           regexp = "Removing spaces")
  testthat::expect_equal(subsetNames(es), "mysubset")
  testthat::expect_equal(subsetDim(es, "mysubset"), c(5L, 5L))
  es <- setSubsetAssay(es, "mysubset", assay(es, "mysubset"), "stored")
  testthat::expect_equal(subsetParent(es, "stored"),
                         list("mysubset", "counts"))
  # the cleaned name can't overwrite an existing subset
  testthat::expect_error(suppressWarnings(createSubset(es, spaced, 1:2, 1:2)),
                         regexp = "already exists")
})

testthat::test_that("Testing createSubset with multiple paramter options",{
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_error(es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "none"), "Input parentAssay does not exist.")
  
  #colnames
  es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109, 200))
  testthat::expect_equal(es@subsets$subset1@colIndices, match(c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), subsetColnames(es, subsetName = "subset1")))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  #null rows and null cols
  es <- createSubset(es, "subset2", parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset2@subsetName, "subset2")
  testthat::expect_equal(es@subsets$subset2@rowIndices, c(1:subsetDim(es, "subset2")[1]))
  testthat::expect_equal(es@subsets$subset2@colIndices, c(1:subsetDim(es, "subset2")[2]))
  testthat::expect_equal(es@subsets$subset2@parentAssay, "counts")
  
  #nas introduced in rows or cols : update: this error is changed and does not appear in this way
  # testthat::expect_error(es <- createSubset(es, "subset22", rows = c("as", "asd", "asd"), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  # testthat::expect_error(es <- createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  
})

testthat::test_that("Completing test coverage, move this into relevant places later",{
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_warning(es <- createSubset(es, "subset 1", rows = c(1:5), cols = (1:5)),
                           "Removing spaces from the specified subsetName.")
  
  es <- ExperimentSubset(list(counts = assay(sce_chcl, "counts")))
  
})