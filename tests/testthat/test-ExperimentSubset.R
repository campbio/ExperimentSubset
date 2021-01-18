library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
library(scds)
library(TENxPBMCData)
library(scater)
library(scran)
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
  testthat::expect_equal(ncol(colData(es, subsetName = "s2")), 10)
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
  expect_equal(ncol(rowData(es, subsetName = "subset1")), 2)
  expect_equal(nrow(rowData(es, subsetName = "subset1")), subsetDim(es, subsetName = "subset1")[1])
  expect_equal(ncol(colData(es, subsetName = "subset1")), 12)
  expect_equal(nrow(rowData(es, subsetName = "subset1")), subsetDim(es, subsetName = "subset1")[2])
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
assay(es, "filteredAssay", withDimnames = FALSE, subsetAssayName = "filteredAssayNormalized") <- scater::normalizeCounts(assay(es, "filteredAssay"))
es <- createSubset(es, "hvg1000", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 1000), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")

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