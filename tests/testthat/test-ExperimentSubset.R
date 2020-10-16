library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
library(scds)
library(TENxPBMCData)
library(scater)
library(scran)

#load data once
data(sce_chcl, package = "scds")
tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")

context("Testing ExperimentSubset functions")

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SingleCellExperiment object",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "ExperimentSubset")
})

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SummarizedExperiment object",{
  se <- SummarizedExperiment::SummarizedExperiment(sce_chcl)
  es <- ExperimentSubset::ExperimentSubset(se)
  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "ExperimentSubset")
})

testthat::test_that("Testing ExperimentSubset constructor with subset & createSubset without parentAssay",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl, subset = list(subsetName = "subset10", rows = c(1:10), cols = c(1:2), parentAssay = "counts"))
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset10")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(1:10))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(1:2))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109, 200))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
})

testthat::test_that("Testing createSubset and storeSubset",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  #createSubset
  es <- ExperimentSubset::createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109,200))
  testthat::expect_equal(es@subsets$subset1@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  es <- ExperimentSubset::createSubset(es, "subset2",
                                       rows = c("OPCML", "RP11-338K17.8", "RP11-770J1.3", "NPSR1", "AC060835.1", "GABRA2", "PRDM16", "CAB39", "LYRM1", "RP1-178F10.3"),
                                       cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
                                       parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset2@subsetName, "subset2")
  testthat::expect_equal(es@subsets$subset2@rowIndices, c(10,11,50,56,98,99,102,105,109,200))
  testthat::expect_equal(es@subsets$subset2@colIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset2@parentAssay, "counts")
  
  expect_error(es <- ExperimentSubset::createSubset(es, "subset3",
                                                    rows = c("OPCML", "RP11-338K17.8", 50,56,98,99,102,105,109,200),
                                                    cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
                                                    parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  
  #storeSubset
  counts1p <- ExperimentSubset::assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  expect_error(es <- ExperimentSubset::storeSubset(es, "subset4", counts1p, "scaledSubset1"),
               "subset4 does not exist in the subsets slot of the object.")
  es <- ExperimentSubset::storeSubset(es, "subset1", counts1p, "scaledSubset1")
  expect_true("scaledSubset1" %in% assayNames(es@subsets$subset1@internalAssay))
  expect_error(es <- ExperimentSubset::storeSubset(es, "subset1", counts1p, subsetAssayName = NULL),
               "subset1 already exists. Please choose a different subsetName parameter.")
  es <- ExperimentSubset::storeSubset(es, "subset4", counts1p, subsetAssayName = NULL)
  expect_equal(rownames(counts1p), rownames(es, "subset4"))
  expect_equal(colnames(counts1p), colnames(es, "subset4"))
})

testthat::test_that("Testing rownames, colnames, rowData and colData",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  es <- ExperimentSubset::createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  counts1p <- ExperimentSubset::assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- ExperimentSubset::storeSubset(es, "subset1", counts1p, "scaledSubset1")
  
  expect_equal(ExperimentSubset::rownames(es, subsetName = "scaledSubset1"), rownames(counts1p))
  expect_equal(ExperimentSubset::colnames(es, subsetName = "scaledSubset1"), colnames(counts1p))
  expect_equal(ExperimentSubset::rownames(es, subsetName = "subset1"), rownames(es)[c(10,11,50,56,98,99,102,105,109,200)])
  expect_equal(ExperimentSubset::colnames(es, subsetName = "subset1"), colnames(es)[c(20,21,40,45,90,99,100,123,166,299)])
  testthat::expect_error(ExperimentSubset::colData(es, subsetName = "subset1Internal"),
                         "Neither a subset nor a subsetAssay.")
  testthat::expect_error(ExperimentSubset::rowData(es, subsetName = "subset1Internal"),
                         "Neither a subset nor a subsetAssay.")
  ExperimentSubset::rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[1])))
  ExperimentSubset::colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[2])))
  expect_equal(ncol(ExperimentSubset::rowData(es, subsetName = "subset1")), 2)
  expect_equal(nrow(ExperimentSubset::rowData(es, subsetName = "subset1")), ExperimentSubset::subsetDim(es, subsetName = "subset1")[1])
  expect_equal(ncol(ExperimentSubset::colData(es, subsetName = "subset1")), 12)
  expect_equal(nrow(ExperimentSubset::rowData(es, subsetName = "subset1")), ExperimentSubset::subsetDim(es, subsetName = "subset1")[2])
})

testthat::test_that("Testing subset helper/supplementary functions",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  es <- ExperimentSubset::createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  counts1p <- ExperimentSubset::assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- ExperimentSubset::storeSubset(es, "subset1", counts1p, "scaledSubset1")
  expect_equal(ExperimentSubset::subsetCount(es), 1)
  expect_equal(ExperimentSubset::subsetAssayCount(es), 2)
  expect_equal(ExperimentSubset::subsetAssayNames(es), c("subset1", "scaledSubset1"))
  expect_equal(ExperimentSubset::subsetNames(es), "subset1")
  ExperimentSubset::subsetSummary(es)
  expect_equal(ExperimentSubset::subsetParent(es, "scaledSubset1"), list("subset1", "counts"))
  ExperimentSubset::show(es)
})

testthat::test_that("Testing assay function",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  es <- ExperimentSubset::createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  es@subsets$subset1@parentAssay <- "counts"
  subset1Assay <- ExperimentSubset::assay(es, "subset1")
  expect_true(all.equal(subset1Assay, ExperimentSubset::assay(es, "counts")[c(10,11,50,56,98,99,102,105,109,200), c(20,21,40,45,90,99,100,123,166,299)]))
  ExperimentSubset::assay(es, "counts2") <- ExperimentSubset::assay(es, "counts")
  expect_true(all.equal(ExperimentSubset::assay(es, "counts"), ExperimentSubset::assay(es, "counts2")))
})

#create es with pbmc4k data, compute colsums and rowsums, create a hvg1000 subset for further testing
es <- ExperimentSubset::ExperimentSubset(tenx_pbmc4k)
ExperimentSubset::colData(es) <- cbind(ExperimentSubset::colData(es), colSums = colSums(assay(es, "counts")))
ExperimentSubset::rowData(es) <- cbind(ExperimentSubset::rowData(es), rowSums = rowSums(assay(es, "counts")))
es <- ExperimentSubset::createSubset(es, "filteredAssay", rows = rownames(es), cols = which(colData(es)$colSums > mean(colData(es)$colSums)), parentAssay = "counts")
ExperimentSubset::assay(es, "filteredAssay", withDimnames = FALSE, subsetAssay = "filteredAssayNormalized") <- scater::normalizeCounts(assay(es, "filteredAssay"))
es <- ExperimentSubset::createSubset(es, "hvg1000", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 1000), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")

testthat::test_that("Testing reducedDims",{
  set.seed(20)
  ExperimentSubset::colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
  expect_true("cluster" %in% colnames(ExperimentSubset::colData(es, subsetName = "hvg1000")))

  ExperimentSubset::reducedDim(es, type = "PCA", subsetName = "hvg1000") <- scater::calculatePCA(assay(es, "hvg1000"))
  expect_true("PCA" %in% ExperimentSubset::reducedDimNames(es, "hvg1000"))
  ExperimentSubset::reducedDims(es, subsetName = "hvg1000") <- list(a = scater::calculatePCA(assay(es, "hvg1000")), b = scater::calculatePCA(assay(es, "hvg1000")))
  expect_true(all.equal(ExperimentSubset::reducedDimNames(es, "hvg1000"), c("a","b")))
  
  ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1", "PCA_2")
  expect_true(all.equal(ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000"), c("PCA_1", "PCA_2")))

  testthat::expect_error(ExperimentSubset::reducedDimNames(es, subsetName = "hvg123"),
                         "hvg123 does not exist in the subsets slot of the object.")
  
  testthat::expect_error(ExperimentSubset::reducedDimNames(es, subsetName = "hvg123") <- c("PCA_1"),
                         "hvg123 does not exist in the subsets slot of the object.")
  
  ExperimentSubset::reducedDim(es, type = "PCA") <- scater::calculatePCA(assay(es, "counts"))
  expect_true(!is.null(ExperimentSubset::reducedDim(es, type = "PCA")))

  ExperimentSubset::reducedDims(es) <- list(a = scater::calculatePCA(assay(es, "counts")), b = scater::calculatePCA(assay(es, "counts")))
  expect_true(!is.null(ExperimentSubset::reducedDims(es)))
  expect_true(all.equal(ExperimentSubset::reducedDimNames(es), c("a", "b")))

  ExperimentSubset::reducedDimNames(es) <- c("PCA_Counts_1", "PCA_Counts_2")
  expect_true(all.equal(ExperimentSubset::reducedDimNames(es), c("PCA_Counts_1", "PCA_Counts_2")))

})

testthat::test_that("Testing metadata",{
  ExperimentSubset::metadata(es, subsetName = "hvg1000") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(ExperimentSubset::metadata(es, subsetName = "hvg1000")$meta1, "This is testing meta1")
  expect_equal(ExperimentSubset::metadata(es, subsetName = "hvg1000")$meta2, "This is testing meta2")
  
  testthat::expect_error(ExperimentSubset::metadata(es, subsetName = "hvg") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2"),
                         "hvg does not exist in the subsets slot of the object.")
  
  testthat::expect_error(ExperimentSubset::metadata(es, subsetName = "hvg"),
                         "hvg does not exist in the subsets slot of the object.")
  
  ExperimentSubset::metadata(es) <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(ExperimentSubset::metadata(es)$meta1, "This is testing meta1")
  expect_equal(ExperimentSubset::metadata(es)$meta2, "This is testing meta2")
})

testthat::test_that("Testing altExperiments",{
  ExperimentSubset::altExp(es, "alt1", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  ExperimentSubset::altExp(es, "alt2", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  expect_true(all.equal(ExperimentSubset::altExpNames(es, "hvg1000"), c("alt1", "alt2")))

  testthat::expect_error(ExperimentSubset::altExp(es, "alt1", subsetName = "hvg") <- es@subsets$hvg1000@internalAssay,
                         "hvg does not exist in the subsets slot of the object.")
  
  ExperimentSubset::altExp(es, subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  expect_equal(length(ExperimentSubset::altExpNames(es, "hvg1000")), 2)
  
  expect_true(!is.null(ExperimentSubset::altExp(es, subsetName = "hvg1000")))
  
  ExperimentSubset::altExp(es) <- tenx_pbmc4k
  expect_true(!is.null(ExperimentSubset::altExp(es)))

  ExperimentSubset::altExp(es, "a1") <- tenx_pbmc4k
  expect_true("a1" %in% ExperimentSubset::altExpNames(es))
  expect_true(!is.null(ExperimentSubset::altExp(es, "a1")))
  
  testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "na"),
                         "na does not exist in the subsets slot of the object.")
  
  testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "na") <- c("na2"),
                         "na does not exist in the subsets slot of the object.")
  
  ExperimentSubset::altExpNames(es) <- c("a1", "a2")
  expect_true(all.equal(ExperimentSubset::altExpNames(es), c("a1", "a2")))
  
  testthat::expect_error(ExperimentSubset::altExp(es, "alt1", subsetName = "hvg10002"),
                         "hvg10002 does not exist in the subsets slot of the object.")
  
  ExperimentSubset::altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2")
  expect_true(all.equal(ExperimentSubset::altExpNames(es, subsetName = "hvg1000"), c("a1", "a2")))

  testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "hvg1000x"),
                         "hvg1000x does not exist in the subsets slot of the object.")
  
  ExperimentSubset::altExps(es, subsetName = "hvg1000") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay)
  expect_true(all.equal(ExperimentSubset::altExpNames(es, "hvg1000"), c("a","b")))

  testthat::expect_error(ExperimentSubset::altExps(es, subsetName = "hvg") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay),
                         "hvg does not exist in the subsets slot of the object.")
  
  testthat::expect_error(ExperimentSubset::altExps(es, subsetName = "hvg10002"),
                         "hvg10002 does not exist in the subsets slot of the object.")
  
  ExperimentSubset::altExps(es) <- list(a = tenx_pbmc4k, b = tenx_pbmc4k)
  expect_true(all.equal(ExperimentSubset::altExpNames(es), c("a", "b")))
} 
)

testthat::test_that("Testing createSubset with multiple paramter options",{
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "none"), "Input parentAssay does not exist.")
  
  #colnames
  es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109, 200))
  testthat::expect_equal(es@subsets$subset1@colIndices, match(c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), ExperimentSubset::colnames(es, subsetName = "subset1")))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
  
  #null rows and null cols
  es <- ExperimentSubset::createSubset(es, "subset2", parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset2@subsetName, "subset2")
  testthat::expect_equal(es@subsets$subset2@rowIndices, c(1:subsetDim(es, "subset2")[1]))
  testthat::expect_equal(es@subsets$subset2@colIndices, c(1:subsetDim(es, "subset2")[2]))
  testthat::expect_equal(es@subsets$subset2@parentAssay, "counts")
  
  #nas introduced in rows or cols
  testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset22", rows = c("as", "asd", "asd"), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
  
})
