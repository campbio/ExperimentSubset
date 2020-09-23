library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
library(scds)
library(TENxPBMCData)
library(scater)
library(scran)


context("Testing ExperimentSubset functions")

testthat::test_that("Testing ExperimentSubset constructor by manually specifying the data slots",{
  data(sce_chcl, package = "scds")

  es <- ExperimentSubset::ExperimentSubset(
    assays = list(counts = assay(sce_chcl, "counts"),
                  logcounts = assay(sce_chcl, "logcounts")),
    colData=colData(sce_chcl),
    rowData= rowData(sce_chcl))

  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "ExperimentSubset")
})

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SingleCellExperiment object",{
  data(sce_chcl, package = "scds")

  es <- ExperimentSubset::ExperimentSubset(sce_chcl)

  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "ExperimentSubset")
})

testthat::test_that("Testing ExperimentSubset constructor by implicitly providing a SummarizedExperiment object",{
  data(sce_chcl, package = "scds")

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = assay(sce_chcl, "counts"),
                                                                 logcounts = assay(sce_chcl, "logcounts")),
                                                   colData=colData(sce_chcl),
                                                   rowData= rowData(sce_chcl))

  es <- ExperimentSubset::ExperimentSubset(se)

  testthat::expect_true(validObject(es))
  testthat::expect_equal(class(es)[1], "ExperimentSubset")
})

testthat::test_that("Testing SingleCellSubset constructor",{

  testthat::expect_warning(scs <- ExperimentSubset::SingleCellSubset(
    subsetName = "subset 1",
    rowIndices = c(1:10),
    colIndices = c(1:10),
    parentAssay = "counts"
  ), "Removing spaces from subsetName argument.")

  scs <- ExperimentSubset::SingleCellSubset(
    subsetName = "subset1",
    rowIndices = c(1:10),
    colIndices = c(1:10),
    parentAssay = "counts"
  )

  testthat::expect_true(is.character(scs@subsetName))
  testthat::expect_true(is.character(scs@parentAssay))
  testthat::expect_true(is.numeric(scs@rowIndices))
  testthat::expect_true(is.numeric(scs@colIndices))

})

testthat::test_that("Testing createSubset and storeSubset",{
  data(sce_chcl, package = "scds")
  es <- ExperimentSubset::ExperimentSubset(sce_chcl)
  es <- ExperimentSubset::createSubset(es, "subset1",
                                       rows = c(10,11,50,56,98,99,102,105,109,200),
                                       cols = c(20,21,40,45,90,99,100,123,166,299),
                                       parentAssay = "counts")
  testthat::expect_equal(es@subsets$subset1@subsetName, "subset1")
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(10,11,50,56,98,99,102,105,109,200))
  testthat::expect_equal(es@subsets$subset1@rowIndices, c(20,21,40,45,90,99,100,123,166,299))
  testthat::expect_equal(es@subsets$subset1@parentAssay, "counts")
})

# testthat::test_that("Testing createSubset and storeSubset",{
#   data(sce_chcl, package = "scds")
#

#   ExperimentSubset::showSubsetLink(es)
#
#   testthat::expect_true(validObject(es))
#
#   scaledCounts <- ExperimentSubset::assay(es, "subset1")
#   scaledCounts[,] <- scaledCounts[,] + 1
#
#   es <- ExperimentSubset::storeSubset(es, "subset1", scaledCounts, "scaledSubset1")
#
#   ExperimentSubset::rowData(es, subsetName = "scaledSubset1")
#   ExperimentSubset::colData(es, subsetName = "scaledSubset1")
#   ExperimentSubset::rownames(es, subsetName = "scaledSubset1")
#   ExperimentSubset::colnames(es, subsetName = "scaledSubset1")
#
#   es <- ExperimentSubset::storeSubset(es, "subset1", scaledCounts, subsetAssay = NULL)
#
#   testthat::expect_error(es <- ExperimentSubset::storeSubset(es, "subset0", scaledCounts, "scaledSubset2"),
#                          "subset0 does not exist in the subsets slot of the object.")
#
#
#   ExperimentSubset::subsetCount(es)
#   ExperimentSubset::subsetAssayCount(es)
#   ExperimentSubset::subsetAssayNames(es)
#   ExperimentSubset::subsetNames(es)
#   ExperimentSubset::showSubsetLink(es)
#   ExperimentSubset::subsetParent(es, "subset1")
#   ExperimentSubset::show(es)
#   ExperimentSubset::rownames(es, subsetName = "subset1")
#   ExperimentSubset::colnames(es, subsetName = "subset1")
#
#   ExperimentSubset::rownames(es, subsetName = "subset1Internal")
#   ExperimentSubset::colnames(es, subsetName = "subset1Internal")
#
#   testthat::expect_error(ExperimentSubset::colData(es, subsetName = "subset1Internal"),
#                          "Neither a subset nor a subsetAssay.")
#
#   testthat::expect_error(ExperimentSubset::rowData(es, subsetName = "subset1Internal"),
#                          "Neither a subset nor a subsetAssay.")
#
#
#
#   ExperimentSubset::rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[1])))
#   ExperimentSubset::colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[2])))
#   ExperimentSubset::rowData(es, subsetName = "subset1")
#   ExperimentSubset::colData(es, subsetName = "subset1")
#
#
#   es@subsets$subset1@parentAssay = NULL
#   assay2 <- ExperimentSubset::assay(es, "subset1")
#
#   ExperimentSubset::assay(es, "counts2") <- ExperimentSubset::assay(es, "counts")
#
#   testthat::expect_true(validObject(es))
# })
#
# testthat::test_that("Testing 2",{
#   #1 = get pbmc4k dataset and setup es object
#   library(TENxPBMCData)
#   library(ExperimentSubset)
#   tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
#   es <- ExperimentSubset::ExperimentSubset(assays = list(counts = assay(tenx_pbmc4k, "counts")), colData = colData(tenx_pbmc4k), rowData = rowData(tenx_pbmc4k))
#
#   #2 = perform colsums on counts matrix (4340 cells)
#   ExperimentSubset::colData(es) <- cbind(ExperimentSubset::colData(es), colSums = colSums(assay(es, "counts")))
#   ExperimentSubset::rowData(es) <- cbind(ExperimentSubset::rowData(es), rowSums = rowSums(assay(es, "counts")))
#
#   #3 = filter columns with low colsums and create new subset
#   es <- ExperimentSubset::createSubset(es, "filteredAssay", rows = rownames(es), cols = which(colData(es)$colSums > mean(colData(es)$colSums)), parentAssay = "counts")
#
#   #4 = create normalized assay from subset in #3 using scater function
#   ExperimentSubset::assay(es, "filteredAssay", withDimnames = FALSE, subsetAssay = "filteredAssayNormalized") <- scater::normalizeCounts(assay(es, "filteredAssay"))
#
#   #5 = find variable genes and create subset
#   es <- ExperimentSubset::createSubset(es, "hvg10", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 10), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")
#   es <- ExperimentSubset::createSubset(es, "hvg1000", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 1000), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")
#
#   #6 = run pca on hvg, run clustering and store back cluster labels into subset colData
#   #pca
#   ExperimentSubset::reducedDim(es, type = "PCA", subsetName = "hvg1000") <- scater::calculatePCA(assay(es, "hvg1000"))
#
#   ExperimentSubset::reducedDims(es, subsetName = "hvg1000") <- list(a = scater::calculatePCA(assay(es, "hvg1000")), b = scater::calculatePCA(assay(es, "hvg1000")))
#
#   print(ExperimentSubset::reducedDims(es, subsetName = "hvg1000"))
#
#   #clustering
#   set.seed(20)
#   ExperimentSubset::colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
#   ExperimentSubset::reducedDim(es, type = "a", subsetName = "hvg1000")
#   ExperimentSubset::metadata(es, subsetName = "hvg1000") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
#   ExperimentSubset::metadata(es, subsetName = "hvg1000")
#
#   testthat::expect_error(ExperimentSubset::metadata(es, subsetName = "hvg") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2"),
#                          "hvg does not exist in the subsets slot of the object.")
#
#   testthat::expect_error(ExperimentSubset::metadata(es, subsetName = "hvg"),
#                          "hvg does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::metadata(es) <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
#   ExperimentSubset::metadata(es)
#
#   ExperimentSubset::altExp(es, "alt1", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
#   ExperimentSubset::altExp(es, "alt2", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
#   ExperimentSubset::altExp(es, "alt1", subsetName = "hvg1000")
#
#   ExperimentSubset::showSubsetLink(es)
#
#   testthat::expect_error(ExperimentSubset::altExp(es, "alt1", subsetName = "hvg") <- es@subsets$hvg1000@internalAssay,
#                          "hvg does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::altExp(es, subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
#
#   ExperimentSubset::altExp(es, subsetName = "hvg1000")
#
#   ExperimentSubset::altExp(es) <- tenx_pbmc4k
#
#   ExperimentSubset::altExp(es)
#
#   ExperimentSubset::altExp(es, "a1") <- tenx_pbmc4k
#
#   ExperimentSubset::altExp(es, "a1")
#
#   testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "s12121"),
#                          "s12121 does not exist in the subsets slot of the object.")
#
#   testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "s12121") <- c("a123"),
#                          "s12121 does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::altExpNames(es)
#
#   ExperimentSubset::altExpNames(es) <- c("a2")
#
#   testthat::expect_error(ExperimentSubset::altExp(es, "alt1", subsetName = "hvg10002"),
#                          "hvg10002 does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2")
#   ExperimentSubset::altExpNames(es, subsetName = "hvg1000")
#
#   testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "hvg1000x"),
#                          "hvg1000x does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::altExps(es, subsetName = "hvg1000") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay)
#   ExperimentSubset::altExps(es, subsetName = "hvg1000")
#
#   testthat::expect_error(ExperimentSubset::altExps(es, subsetName = "hvg") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay),
#                          "hvg does not exist in the subsets slot of the object.")
#
#   testthat::expect_error(ExperimentSubset::altExps(es, subsetName = "hvg10002"),
#                          "hvg10002 does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::altExps(es) <- list(a = tenx_pbmc4k, b = tenx_pbmc4k)
#
#   ExperimentSubset::altExps(es)
#
#   ExperimentSubset::subsetDim(es, "hvg1000")
#   ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1")
#   ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000")
#
#   ExperimentSubset::showSubsetLink(es)
#
#   testthat::expect_error(ExperimentSubset::reducedDimNames(es, subsetName = "hvg123"),
#                          "hvg123 does not exist in the subsets slot of the object.")
#
#   testthat::expect_error(ExperimentSubset::reducedDimNames(es, subsetName = "hvg123") <- c("PCA_1"),
#                          "hvg123 does not exist in the subsets slot of the object.")
#
#   ExperimentSubset::reducedDim(es, type = "PCA") <- scater::calculatePCA(assay(es, "counts"))
#   ExperimentSubset::reducedDim(es, type = "PCA")
#
#   ExperimentSubset::reducedDims(es) <- list(a = scater::calculatePCA(assay(es, "counts")), b = scater::calculatePCA(assay(es, "counts")))
#   ExperimentSubset::reducedDims(es)
#
#   ExperimentSubset::reducedDimNames(es)
#   ExperimentSubset::reducedDimNames(es) <- c("PCA_Counts")
#
#   ExperimentSubset::assay(es, "counts")
#
#   testthat::expect_true(validObject(es))
# })
#
# testthat::test_that("Testint 3",{
#
#   data(sce_chcl, package = "scds")
#
#   #constructor subset
#   es <- ExperimentSubset::ExperimentSubset(sce_chcl, subset = list(subsetName = "subset10", rows = c(1:10), cols = c(1:2), parentAssay = "counts"))
#
#   #no parent assay
#   es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299))
#
#   testthat::expect_true(validObject(es))
#
# })
#
#
# testthat::test_that("Testint 4",{
#
#   data(sce_chcl, package = "scds")
#
#   es <- ExperimentSubset::ExperimentSubset(sce_chcl)
#
#   #wrong parent assay
#   testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "none"), "Input parentAssay does not exist.")
#
#
#   #colnames
#   es <- ExperimentSubset::createSubset(es, "subset22", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
#   testthat::expect_true(validObject(es))
#
#   #null rows and null cols
#   es <- ExperimentSubset::createSubset(es, "subset223", parentAssay = "counts")
#   testthat::expect_true(validObject(es))
#
#   #nas introduced in rows or cols
#   testthat::expect_error(  es <- ExperimentSubset::createSubset(es, "subset22", rows = c("as", "asd", "asd"), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
#
#
#   testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")
#
#
#
#
#
# })
