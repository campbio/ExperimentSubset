library(testthat)
library(ExperimentSubset)
library(SingleCellExperiment)
library(SummarizedExperiment)

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
  es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "counts")

  testthat::expect_true(validObject(es))

  scaledCounts <- assay(es, "subset1")
  scaledCounts[,] <- scaledCounts[,] + 1
  es <- ExperimentSubset::storeSubset(es, "subset1", scaledCounts, "scaledSubset1")
  ExperimentSubset::subsetCount(es)
  ExperimentSubset::subsetAssayCount(es)
  ExperimentSubset::subsetAssayNames(es)
  ExperimentSubset::subsetNames(es)
  ExperimentSubset::showSubsetLink(es)
  ExperimentSubset::getFirstParent(es, "subset1")
  ExperimentSubset::show(es)
  ExperimentSubset::rownames(es, subsetName = "subset1")
  ExperimentSubset::colnames(es, subsetName = "subset1")
  ExperimentSubset::assay(es, "subset1", newInternalAssay = "subset1Internal") <- ExperimentSubset::assay(es, "subset1")
  ExperimentSubset::rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[1])))
  ExperimentSubset::colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(ExperimentSubset::subsetDim(es, subsetName = "subset1")[2])))
  ExperimentSubset::rowData(es, subsetName = "subset1")
  ExperimentSubset::colData(es, subsetName = "subset1")

  testthat::expect_true(validObject(es))
})

testthat::test_that("Testing 2",{
  #1 = get pbmc4k dataset and setup es object
  library(TENxPBMCData)
  library(ExperimentSubset)
  tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
  es <- ExperimentSubset::ExperimentSubset(assays = list(counts = assay(tenx_pbmc4k, "counts")), colData = colData(tenx_pbmc4k), rowData = rowData(tenx_pbmc4k))

  #2 = perform colsums on counts matrix (4340 cells)
  colData(es) <- cbind(colData(es), colSums = colSums(assay(es, "counts")))

  #3 = filter columns with low colsums and create new subset
  es <- createSubset(es, "filteredAssay", rows = rownames(es), cols = which(colData(es)$colSums > mean(colData(es)$colSums)), parentAssay = "counts")

  #4 = create normalized assay from subset in #3 using scater function
  assay(es, "filteredAssay", newInternalAssay = "filteredAssayNormalized") <- scater::normalizeCounts(assay(es, "filteredAssay"))

  #5 = find variable genes and create subset
  es <- createSubset(es, "hvg10", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 10), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")
  es <- createSubset(es, "hvg1000", rows = scran::getTopHVGs(scran::modelGeneVar(assay(es, "filteredAssayNormalized")), n = 1000), cols = seq(1:1622), parentAssay = "filteredAssayNormalized")

  #6 = run pca on hvg, run clustering and store back cluster labels into subset colData
  #pca
  reducedDim(es, type = "PCA", subsetName = "hvg1000") <- scater::calculatePCA(assay(es, "hvg1000"))
  #clustering
  set.seed(20)
  colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
  reducedDim(es, type = "PCA", subsetName = "hvg1000")
  metadata(es, subsetName = "hvg1000") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  metadata(es, subsetName = "hvg1000")
  ExperimentSubset::altExp(es, "alt1", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  ExperimentSubset::altExp(es, "alt2", subsetName = "hvg1000") <- es@subsets$hvg1000@internalAssay
  ExperimentSubset::altExp(es, "alt1", subsetName = "hvg1000")

  ExperimentSubset::altExp(es, subsetName = "hvg1000")

  ExperimentSubset::altExp(es) <- tenx_pbmc4k

  ExperimentSubset::altExp(es)

  ExperimentSubset::altExp(es, "a1") <- tenx_pbmc4k

  testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "s12121"),
                         "s12121 does not exist in the subsets slot of the object.")

  ExperimentSubset::altExpNames(es)

  testthat::expect_error(ExperimentSubset::altExp(es, "alt1", subsetName = "hvg10002"),
                         "hvg10002 does not exist in the subsets slot of the object.")

  ExperimentSubset::altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2")
  ExperimentSubset::altExpNames(es, subsetName = "hvg1000")

  testthat::expect_error(ExperimentSubset::altExpNames(es, subsetName = "hvg1000x"),
                         "hvg1000x does not exist in the subsets slot of the object.")

  ExperimentSubset::altExps(es, subsetName = "hvg1000") <- list(a = es@subsets$hvg1000@internalAssay, b = es@subsets$hvg1000@internalAssay)
  ExperimentSubset::altExps(es, subsetName = "hvg1000")

  testthat::expect_error(ExperimentSubset::altExps(es, subsetName = "hvg10002"),
                         "hvg10002 does not exist in the subsets slot of the object.")

  ExperimentSubset::altExps(es)

  ExperimentSubset::subsetDim(es, "hvg1000")
  ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1")
  ExperimentSubset::reducedDimNames(es, subsetName = "hvg1000")

  testthat::expect_error(ExperimentSubset::reducedDimNames(es, subsetName = "hvg123"),
                         "hvg123 does not exist in the subsets slot of the object.")

  ExperimentSubset::reducedDimNames(es)

  ExperimentSubset::assay(es, "counts")

  testthat::expect_true(validObject(es))
})

testthat::test_that("Testint 3",{

  data(sce_chcl, package = "scds")

  #constructor subset
  es <- ExperimentSubset::ExperimentSubset(sce_chcl, subset = list(subsetName = "subset10", rows = c(1:10), cols = c(1:2), parentAssay = "counts"))

  #no parent assay
  es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299))

  testthat::expect_true(validObject(es))

})


testthat::test_that("Testint 4",{

  data(sce_chcl, package = "scds")

  es <- ExperimentSubset::ExperimentSubset(sce_chcl)

  #wrong parent assay
  testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "none"), "Input parentAssay does not exist.")


  #colnames
  es <- ExperimentSubset::createSubset(es, "subset22", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
  testthat::expect_true(validObject(es))

  #null rows and null cols
  es <- ExperimentSubset::createSubset(es, "subset223", parentAssay = "counts")
  testthat::expect_true(validObject(es))

  #nas introduced in rows or cols
  testthat::expect_error(  es <- ExperimentSubset::createSubset(es, "subset22", rows = c("as", "asd", "asd"), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")


  testthat::expect_error(es <- ExperimentSubset::createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"), parentAssay = "counts"), "NAs introduced in input rows or columns. Some or all indicated rows or columns not found in specified parent.")





})

