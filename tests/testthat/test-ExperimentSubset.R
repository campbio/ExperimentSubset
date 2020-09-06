library(testthat)
library(SingleCellExperiment)

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


