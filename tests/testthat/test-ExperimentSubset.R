library(testthat)
context("Testing ExperimentSubset functions")


testthat::test_that("Testing subsetNames",{
  data(sce_chcl, package = "scds")

  es <- ExperimentSubset::ExperimentSubset(
    assays = list(counts = assay(sce_chcl, "counts"),
                  logcounts = assay(sce_chcl, "logcounts")),
    colData=colData(sce_chcl),
    rowData= rowData(sce_chcl))

  subsetNames <- subsetNames(es)

  testthat::expect_null(subsetNames)
})
