library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
# most tests use sce_chcl from scds (Suggests); skip the file without it
testthat::skip_if_not_installed("scds")

#load data once
data(sce_chcl, package = "scds")
# a small simulated count matrix with dimnames, in place of a download
set.seed(1)
simCounts <- matrix(rpois(300 * 120, lambda = 1), nrow = 300,
                    dimnames = list(paste0("gene", 1:300),
                                    paste0("cell", 1:120)))
simSce <- SingleCellExperiment(list(counts = simCounts))

# Checks that a subset made from the main assay "counts" holds the given
# rows and columns of it, through the exported API
expectSubsetOf <- function(es, subsetName, rows, cols) {
  testthat::expect_equal(subsetRownames(es, subsetName), rownames(es)[rows])
  testthat::expect_equal(subsetColnames(es, subsetName), colnames(es)[cols])
  testthat::expect_equal(subsetParent(es, subsetName), list("counts"))
  testthat::expect_equal(assay(es, subsetName),
                         assay(es, "counts")[rows, cols, drop = FALSE])
}

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
  expectSubsetOf(es, "s1", rows = 1:10, cols = seq_len(ncol(se)))
  testthat::expect_output(show(es), "subsets\\(1\\): s1")
  s1Plus <- assay(es, "s1") + 1
  es <- setSubsetAssay(es, "s1", s1Plus, "newS1")
  testthat::expect_true("s1/newS1" %in% subsetAssayNames(es))
  testthat::expect_false("newS1" %in% subsetNames(es))
  testthat::expect_equal(getSubsetAssay(es, "s1", subsetAssayName = "newS1"),
                         s1Plus)
  testthat::expect_output(show(es), "subsetAssays\\(2\\): s1 s1/newS1")
  assay(es, "assay2") <- assay(es, "counts")
  testthat::expect_true("assay2" %in% assayNames(es))
  testthat::expect_equal(assay(es, "assay2"), assay(es, "counts"))
  testthat::expect_equal(subsetParent(es, "s1"), list("counts"))
  testthat::expect_equal(subsetParent(es, "s1/newS1"), list("s1", "counts"))
  testthat::expect_output(subsetSummary(es), "s1 +10, 2000 +counts +newS1")
  testthat::expect_equal(subsetDim(es, "s1")[1], 10)
  testthat::expect_equal(subsetDim(es, "s1")[2], 2000)
  testthat::expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  testthat::expect_equal(length(metadata(es)), 1)
  testthat::expect_s4_class(getSubsetAssay(es, "s1"), "dgCMatrix")
  testthat::expect_equal(getSubsetAssay(es, "s1"),
                         assay(se, "counts")[1:10, ])
  subsetRownames(es, "s1") <- subsetRownames(es, "s1")
  testthat::expect_equal(subsetRownames(es, "s1"), rownames(se)[1:10])
  subsetColnames(es, "s1") <- subsetColnames(es, "s1")
  testthat::expect_equal(subsetColnames(es, "s1"), colnames(se))
  subsetRownames(es, "s1/newS1") <- subsetRownames(es, "s1/newS1")
  testthat::expect_equal(subsetRownames(es, "s1/newS1"), rownames(se)[1:10])
  subsetColnames(es, "s1/newS1") <- subsetColnames(es, "s1/newS1")
  testthat::expect_equal(subsetColnames(es, "s1/newS1"), colnames(se))
  colData(es, subsetName = "s1") <- cbind(colData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[2]))))
  testthat::expect_equal(ncol(colData(es, subsetName = "s1")), 1)
  testthat::expect_equal(colData(es, subsetName = "s1")$col1, seq_len(2000))
  rowData(es, subsetName = "s1") <- cbind(rowData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[1]))))
  testthat::expect_equal(ncol(rowData(es, subsetName = "s1")), 1)
  testthat::expect_equal(rowData(es, subsetName = "s1")$col1, seq_len(10))
  
  
  testthat::skip_if_not_installed("airway")
  data(airway, package="airway")
  rse <- airway
  testthat::expect_s4_class(es <- ExperimentSubset(rse), "SubsetRangedSummarizedExperiment")
  testthat::expect_output(subsetSummary(es),
                          "Main assay\\(s\\):\n counts.*Subset\\(s\\):\nNULL")
  testthat::expect_equal(length(metadata(es)), 1)

  es <- createSubset(es, "s2", rows = c(1:15))
  testthat::expect_true("s2" %in% subsetNames(es))
  testthat::expect_true("s2" %in% subsetAssayNames(es))
  expectSubsetOf(es, "s2", rows = 1:15, cols = seq_len(ncol(rse)))
  testthat::expect_output(show(es), "subsets\\(1\\): s2")
  s2Plus <- assay(es, "s2") + 1
  es <- setSubsetAssay(es, "s2", s2Plus, "newS2")
  testthat::expect_true("s2/newS2" %in% subsetAssayNames(es))
  testthat::expect_false("newS2" %in% subsetNames(es))
  testthat::expect_equal(getSubsetAssay(es, "s2", subsetAssayName = "newS2"),
                         s2Plus)
  assay(es, "assay2") <- assay(es, "counts")
  testthat::expect_true("assay2" %in% assayNames(es))
  testthat::expect_equal(assay(es, "assay2"), assay(es, "counts"))
  testthat::expect_equal(subsetParent(es, "s2"), list("counts"))
  testthat::expect_output(subsetSummary(es), "s2 +15, 8 +counts +newS2")
  testthat::expect_equal(subsetDim(es, "s2")[1], 15)
  testthat::expect_equal(subsetDim(es, "s2")[2], 8)
  testthat::expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  testthat::expect_equal(length(metadata(es)), 1)

  subsetRownames(es, "s2") <- subsetRownames(es, "s2")
  testthat::expect_equal(subsetRownames(es, "s2"), rownames(rse)[1:15])
  subsetColnames(es, "s2") <- subsetColnames(es, "s2")
  testthat::expect_equal(subsetColnames(es, "s2"), colnames(rse))
  subsetRownames(es, "s2/newS2") <- subsetRownames(es, "s2/newS2")
  testthat::expect_equal(subsetRownames(es, "s2/newS2"),
                         rownames(rse)[1:15])
  subsetColnames(es, "s2/newS2") <- subsetColnames(es, "s2/newS2")
  testthat::expect_equal(subsetColnames(es, "s2/newS2"), colnames(rse))
  colData(es, subsetName = "s2") <- cbind(colData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[2]))))
  testthat::expect_equal(ncol(colData(es, subsetName = "s2", parentColData = TRUE)), 10)
  testthat::expect_equal(colData(es, subsetName = "s2")$col1, seq_len(8))
  rowData(es, subsetName = "s2") <- cbind(rowData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[1]))))
  testthat::expect_equal(ncol(rowData(es, subsetName = "s2")), 1)
  testthat::expect_equal(rowData(es, subsetName = "s2")$col1, seq_len(15))
  
  # compare the whole message, so it is one line with no runs of spaces
  testthat::expect_equal(
    tryCatch(createSubset(es, "s2", rows = c(1:15)),
             error = conditionMessage),
    paste0("'subsetName' must differ from the names of existing subsets; ",
           "'s2' is already a subset name.")
  )
  testthat::expect_error(es <- createSubset(es, subsetName = 12, rows = c(1:15)),
                         "'subsetName' must be a single, non-empty")
  testthat::expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", rows = list()),
             error = conditionMessage),
    paste0("'rows' must be NULL, to include all rows, or a numeric or ",
           "character vector of the rows to include.")
  )
  testthat::expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", cols = list()),
             error = conditionMessage),
    paste0("'cols' must be NULL, to include all columns, or a numeric or ",
           "character vector of the columns to include.")
  )
  testthat::expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", parentAssay = 1),
             error = conditionMessage),
    paste0("'parentAssay' must be NULL, to use the default assay, or a ",
           "single character string naming the assay or subset to make ",
           "the subset from.")
  )
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
  testthat::expect_equal(subsetNames(es), "subset10")
  expectSubsetOf(es, "subset10", rows = 1:10, cols = 1:2)
  
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1)
  testthat::expect_equal(subsetNames(es), c("subset10", "subset1"))
  expectSubsetOf(es, "subset1", rows = rows1, cols = cols1)
  expectSubsetOf(es, "subset10", rows = 1:10, cols = 1:2)
})

testthat::test_that("Testing createSubset and setSubsetAssay",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  #createSubset
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  expectSubsetOf(es, "subset1", rows = rows1, cols = cols1)
  
  subsetRownames(es, "subset1") <- subsetRownames(es, "subset1")
  subsetColnames(es, "subset1") <- subsetColnames(es, "subset1")
  testthat::expect_equal(subsetRownames(es, "subset1"), rownames(es)[rows1])
  testthat::expect_equal(subsetColnames(es, "subset1"), colnames(es)[cols1])
  
  # the same rows and columns, given by name
  es <- createSubset(es, "subset2",
                                       rows = c("OPCML", "RP11-338K17.8", "RP11-770J1.3", "NPSR1", "AC060835.1", "GABRA2", "PRDM16", "CAB39", "LYRM1", "RP1-178F10.3"),
                                       cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
                                       parentAssay = "counts")
  expectSubsetOf(es, "subset2", rows = rows1, cols = cols1)
  
  # numbers given as names are not row names
  testthat::expect_error(
    createSubset(es, "subset3", rows = c("OPCML", "RP11-338K17.8", "50"),
                 cols = 1:2, parentAssay = "counts"),
    regexp = "Selected rows not available in the specified assay."
  )
  
  #setSubsetAssay
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  expect_error(es <- setSubsetAssay(es, "subset4", counts1p, "scaledSubset1"),
               "subset4 does not exist in the subsets slot of the object.")
  # the whole message is one line with no runs of spaces
  testthat::expect_equal(
    tryCatch(setSubsetAssay(es, "subset4", counts1p, "scaledSubset1"),
             error = conditionMessage),
    paste0("subset4 does not exist in the subsets slot of the object. ",
           "You need to create a new subset with createSubset() function.")
  )
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  testthat::expect_true("subset1/scaledSubset1" %in% subsetAssayNames(es))
  testthat::expect_equal(
    getSubsetAssay(es, "subset1", subsetAssayName = "scaledSubset1"), counts1p
  )
  # the subset's own data is unchanged
  testthat::expect_equal(getSubsetAssay(es, "subset1"),
                         assay(es, "counts")[rows1, cols1])
  # an assay stored in a subset needs a name
  # NULL, NA and 1L failed at dispatch before reaching this check (#97)
  for (badName in list(NA_character_, "", c("a", "b"), NULL, NA, 1L)) {
    testthat::expect_error(
      setSubsetAssay(es, "subset1", counts1p, subsetAssayName = badName),
      regexp = "'subsetAssayName' must be a single, non-empty character",
      info = deparse(badName)
    )
  }
})

testthat::test_that("Testing rownames, colnames, rowData and colData",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  
  testthat::expect_equal(subsetRownames(es, "subset1/scaledSubset1"),
                         rownames(es)[rows1])
  testthat::expect_equal(subsetColnames(es, "subset1/scaledSubset1"),
                         colnames(es)[cols1])
  testthat::expect_equal(subsetRownames(es, "subset1"), rownames(es)[rows1])
  testthat::expect_equal(subsetColnames(es, "subset1"), colnames(es)[cols1])
  testthat::expect_error(colData(es, subsetName = "subset1Internal"),
                         "subset1Internal subset does not exist.")
  testthat::expect_error(rowData(es, subsetName = "subset1Internal"),
                         "subset1Internal subset does not exist.")
  rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[1])))
  colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[2])))
  expect_equal(ncol(rowData(es, subsetName = "subset1", parentRowData = TRUE)), 2)
  expect_equal(nrow(rowData(es, subsetName = "subset1", parentRowData = TRUE)), subsetDim(es, subsetName = "subset1")[1])
  expect_equal(ncol(colData(es, subsetName = "subset1", parentColData = TRUE)), 12)
  expect_equal(nrow(colData(es, subsetName = "subset1", parentColData = TRUE)),
               subsetDim(es, subsetName = "subset1")[2])
  testthat::expect_equal(rowData(es, subsetName = "subset1")$col1,
                         seq_along(rows1))
  testthat::expect_equal(colData(es, subsetName = "subset1")$col1,
                         seq_along(cols1))
  # the subset is square, so check names to tell its rows from its columns
  testthat::expect_equal(
    rownames(rowData(es, subsetName = "subset1", parentRowData = TRUE)),
    rownames(es)[rows1]
  )
  testthat::expect_equal(
    rownames(colData(es, subsetName = "subset1", parentColData = TRUE)),
    colnames(es)[cols1]
  )
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
  expect_equal(subsetAssayNames(es), c("subset1", "subset1/scaledSubset1"))
  expect_equal(subsetNames(es), "subset1")
  testthat::expect_output(subsetSummary(es),
                          "subset1 +10, 10 +counts +scaledSubset1")
  expect_equal(subsetParent(es, "subset1/scaledSubset1"),
               list("subset1", "counts"))
  testthat::expect_output(show(es),
                          "subsetAssays\\(2\\): subset1 subset1/scaledSubset1")
})

testthat::test_that("Testing assay function",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  subset1Assay <- assay(es, "subset1")
  testthat::expect_equal(subset1Assay, assay(es, "counts")[rows1, cols1])
  assay(es, "counts2") <- assay(es, "counts")
  testthat::expect_equal(assay(es, "counts2"), assay(es, "counts"))
})

# Builds the object the reducedDims, metadata and altExp tests share, as
# the vignette does: keep cells with more counts than average, normalize
# and choose highly variable genes with scrapper, then subset those genes
# by name. The subset keeps its historical name, hvg1000.
makeHvgFixture <- function() {
  es <- ExperimentSubset(simSce)
  colData(es) <- cbind(colData(es), colSums = colSums(assay(es, "counts")))
  rowData(es) <- cbind(rowData(es), rowSums = rowSums(assay(es, "counts")))
  keptCells <- which(colData(es)$colSums > mean(colData(es)$colSums))
  es <- createSubset(es, "filteredAssay", rows = rownames(es),
                     cols = keptCells, parentAssay = "counts")
  filteredCounts <- assay(es, "filteredAssay")
  filteredSizeFactors <- scrapper::centerSizeFactors(colSums(filteredCounts))
  assay(
    es, "filteredAssay", withDimnames = FALSE,
    subsetAssayName = "filteredAssayNormalized"
  ) <- scrapper::normalizeCounts(
    filteredCounts, size.factors = filteredSizeFactors
  )
  geneVariances <- scrapper::modelGeneVariances(
    assay(es, "filteredAssay/filteredAssayNormalized")
  )
  hvgIndices <- scrapper::chooseHighlyVariableGenes(
    geneVariances$statistics$residuals, top = 50
  )
  #pass gene names to keep character rows with a subset parent covered
  normalized <- assay(es, "filteredAssay/filteredAssayNormalized")
  hvgNames <- rownames(normalized)[hvgIndices]
  stopifnot(length(hvgNames) > 0, !anyNA(hvgNames))
  createSubset(
    es, "hvg1000", rows = hvgNames, cols = seq_along(keptCells),
    parentAssay = "filteredAssay/filteredAssayNormalized"
  )
}

testthat::test_that("Testing reducedDims",{
  testthat::skip_if_not_installed("scrapper")
  testthat::skip_if_not_installed("scater")
  es <- makeHvgFixture()
  set.seed(20)
  colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
  expect_true("cluster" %in% colnames(colData(es, subsetName = "hvg1000")))

  hvgPca <- scater::calculatePCA(assay(es, "hvg1000"))
  reducedDim(es, "PCA", subsetName = "hvg1000") <- hvgPca
  expect_true("PCA" %in% reducedDimNames(es, subsetName = "hvg1000"))
  pcaList <- list(a = hvgPca, b = hvgPca)
  reducedDims(es, subsetName = "hvg1000") <- pcaList
  testthat::expect_equal(reducedDimNames(es, subsetName = "hvg1000"),
                         c("a", "b"))
  testthat::expect_equal(reducedDim(es, "b", subsetName = "hvg1000"), hvgPca)
  
  reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1", "PCA_2")
  testthat::expect_equal(reducedDimNames(es, subsetName = "hvg1000"),
                         c("PCA_1", "PCA_2"))
  testthat::expect_equal(reducedDims(es, subsetName = "hvg1000"),
                         SimpleList(PCA_1 = hvgPca, PCA_2 = hvgPca))

  testthat::expect_error(reducedDimNames(es, subsetName = "hvg123"),
                         "hvg123 subset does not exist.")
  
  testthat::expect_error(reducedDimNames(es, subsetName = "hvg123") <- c("PCA_1"),
                         "hvg123 subset does not exist.")
  
  countsPca <- scater::calculatePCA(assay(es, "counts"))
  reducedDim(es, type = "PCA") <- countsPca
  testthat::expect_equal(reducedDim(es, type = "PCA"), countsPca)

  reducedDims(es) <- list(a = countsPca, b = countsPca * 2)
  testthat::expect_equal(reducedDims(es),
                         SimpleList(a = countsPca, b = countsPca * 2))
  testthat::expect_equal(reducedDimNames(es), c("a", "b"))

  reducedDimNames(es) <- c("PCA_Counts_1", "PCA_Counts_2")
  testthat::expect_equal(reducedDimNames(es), c("PCA_Counts_1", "PCA_Counts_2"))
  testthat::expect_equal(reducedDim(es, "PCA_Counts_2"), countsPca * 2)

})

testthat::test_that("Testing metadata",{
  testthat::skip_if_not_installed("scrapper")
  es <- makeHvgFixture()
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
  testthat::skip_if_not_installed("scrapper")
  es <- makeHvgFixture()
  # altExps for the subset's cells, told apart by their counts
  hvgCells <- subsetColnames(es, "hvg1000")
  makeHvgAlt <- function(value) {
    counts <- matrix(value, 2, length(hvgCells),
                     dimnames = list(c("f1", "f2"), hvgCells))
    SingleCellExperiment(list(counts = counts))
  }
  hvgAlt <- makeHvgAlt(1)
  altExp(es, "alt1", subsetName = "hvg1000") <- hvgAlt
  altExp(es, "alt2", subsetName = "hvg1000") <- hvgAlt
  testthat::expect_equal(altExpNames(es, subsetName = "hvg1000"),
                         c("alt1", "alt2"))

  testthat::expect_error(altExp(es, "alt1", subsetName = "hvg") <- hvgAlt,
                         "hvg subset does not exist.")
  
  # without a name, altExp<- replaces the first altExp
  altExp(es, subsetName = "hvg1000") <- makeHvgAlt(2)
  expect_equal(length(altExpNames(es, subsetName = "hvg1000")), 2)
  
  testthat::expect_equal(
    assay(altExp(es, subsetName = "hvg1000"), "counts"),
    assay(makeHvgAlt(2), "counts")
  )
  testthat::expect_equal(
    assay(altExp(es, "alt2", subsetName = "hvg1000"), "counts"),
    assay(hvgAlt, "counts")
  )
  
  altExp(es) <- simSce
  testthat::expect_equal(assay(altExp(es), "counts"), simCounts)

  altExp(es, "a1") <- simSce
  expect_true("a1" %in% altExpNames(es))
  testthat::expect_equal(assay(altExp(es, "a1"), "counts"), simCounts)
  
  testthat::expect_error(altExpNames(es, subsetName = "na"),
                         "na subset does not exist.")
  
  testthat::expect_error(altExpNames(es, subsetName = "na") <- c("na2"),
                         "na subset does not exist.")
  
  altExpNames(es) <- c("a1", "a2")
  testthat::expect_equal(altExpNames(es), c("a1", "a2"))
  
  testthat::expect_error(altExp(es, "alt1", subsetName = "hvg10002"),
                         "hvg10002 subset does not exist.")
  
  # the subset sees the main object's altExps a1 and a2, then its own
  # alt1 and alt2, and renames only its own (#58)
  altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2", "b1", "b2")
  testthat::expect_equal(altExpNames(es, subsetName = "hvg1000"),
                         c("a1", "a2", "b1", "b2"))

  testthat::expect_error(altExpNames(es, subsetName = "hvg1000x"),
                         "hvg1000x subset does not exist.")
  
  altExps(es, subsetName = "hvg1000") <- list(a = hvgAlt, b = makeHvgAlt(3))
  # its own altExps are replaced; the main object's stay visible (#58)
  testthat::expect_equal(altExpNames(es, subsetName = "hvg1000"),
                         c("a1", "a2", "a", "b"))
  testthat::expect_equal(
    assay(altExp(es, "b", subsetName = "hvg1000"), "counts"),
    assay(makeHvgAlt(3), "counts")
  )

  testthat::expect_error(altExps(es, subsetName = "hvg") <- list(a = hvgAlt),
                         "hvg subset does not exist.")
  
  testthat::expect_error(altExps(es, subsetName = "hvg10002"),
                         "hvg10002 subset does not exist.")
  
  altExps(es) <- list(a = simSce, b = simSce)
  testthat::expect_equal(altExpNames(es), c("a", "b"))
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
  es <- createSubset(es, "s4", rows = 5:6, cols = 2:3,
                     parentAssay = "s1/s1copy")

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

testthat::test_that("Subset rowData and colData keep the subset's names", {
  sce <- sce_chcl
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  genes <- rownames(sce)[1:5]
  cells <- colnames(sce)[1:8]

  # getters return the subset's row and column names
  testthat::expect_equal(rownames(rowData(es, subsetName = "s1")), genes)
  testthat::expect_equal(rownames(colData(es, subsetName = "s1")), cells)
  s1ColData <- subsetColData(es, "s1", parentColData = TRUE)
  testthat::expect_equal(rownames(s1ColData), cells)
  testthat::expect_equal(s1ColData$cellIndex, 1:8)

  # get, modify and set back keeps the cell names
  s1ColData <- subsetColData(es, "s1", parentColData = FALSE)
  s1ColData$score <- seq_len(8)
  subsetColData(es, "s1") <- s1ColData
  testthat::expect_equal(subsetColnames(es, "s1"), cells)
  testthat::expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  # so a reducedDim with the cell names can be stored without a warning
  umap <- matrix(0, 8, 2, dimnames = list(cells, c("U1", "U2")))
  testthat::expect_warning(reducedDim(es, "umap", subsetName = "s1") <- umap,
                           NA)

  # a value without rownames takes the subset's names
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:8])
  testthat::expect_equal(subsetColnames(es, "s1"), cells)
  rowData(es, subsetName = "s1") <- DataFrame(rank = 5:1)
  testthat::expect_equal(subsetRownames(es, "s1"), genes)

  # colData with other rownames, or of the wrong size, is an error
  otherCells <- DataFrame(group = letters[1:8], row.names = LETTERS[1:8])
  testthat::expect_error(subsetColData(es, "s1") <- otherCells,
                         regexp = "rownames of 'value'")
  testthat::expect_error(colData(es, subsetName = "s1") <- otherCells,
                         regexp = "rownames of 'value'")
  testthat::expect_error(subsetColData(es, "s1") <- DataFrame(a = 1:5),
                         regexp = "nrow")
  testthat::expect_equal(subsetColnames(es, "s1"), cells)
  # rowData with other rownames is an error too, so values meant for other
  # rows aren't stored by position (#100)
  rowData(es, subsetName = "s1") <- DataFrame(rank = 11:15, row.names = genes)
  testthat::expect_equal(rowData(es, subsetName = "s1")$rank, 11:15)
  rowMsg <- paste0("The rownames of 'value' must be NULL or the row names ",
                   "of subset 's1', in the same order.")
  otherGenes <- DataFrame(rank = 1:5, row.names = LETTERS[1:5])
  testthat::expect_error(subsetRowData(es, "s1") <- otherGenes,
                         regexp = rowMsg, fixed = TRUE)
  testthat::expect_error(rowData(es, subsetName = "s1") <- otherGenes,
                         regexp = rowMsg, fixed = TRUE)
  reordered <- DataFrame(rank = 1:5, row.names = rev(genes))
  testthat::expect_error(subsetRowData(es, "s1") <- reordered,
                         regexp = rowMsg, fixed = TRUE)
  otherFrame <- data.frame(rank = 1:5, row.names = paste0("x", 1:5))
  testthat::expect_error(subsetRowData(es, "s1") <- otherFrame,
                         regexp = rowMsg, fixed = TRUE)
  testthat::expect_equal(subsetRownames(es, "s1"), genes)
  testthat::expect_equal(rowData(es, subsetName = "s1")$rank, 11:15)
  # a value of the wrong size gets the size error, even with rownames
  tooFewGenes <- DataFrame(rank = 1:3, row.names = genes[1:3])
  tooFewCells <- DataFrame(group = letters[1:5], row.names = cells[1:5])
  testthat::expect_error(subsetRowData(es, "s1") <- tooFewGenes,
                         regexp = "3 rows in value to replace 5")
  testthat::expect_error(subsetColData(es, "s1") <- tooFewCells,
                         regexp = "nrow")
  # the colData message names the subset's column names
  testthat::expect_error(subsetColData(es, "s1") <- otherCells,
                         regexp = paste0("The rownames of 'value' must be ",
                                         "NULL or the column names of ",
                                         "subset 's1', in the same order."),
                         fixed = TRUE)
  # a subset without row names has none to check against
  unnamed <- ExperimentSubset(list(counts = matrix(1:20, nrow = 5)))
  unnamed <- createSubset(unnamed, "s1", rows = 1:3, cols = 1:2)
  rowData(unnamed, subsetName = "s1") <- otherGenes[1:3, , drop = FALSE]
  testthat::expect_equal(rowData(unnamed, subsetName = "s1")$rank, 1:3)
  testthat::expect_null(subsetRownames(unnamed, "s1"))

  # with the parent's data, the result has the subset's renamed names
  subsetColnames(es, "s1") <- paste0("cell", 1:8)
  subsetRownames(es, "s1") <- paste0("gene", 1:5)
  s1ColData <- subsetColData(es, "s1", parentColData = TRUE)
  s1RowData <- subsetRowData(es, "s1", parentRowData = TRUE)
  testthat::expect_equal(rownames(s1ColData), paste0("cell", 1:8))
  testthat::expect_equal(rownames(s1RowData), paste0("gene", 1:5))
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
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
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

testthat::test_that("Every Subset* class is an ExperimentSubset (#5)", {
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  se <- SummarizedExperiment(assays = list(counts = counts))
  rse <- as(se, "RangedSummarizedExperiment")
  sce <- SingleCellExperiment(assays = list(counts = counts))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts))
  spe <- SpatialExperiment(assays = list(counts = counts))
  inputs <- list(se = se, rse = rse, sce = sce, tse = tse, spe = spe,
                 list = list(counts = counts))
  expectedClasses <- c(se = "SubsetSummarizedExperiment",
                       rse = "SubsetRangedSummarizedExperiment",
                       sce = "SubsetSingleCellExperiment",
                       tse = "SubsetTreeSummarizedExperiment",
                       spe = "SubsetSpatialExperiment",
                       list = "SubsetSingleCellExperiment")
  for (name in names(inputs)) {
    es <- ExperimentSubset(inputs[[name]])
    testthat::expect_equal(class(es)[1], expectedClasses[[name]],
                           info = name)
    testthat::expect_true(is(es, "ExperimentSubset"), info = name)
  }

  # the issue's case: a SubsetSpatialExperiment isn't a
  # SubsetSingleCellExperiment, but it is an ExperimentSubset
  testthat::expect_false(is(ExperimentSubset(spe),
                            "SubsetSingleCellExperiment"))
  testthat::expect_false(is(sce, "ExperimentSubset"))
  testthat::expect_false(is(se, "ExperimentSubset"))
  testthat::expect_false(is(spe, "ExperimentSubset"))

  # a new Subset* class must be added to the union in allClasses.R
  packageClasses <- getClasses(asNamespace("ExperimentSubset"))
  subsetClasses <- grep("^Subset", packageClasses, value = TRUE)
  testthat::expect_length(subsetClasses, 5)
  for (subsetClass in subsetClasses) {
    testthat::expect_true(extends(subsetClass, "ExperimentSubset"),
                          info = subsetClass)
  }
})

# Runs the methods every Subset*Experiment class shares on es, whose main
# assay "counts" has row and column names and at least 5 rows and 6
# columns. Returns es with a subset "s1" and a stored assay "doubled".
checkSharedMethods <- function(es) {
  counts <- as.matrix(assay(es, "counts"))
  es <- createSubset(es, "s1", rows = 2:5, cols = 3:6)
  testthat::expect_equal(subsetNames(es), "s1")
  testthat::expect_equal(subsetCount(es), 1L)
  testthat::expect_equal(subsetDim(es, "s1"), c(4L, 4L))
  testthat::expect_equal(subsetRownames(es, "s1"), rownames(counts)[2:5])
  testthat::expect_equal(subsetColnames(es, "s1"), colnames(counts)[3:6])
  testthat::expect_equal(as.matrix(assay(es, "s1")), counts[2:5, 3:6])

  doubled <- counts[2:5, 3:6] * 2
  es <- setSubsetAssay(es, "s1", doubled, "doubled")
  storedDoubled <- getSubsetAssay(es, "s1", subsetAssayName = "doubled")
  testthat::expect_equal(as.matrix(storedDoubled), doubled)
  testthat::expect_equal(subsetAssayNames(es), c("s1", "s1/doubled"))
  testthat::expect_equal(subsetAssayCount(es), 2L)
  testthat::expect_equal(subsetParent(es, "s1/doubled"), list("s1", "counts"))
  testthat::expect_output(show(es), "subsets\\(1\\): s1")
  testthat::expect_output(subsetSummary(es), "doubled")

  metadata(es, subsetName = "s1") <- list(note = "kept")
  testthat::expect_equal(metadata(es, subsetName = "s1")$note, "kept")

  genes <- subsetRownames(es, "s1")
  cells <- subsetColnames(es, "s1")
  rowData(es, subsetName = "s1") <- DataFrame(score = 1:4, row.names = genes)
  testthat::expect_equal(rowData(es, subsetName = "s1")$score, 1:4)
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4],
                                              row.names = cells)
  testthat::expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  subsetRowData(es, "s1") <- DataFrame(rank = 4:1, row.names = genes)
  testthat::expect_equal(subsetRowData(es, "s1", FALSE)$rank, 4:1)
  subsetColData(es, "s1") <- DataFrame(batch = 1:4, row.names = cells)
  testthat::expect_equal(subsetColData(es, "s1", FALSE)$batch, 1:4)
  parentCols <- subsetColData(es, "s1", parentColData = TRUE)
  testthat::expect_equal(rownames(parentCols), colnames(counts)[3:6])

  # reducedDims and altExps, which these classes get from SCE
  pca <- matrix(as.numeric(1:8), nrow = 4, dimnames = list(cells, NULL))
  reducedDim(es, "PCA", subsetName = "s1") <- pca
  testthat::expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca)
  testthat::expect_equal(reducedDimNames(es, subsetName = "s1"), "PCA")
  testthat::expect_equal(names(reducedDims(es, subsetName = "s1")), "PCA")
  altCounts <- matrix(1, 2, 4, dimnames = list(NULL, cells))
  alt <- SingleCellExperiment(list(counts = altCounts))
  altExp(es, "alt", subsetName = "s1") <- alt
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), "alt")
  testthat::expect_equal(dim(altExp(es, "alt", subsetName = "s1")), c(2L, 4L))
  testthat::expect_equal(names(altExps(es, subsetName = "s1")), "alt")

  # access by position or default, and the remaining setters
  testthat::expect_equal(reducedDim(es, 1, subsetName = "s1"), pca)
  reducedDim(es, 1, subsetName = "s1") <- pca * 2
  testthat::expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca * 2)
  testthat::expect_equal(dim(altExp(es, 1, subsetName = "s1")), c(2L, 4L))
  testthat::expect_equal(dim(altExp(es, subsetName = "s1")), c(2L, 4L))
  altExp(es, subsetName = "s1") <- alt
  altExp(es, 1, subsetName = "s1") <- alt
  altExps(es, subsetName = "s1") <- list(first = alt, second = alt)
  testthat::expect_equal(altExpNames(es, subsetName = "s1"),
                         c("first", "second"))
  altExpNames(es, subsetName = "s1") <- c("x1", "x2")
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), c("x1", "x2"))
  tripled <- counts[2:5, 3:6] * 3
  suppressMessages(assay(es, "s1", subsetAssayName = "tripled") <- tripled)
  testthat::expect_equal(as.matrix(assay(es, "s1/tripled")), tripled)

  # renaming comes last, since later lookups would see the new names
  subsetRownames(es, "s1") <- paste0("g", 1:4)
  subsetColnames(es, "s1") <- paste0("c", 1:4)
  testthat::expect_equal(subsetRownames(es, "s1"), paste0("g", 1:4))
  testthat::expect_equal(subsetColnames(es, "s1"), paste0("c", 1:4))
  es
}

testthat::test_that("TreeSummarizedExperiment subsets support all methods", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  testthat::expect_s4_class(es, "SubsetTreeSummarizedExperiment")
  es <- checkSharedMethods(es)

  links <- subsetRowLinks(es, "s1", parentRowLinkData = TRUE)
  testthat::expect_equal(links$nodeLab, rowLinks(tse)$nodeLab[2:5])
  # the subset itself stores no links
  testthat::expect_null(subsetRowLinks(es, "s1", parentRowLinkData = FALSE))
  testthat::expect_null(subsetColLinks(es, "s1", parentColLinkData = FALSE))
  # parentColLinkData = TRUE errors without a column tree (#44)
})

testthat::test_that("SpatialExperiment subsets support all methods", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  coords <- cbind(x = as.numeric(1:6), y = as.numeric(1:6) * 10)
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = coords)
  es <- ExperimentSubset(spe)
  testthat::expect_s4_class(es, "SubsetSpatialExperiment")
  es <- checkSharedMethods(es)

  spotCoords <- spatialCoords(es, subsetName = "s1")
  testthat::expect_equal(unname(spotCoords[, "y"]), c(30, 40, 50, 60))
  testthat::expect_equal(subsetSpatialCoords(es, "s1"), spotCoords)
  # the main object's coordinates are unchanged
  testthat::expect_equal(spatialCoords(es), spatialCoords(spe))
})

testthat::test_that("Subset spatialData is deprecated, its setter defunct", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  # SpatialExperiment's own deprecation note is a message
  suppressMessages(spatialData(spe) <- DataFrame(region = letters[1:6]))
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = 2:5)

  # the getters warn, then return the main object's columns as before
  testthat::expect_warning(
    s1Data <- suppressMessages(
      subsetSpatialData(es, "s1", parentSpatialData = TRUE)
    ),
    "deprecated"
  )
  testthat::expect_equal(s1Data$region, letters[2:5])
  testthat::expect_warning(
    s1Data <- suppressMessages(
      spatialData(es, subsetName = "s1", parentSpatialData = TRUE)
    ),
    "deprecated"
  )
  testthat::expect_equal(s1Data$region, letters[2:5])
  # the main object's spatialData is SpatialExperiment's, with no warning
  testthat::expect_warning(suppressMessages(spatialData(es)), NA)
  # one warning, without SpatialExperiment's own notes as well
  testthat::expect_message(
    suppressWarnings(subsetSpatialData(es, "s1", parentSpatialData = TRUE)),
    NA
  )
  # a missing subset is an error, before any deprecation warning
  testthat::expect_warning(
    testthat::expect_error(subsetSpatialData(es, "s9"), "does not exist"),
    NA
  )

  # the setter fails clearly, for both kinds of data frame
  testthat::expect_error(
    subsetSpatialData(es, "s1") <- data.frame(region = letters[1:4]),
    "defunct"
  )
  testthat::expect_error(
    subsetSpatialData(es, "s1") <- DataFrame(region = letters[1:4]),
    "defunct"
  )
  # and for a subset given by number
  testthat::expect_error(
    subsetSpatialData(es, 1) <- data.frame(region = letters[1:4]),
    "defunct"
  )
})

testthat::test_that("reducedDim with subsetName works without a type", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  pca <- matrix(1, nrow = 4, ncol = 2)
  reducedDim(es, "PCA", subsetName = "s1") <- pca
  testthat::expect_equal(dim(reducedDim(es, subsetName = "s1")), c(4L, 2L))
  # without a type, the setter replaces the first reducedDim
  reducedDim(es, subsetName = "s1") <- pca * 2
  testthat::expect_equal(unname(reducedDim(es, "PCA", subsetName = "s1")),
                         pca * 2)
})

testthat::test_that("reducedDim on the main object works without a type", {
  # (#77) each class with reducedDims, getter and setter, no subsetName
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  objects <- list(
    sce = SingleCellExperiment(list(counts = counts)),
    spe = SpatialExperiment(assays = list(counts = counts),
                            spatialCoords = cbind(x = 1:6, y = 1:6)),
    tse = TreeSummarizedExperiment(assays = list(counts = counts),
                                   rowTree = tinyTree)
  )
  pca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    testthat::expect_s4_class(es, paste0("Subset", class(objects[[name]])))
    reducedDim(es, "PCA") <- pca
    reducedDim(es, "UMAP") <- pca * 10
    # without a type, the getter returns the first reducedDim
    testthat::expect_equal(unname(reducedDim(es)), pca, info = name)
    # without a type, the setter replaces the first reducedDim
    reducedDim(es) <- pca * 2
    testthat::expect_equal(reducedDimNames(es), c("PCA", "UMAP"),
                           info = name)
    testthat::expect_equal(unname(reducedDim(es, "PCA")), pca * 2,
                           info = name)
    testthat::expect_equal(unname(reducedDim(es, "UMAP")), pca * 10,
                           info = name)
  }
})

# (#26) a plain SCE, SpE and TreeSE with 10 genes and 6 cells
makeReducedDimObjects <- function() {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  list(
    sce = SingleCellExperiment(list(counts = counts)),
    spe = SpatialExperiment(assays = list(counts = counts),
                            spatialCoords = cbind(x = 1:6, y = 1:6)),
    tse = TreeSummarizedExperiment(assays = list(counts = counts),
                                   rowTree = tinyTree)
  )
}

testthat::test_that("reducedDims<- is SingleCellExperiment's generic", {
  # (#26) a package generic of the same name masked SCE's
  testthat::expect_equal(attr(ExperimentSubset::`reducedDims<-`, "package"),
                         "SingleCellExperiment")
})

testthat::test_that("reducedDims<- takes withDimnames", {
  # (#26) the package generic had no withDimnames and no ...
  objects <- makeReducedDimObjects()
  pca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  subsetPca <- pca[2:4, ]
  for (name in names(objects)) {
    plain <- objects[[name]]
    testthat::expect_no_error(
      reducedDims(plain, withDimnames = FALSE) <- list(PCA = pca)
    )
    testthat::expect_equal(unname(reducedDim(plain, "PCA")), pca,
                           info = name)

    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    testthat::expect_no_error(
      reducedDims(es, withDimnames = FALSE) <- list(PCA = pca)
    )
    testthat::expect_equal(unname(reducedDim(es, "PCA")), pca, info = name)
    testthat::expect_no_error(
      reducedDims(es, withDimnames = FALSE, subsetName = "s1") <-
        list(UMAP = subsetPca)
    )
    testthat::expect_equal(reducedDimNames(es, subsetName = "s1"), "UMAP",
                           info = name)
    testthat::expect_equal(
      unname(reducedDim(es, "UMAP", subsetName = "s1")), subsetPca,
      info = name
    )
    # the main object's reducedDims are unchanged
    testthat::expect_equal(reducedDimNames(es), "PCA", info = name)
  }
})

testthat::test_that("reducedDims<- needs subsetName by name", {
  # (#26) the old generic took subsetName second; that is now withDimnames
  objects <- makeReducedDimObjects()
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    testthat::expect_error(
      reducedDims(es, "s1") <- list(PCA = subsetPca),
      "name the subset with subsetName =", info = name
    )
    testthat::expect_equal(reducedDimNames(es, subsetName = "s1"),
                           character(0), info = name)
  }
})

testthat::test_that("reducedDims<- on a plain SCE ignores subsetName", {
  # (#26) as for altExps<-: SingleCellExperiment's method drops it
  sce <- makeReducedDimObjects()$sce
  pca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  reducedDims(sce, subsetName = "s1") <- list(PCA = pca)
  testthat::expect_equal(reducedDimNames(sce), "PCA")
  testthat::expect_equal(unname(reducedDim(sce, "PCA")), pca)
})

testthat::test_that("reducedDimNames() takes only subsetName, by name", {
  # (#82) other arguments, or a subset name by position, were dropped
  objects <- makeReducedDimObjects()
  mainPca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    reducedDim(es, "MAIN") <- mainPca
    reducedDim(es, "PCA", subsetName = "s1") <- subsetPca
    testthat::expect_error(reducedDimNames(es, "s1"),
                           "name the subset with subsetName =", info = name)
    testthat::expect_error(reducedDimNames(es, foo = 1),
                           "takes only 'x' and 'subsetName'", info = name)
    testthat::expect_error(reducedDimNames(es, subsetName = "s1", foo = 1),
                           "takes only 'x' and 'subsetName'", info = name)
    testthat::expect_error(
      reducedDimNames(es, subsetName = "s1", subsetName = "s1"),
      "takes one 'subsetName'", info = name
    )
    # an empty name can't be given as subsetName =
    testthat::expect_error(reducedDimNames(es, character(0)),
                           "takes only 'x' and 'subsetName'", info = name)
    # calls without other arguments are unchanged
    testthat::expect_equal(reducedDimNames(es), "MAIN", info = name)
    testthat::expect_equal(reducedDimNames(es, subsetName = "s1"), "PCA",
                           info = name)
  }
  # a plain SCE has no subsets, so say that rather than how to name one
  sce <- objects$sce
  testthat::expect_error(reducedDimNames(sce, "s1"), "which has no subsets")
  testthat::expect_equal(reducedDimNames(sce), character(0))
})

testthat::test_that("subsetName on a plain SCE says it has no subsets", {
  # (#82) .isSubset() failed in subsetNames() dispatch for these objects
  sce <- makeReducedDimObjects()$sce
  reducedDim(sce, "PCA") <- matrix(as.numeric(seq_len(12)), nrow = 6)
  testthat::expect_error(reducedDimNames(sce, subsetName = "s1"),
                         "which has no subsets")
  testthat::expect_error(reducedDimNames(sce, subsetName = "s1") <- "X",
                         "which has no subsets")
  testthat::expect_error(reducedDimNames(sce, "s1") <- "X",
                         "which has no subsets")
  # without subsetName, the main object is unchanged
  testthat::expect_equal(reducedDimNames(sce), "PCA")
})

testthat::test_that("reducedDims and their names check the subset name", {
  # (#26) each class's methods, not only the last file in Collate
  objects <- makeReducedDimObjects()
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    reducedDims(es, subsetName = "s1") <- list(PCA = subsetPca)
    reducedDimNames(es, subsetName = "s1") <- "renamed"
    testthat::expect_equal(reducedDimNames(es, subsetName = "s1"), "renamed",
                           info = name)
    testthat::expect_equal(reducedDimNames(es), character(0), info = name)
    testthat::expect_error(
      reducedDims(es, subsetName = "nope") <- list(PCA = subsetPca),
      "nope subset does not exist", info = name
    )
    testthat::expect_error(
      reducedDimNames(es, subsetName = "nope"),
      "nope subset does not exist", info = name
    )
    testthat::expect_error(
      reducedDimNames(es, subsetName = "nope") <- "renamed",
      "nope subset does not exist", info = name
    )
  }
})

testthat::test_that("Links and spatialData of a subset keep its names", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  # a row tree only, so the main object has no column links
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  es <- createSubset(es, "s1", rows = 2:5, cols = 1:3)
  links <- subsetRowLinks(es, "s1", parentRowLinkData = TRUE)
  testthat::expect_equal(rownames(links), tinyTree$tip.label[2:5])
  testthat::expect_null(subsetColLinks(es, "s1", parentColLinkData = TRUE))

  spots <- paste0("spot", 1:6)
  spe <- SpatialExperiment(assays = list(counts = counts[, 1:6]),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  colnames(spe) <- spots
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = 2:4)
  # subsetSpatialData() is deprecated (#34)
  testthat::expect_warning(
    s1Data <- subsetSpatialData(es, "s1", parentSpatialData = FALSE),
    "deprecated"
  )
  testthat::expect_equal(rownames(s1Data), spots[2:4])
})

testthat::test_that("withDimnames, withColData and use.names reach SCE", {
  # names set on the object, not on the stored matrix, so the arguments
  # change what comes back
  sce <- SingleCellExperiment(list(counts = matrix(1, 4, 3)))
  rownames(sce) <- paste0("gene", 1:4)
  colnames(sce) <- paste0("cell", 1:3)
  small <- ExperimentSubset(sce)
  testthat::expect_null(rownames(assay(small, "counts", withDimnames = FALSE)))
  # SingleCellExperiment keeps row names here; SummarizedExperiment drops them
  se <- SummarizedExperiment(list(counts = matrix(1, 4, 3)))
  rownames(se) <- paste0("gene", 1:4)
  seRowData <- rowData(ExperimentSubset(se), use.names = FALSE)
  testthat::expect_null(rownames(seRowData))

  es <- ExperimentSubset(sce_chcl)
  counts <- assay(es, "counts")
  # values with other dimnames are accepted when withDimnames = FALSE
  renamed <- counts
  rownames(renamed) <- paste0("other", seq_len(nrow(renamed)))
  testthat::expect_error(assay(es, "renamed", withDimnames = FALSE) <- renamed,
                         NA)
  otherCells <- paste0("other", seq_len(ncol(es)))
  pca <- matrix(1, nrow = ncol(es), ncol = 2,
                dimnames = list(otherCells, NULL))
  testthat::expect_warning(reducedDim(es, "PCA", withDimnames = FALSE) <- pca,
                           NA)
  altCounts <- matrix(1, nrow = 2, ncol = ncol(es),
                      dimnames = list(NULL, otherCells))
  alt <- SingleCellExperiment(list(counts = altCounts))
  testthat::expect_warning(altExp(es, "alt", withDimnames = FALSE) <- alt, NA)
  altList <- list(alt = alt)
  testthat::expect_warning(altExps(es, withDimnames = FALSE) <- altList, NA)

  # and on a subset
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  subsetPca <- pca[1:4, ]
  testthat::expect_warning(reducedDim(es, "PCA", subsetName = "s1",
                                      withDimnames = FALSE) <- subsetPca,
                           NA)
  storedPca <- reducedDim(es, "PCA", subsetName = "s1", withDimnames = FALSE)
  testthat::expect_equal(rownames(storedPca), otherCells[1:4])
})

testthat::test_that("Writing back parent rowData or colData stores no copies", {
  sce <- sce_chcl
  colData(sce)$batch <- rep(c("A", "B"), length.out = ncol(sce))
  rowData(sce)$symbol <- rownames(sce)
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)

  s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE)
  s1ColData$score <- 1:4
  testthat::expect_message(colData(es, subsetName = "s1") <- s1ColData,
                           regexp = "batch")
  testthat::expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  withParent <- colData(es, subsetName = "s1", parentColData = TRUE)
  testthat::expect_equal(sum(colnames(withParent) == "batch"), 1L)

  # a column the user adds with a parent column's name is stored
  own <- colData(es, subsetName = "s1")
  own$batch <- letters[1:4]
  testthat::expect_message(colData(es, subsetName = "s1") <- own, NA)
  testthat::expect_equal(colData(es, subsetName = "s1")$batch, letters[1:4])

  s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE)
  s1RowData$rank <- 5:1
  testthat::expect_message(rowData(es, subsetName = "s1") <- s1RowData,
                           regexp = "symbol")
  testthat::expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")
})

testthat::test_that("A subset's own column hides the main object's (#57)", {
  sce <- simSce
  colData(sce)$score <- 0
  colData(sce)$batch <- rep(c("A", "B"), length.out = ncol(sce))
  rowData(sce)$rank <- 0
  rowData(sce)$symbol <- rownames(sce)
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8,
                     parentAssay = "counts")
  colData(es, subsetName = "s1") <- DataFrame(score = 1:8)
  rowData(es, subsetName = "s1") <- DataFrame(rank = 5:1)

  # different values: one column, the subset's, with a message
  testthat::expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score.*colData\\(x\\)"
  )
  testthat::expect_equal(colnames(s1ColData), c("batch", "score"))
  testthat::expect_equal(s1ColData$score, 1:8)
  testthat::expect_equal(s1ColData$batch, colData(sce)$batch[1:8])
  testthat::expect_message(
    s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE),
    regexp = "rank.*rowData\\(x\\)"
  )
  testthat::expect_equal(colnames(s1RowData), c("symbol", "rank"))
  testthat::expect_equal(s1RowData$rank, 5:1)

  # writing back keeps the subset's values and leaves out the main object's
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  testthat::expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  testthat::expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  testthat::expect_equal(colData(es)$score, rep(0, ncol(sce)))
  suppressMessages(rowData(es, subsetName = "s1") <- s1RowData)
  testthat::expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")
  testthat::expect_equal(rowData(es, subsetName = "s1")$rank, 5:1)
  testthat::expect_equal(rowData(es)$rank, rep(0, nrow(sce)))

  # the extracted object has the same columns, and stores back the same way
  s1 <- suppressMessages(getSubsetExperiment(es, "s1"))
  testthat::expect_equal(colnames(colData(s1)), c("batch", "score"))
  testthat::expect_equal(s1$score, 1:8)
  testthat::expect_equal(colnames(rowData(s1)), c("symbol", "rank"))
  testthat::expect_equal(rowData(s1)$rank, 5:1)
  es <- suppressMessages(setSubsetExperiment(es, "s1", s1))
  testthat::expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  testthat::expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  testthat::expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")

  # a subset of a subset: its own column hides the main object's too
  es <- createSubset(es, "s2", cols = 1:4, parentAssay = "s1")
  colData(es, subsetName = "s2") <- DataFrame(score = 11:14)
  s2ColData <- suppressMessages(
    colData(es, subsetName = "s2", parentColData = TRUE)
  )
  testthat::expect_equal(colnames(s2ColData), c("batch", "score"))
  testthat::expect_equal(s2ColData$score, 11:14)

  # the same values: one column, and no message
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0, 8))
  testthat::expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  testthat::expect_equal(colnames(s1ColData), c("batch", "score"))
  rowData(es, subsetName = "s1") <- DataFrame(rank = rep(0, 5))
  testthat::expect_message(
    s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE),
    NA
  )
  testthat::expect_equal(colnames(s1RowData), c("symbol", "rank"))

  # equal values of another type count as different: integer, not double
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0L, 8))
  testthat::expect_message(
    colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score"
  )
})

testthat::test_that("A main object's repeated column name is checked (#57)", {
  sce <- simSce
  sameScores <- DataFrame(score = rep(0, ncol(sce)))
  otherScores <- DataFrame(score = rep(9, ncol(sce)))
  colData(sce) <- cbind(sameScores, otherScores)
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", cols = 1:8)
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0, 8))
  # the second copy differs from the subset's column
  testthat::expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score"
  )
  testthat::expect_equal(colnames(s1ColData), "score")
})

testthat::test_that("SpatialExperiment subsets have one sample_id (#57)", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = 2:5)
  testthat::expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  testthat::expect_equal(sum(colnames(s1ColData) == "sample_id"), 1L)
  testthat::expect_equal(s1ColData$sample_id, c("A", "A", "B", "B"))
  s1 <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(sum(colnames(colData(s1)) == "sample_id"), 1L)

  # the main object's current sample_id, not the copy createSubset() made
  colData(es)$sample_id <- rep(c("C", "D"), each = 3)
  testthat::expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  testthat::expect_equal(s1ColData$sample_id, c("C", "C", "D", "D"))
  s1 <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(s1$sample_id, c("C", "C", "D", "D"))

  # writing back still works
  s1ColData$group <- letters[1:4]
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  testthat::expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  s1$group <- letters[5:8]
  es <- suppressMessages(setSubsetExperiment(es, "s1", s1))
  testthat::expect_equal(colData(es, subsetName = "s1")$group, letters[5:8])
})

testthat::test_that("Subset data setters take NULL or a data.frame", {
  sce <- sce_chcl
  rowData(sce)$symbol <- rownames(sce)
  colData(sce)$batch <- rep(c("A", "B"), length.out = ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)

  rowData(es, subsetName = "s1") <- data.frame(rank = 5:1)
  testthat::expect_equal(rowData(es, subsetName = "s1")$rank, 5:1)
  rowData(es, subsetName = "s1") <- NULL
  testthat::expect_equal(ncol(rowData(es, subsetName = "s1")), 0L)
  testthat::expect_equal(rowData(es)$symbol, rownames(sce))

  colData(es, subsetName = "s1") <- data.frame(group = letters[1:4])
  testthat::expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  colData(es, subsetName = "s1") <- NULL
  testthat::expect_equal(ncol(colData(es, subsetName = "s1")), 0L)
  testthat::expect_equal(subsetColnames(es, "s1"), colnames(sce)[1:4])
  testthat::expect_equal(colData(es)$batch, colData(sce)$batch)
  # the main object accepts only a DataFrame or NULL, as SummarizedExperiment
  testthat::expect_error(colData(es) <- data.frame(a = seq_len(ncol(es))),
                         regexp = "DataFrame or NULL")

  # subsetRowData<- and subsetColData<- take the same values (#97)
  subsetRowData(es, "s1") <- DataFrame(rank = 5:1)
  fromDataFrame <- subsetRowData(es, "s1")
  subsetRowData(es, "s1") <- data.frame(rank = 5:1)
  testthat::expect_equal(subsetRowData(es, "s1"), fromDataFrame)
  testthat::expect_equal(subsetRowData(es, "s1")$rank, 5:1)
  testthat::expect_equal(subsetRownames(es, "s1"), rownames(sce)[1:5])
  subsetRowData(es, "s1") <- NULL
  testthat::expect_equal(ncol(subsetRowData(es, "s1")), 0L)
  subsetColData(es, "s1") <- DataFrame(group = letters[1:4])
  fromDataFrame <- subsetColData(es, "s1")
  subsetColData(es, "s1") <- data.frame(group = letters[1:4])
  testthat::expect_equal(subsetColData(es, "s1"), fromDataFrame)
  testthat::expect_equal(subsetColData(es, "s1")$group, letters[1:4])
  testthat::expect_equal(subsetColnames(es, "s1"), colnames(sce)[1:4])
  subsetColData(es, "s1") <- NULL
  testthat::expect_equal(ncol(subsetColData(es, "s1")), 0L)
  valueMsg <- "'value' must be a DataFrame, a data.frame or NULL."
  testthat::expect_error(subsetRowData(es, "s1") <- 1:5, valueMsg,
                         fixed = TRUE)
  testthat::expect_error(subsetColData(es, "s1") <- 1:4, valueMsg,
                         fixed = TRUE)
  testthat::expect_error(rowData(es, subsetName = "s1") <- 1:5, valueMsg,
                         fixed = TRUE)

  # the parent-column mark is never stored
  s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE)
  rowData(es) <- rowData(es)
  testthat::expect_null(S4Vectors::mcols(rowData(es))$fromParentData)
  s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE)
  s1ColData$score <- 1:4
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  storedInfo <- S4Vectors::mcols(colData(es, subsetName = "s1"))
  testthat::expect_null(storedInfo$fromParentData)
})

testthat::test_that("colData<- rejects other values in every Subset class", {
  # (#99) other values failed at S4 dispatch, before the package's check
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  se <- SummarizedExperiment(list(counts = counts))
  objects <- list(
    se = se,
    rse = as(se, "RangedSummarizedExperiment"),
    sce = SingleCellExperiment(list(counts = counts)),
    tse = TreeSummarizedExperiment(assays = list(counts = counts),
                                   rowTree = tinyTree),
    spe = SpatialExperiment(assays = list(counts = counts),
                            spatialCoords = cbind(x = 1:6, y = 1:6))
  )
  subsetMsg <- "'value' must be a DataFrame, a data.frame or NULL."
  mainMsg <- "'value' must be a DataFrame or NULL."
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
    mainColData <- colData(es)
    testthat::expect_error(colData(es, subsetName = "s1") <- 1:4,
                           subsetMsg, fixed = TRUE, info = name)
    testthat::expect_error(colData(es, subsetName = "s1") <- "a",
                           subsetMsg, fixed = TRUE, info = name)
    testthat::expect_error(colData(es) <- seq_len(ncol(es)),
                           mainMsg, fixed = TRUE, info = name)
    # a DataFrame or NULL with subsetName still sets the subset's colData,
    # not the main object's
    colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4])
    testthat::expect_equal(colData(es, subsetName = "s1")$group,
                           letters[1:4], info = name)
    testthat::expect_equal(colData(es), mainColData, info = name)
    colData(es, subsetName = "s1") <- NULL
    # a SpatialExperiment subset keeps its sample_id column
    keptColumns <- if (name == "spe") "sample_id" else character(0)
    testthat::expect_equal(colnames(colData(es, subsetName = "s1")),
                           keptColumns, info = name)
    testthat::expect_equal(colData(es), mainColData, info = name)
    # a data.frame with subsetName sets the subset's colData too
    colData(es, subsetName = "s1") <- data.frame(group = letters[5:8])
    testthat::expect_equal(colData(es, subsetName = "s1")$group,
                           letters[5:8], info = name)
    testthat::expect_equal(colData(es), mainColData, info = name)
  }
})

testthat::test_that("altExps and reducedDims of a subset pass arguments on", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  cells <- subsetColnames(es, "s1")
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4],
                                              row.names = cells)
  altCounts <- matrix(1, nrow = 2, ncol = 4, dimnames = list(NULL, cells))
  altExp(es, "alt", subsetName = "s1") <-
    SingleCellExperiment(list(counts = altCounts))
  withGroup <- altExps(es, subsetName = "s1", withColData = TRUE)
  testthat::expect_equal(colData(withGroup$alt)$group, letters[1:4])
  withoutGroup <- altExps(es, subsetName = "s1")
  testthat::expect_null(colData(withoutGroup$alt)$group)

  pca <- matrix(1, nrow = 4, ncol = 2,
                dimnames = list(paste0("other", 1:4), NULL))
  reducedDim(es, "PCA", subsetName = "s1", withDimnames = FALSE) <- pca
  storedPcas <- reducedDims(es, subsetName = "s1", withDimnames = FALSE)
  testthat::expect_equal(rownames(storedPcas$PCA), paste0("other", 1:4))
})

testthat::test_that("Subset getters use renamed and repeated names", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  es <- createSubset(es, "s1", rows = 2:5, cols = 1:3)
  subsetRownames(es, "s1") <- paste0("g", 1:4)
  links <- subsetRowLinks(es, "s1", parentRowLinkData = TRUE)
  testthat::expect_equal(rownames(links), paste0("g", 1:4))

  # spot barcodes repeat across the two samples
  barcodes <- rep(paste0("b", 1:3), 2)
  spe <- SpatialExperiment(assays = list(counts = counts[, 1:6]),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  colnames(spe) <- barcodes
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = c(1, 2, 4, 5))
  # subsetSpatialData() is deprecated (#34)
  testthat::expect_warning(s1Data <- subsetSpatialData(es, "s1", FALSE),
                           "deprecated")
  testthat::expect_equal(rownames(s1Data), c("b1", "b2", "b1", "b2"))
  subsetColnames(es, "s1") <- paste0("spot", 1:4)
  testthat::expect_warning(s1Data <- subsetSpatialData(es, "s1", TRUE),
                           "deprecated")
  testthat::expect_equal(rownames(s1Data), paste0("spot", 1:4))
})

testthat::test_that("subsetParent lists each parent up to the main assay", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:10, cols = 1:10)
  es <- createSubset(es, "s2", rows = 1:5, cols = 1:5, parentAssay = "s1")
  testthat::expect_equal(subsetParent(es, "s2"), list("s1", "counts"))
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "stored")
  testthat::expect_equal(subsetParent(es, "s1/stored"), list("s1", "counts"))
  es <- createSubset(es, "s3", rows = 1:2, cols = 1:2,
                     parentAssay = "s1/stored")
  testthat::expect_equal(subsetParent(es, "s3"),
                         list("s1/stored", "s1", "counts"))
})

testthat::test_that("Subsets see their parent's assays as linked views", {
  es <- ExperimentSubset(sce_chcl)
  assay(es, "logcounts") <- log1p(assay(es, "counts"))
  es <- createSubset(es, "filtered", rows = 1:20, cols = 1:10)
  es <- createSubset(es, "tcells", rows = 1:5, cols = 2:4,
                     parentAssay = "filtered")
  mainLog <- as.matrix(assay(es, "logcounts"))

  # a main assay is linked into every subset under the same name
  filteredLog <- as.matrix(assay(es, "filtered/logcounts"))
  testthat::expect_equal(filteredLog, mainLog[1:20, 1:10])
  tcellLog <- as.matrix(assay(es, "tcells", subsetAssayName = "logcounts"))
  testthat::expect_equal(tcellLog, mainLog[1:5, 2:4])
  testthat::expect_equal(subsetAssayNames(es, "tcells"),
                         c("counts", "logcounts"))

  # an assay stored in filtered is linked into tcells
  filteredScaled <- filteredLog * 2
  es <- setSubsetAssay(es, "filtered", filteredScaled, "scaled")
  tcellView <- as.matrix(assay(es, "tcells/scaled"))
  testthat::expect_equal(tcellView, filteredScaled[1:5, 2:4])

  # storing over a linked name needs unlink = TRUE
  tcellScaled <- tcellLog + 1
  testthat::expect_error(setSubsetAssay(es, "tcells", tcellScaled, "scaled"),
                         regexp = "linked to 'filtered/scaled'")
  es <- setSubsetAssay(es, "tcells", tcellScaled, "scaled", unlink = TRUE)
  testthat::expect_equal(as.matrix(assay(es, "tcells/scaled")), tcellScaled)
  # the parent's assay is unchanged
  testthat::expect_equal(as.matrix(assay(es, "filtered/scaled")),
                         filteredScaled)
  # a subset of tcells sees tcells' own scaled
  es <- createSubset(es, "cd4", rows = 1:2, cols = 1:2,
                     parentAssay = "tcells")
  testthat::expect_equal(as.matrix(assay(es, "cd4/scaled")),
                         tcellScaled[1:2, 1:2])
  # main assays are linked too
  filteredCounts <- as.matrix(assay(es, "filtered"))
  counts <- filteredCounts
  testthat::expect_error(setSubsetAssay(es, "filtered", counts, "logcounts"),
                         regexp = "linked to the main assay 'logcounts'")
})

testthat::test_that("Sibling subsets can store the same assay name", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  es <- createSubset(es, "s2", rows = 6:10, cols = 5:8)
  s1Norm <- as.matrix(assay(es, "s1")) + 1
  s2Norm <- as.matrix(assay(es, "s2")) + 2
  es <- setSubsetAssay(es, "s1", s1Norm, "norm")
  es <- setSubsetAssay(es, "s2", s2Norm, "norm")

  # a qualified name or subsetAssayName picks one
  testthat::expect_equal(as.matrix(assay(es, "s1/norm")), s1Norm)
  s2Stored <- assay(es, "s2", subsetAssayName = "norm")
  testthat::expect_equal(as.matrix(s2Stored), s2Norm)
  s1Stored <- getSubsetAssay(es, "s1", subsetAssayName = "norm")
  testthat::expect_equal(as.matrix(s1Stored), s1Norm)
  testthat::expect_equal(as.matrix(getSubsetAssay(es, "s2/norm")), s2Norm)
  # a bare name means a main assay or a subset
  testthat::expect_error(assay(es, "norm"),
                         regexp = "stored in subset 's1', 's2'")
  es <- setSubsetAssay(es, "s1", s1Norm, "onlyInS1")
  testthat::expect_error(assay(es, "onlyInS1"), regexp = "s1/onlyInS1")
  testthat::expect_error(createSubset(es, "s3", 1:2, parentAssay = "norm"),
                         regexp = "stored in subset")

  # a qualified parentAssay makes a subset of that stored assay
  es <- createSubset(es, "s3", rows = 1:2, cols = 1:2,
                     parentAssay = "s2/norm")
  testthat::expect_equal(subsetRownames(es, "s3"), rownames(sce_chcl)[6:7])
  testthat::expect_equal(as.matrix(assay(es, "s3")), s2Norm[1:2, 1:2])
  testthat::expect_equal(subsetParent(es, "s3"),
                         list("s2/norm", "s2", "counts"))
  testthat::expect_true("s2/norm" %in% subsetAssayNames(es))

  # a later sibling storing "norm" doesn't change what s3 points to
  es <- createSubset(es, "s4", rows = 11:15, cols = 9:12)
  es <- setSubsetAssay(es, "s4", as.matrix(assay(es, "s4")), "norm")
  testthat::expect_equal(as.matrix(assay(es, "s3")), s2Norm[1:2, 1:2])

  # assay<- with a qualified name stores into that subset
  assay(es, "s1/extra") <- s1Norm
  testthat::expect_equal(as.matrix(assay(es, "s1", subsetAssayName = "extra")),
                         s1Norm)
})

testthat::test_that("A subset's own assay is used by it and its descendants", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "P", rows = 1:5, cols = 1:4)
  es <- createSubset(es, "Q", rows = 1:5, cols = 1:4)
  mainCounts <- as.matrix(assay(es, "counts"))
  pCounts <- as.matrix(assay(es, "P")) + 100
  es <- setSubsetAssay(es, "P", pCounts, "counts", unlink = TRUE)
  # P's primary assay is now its own counts
  testthat::expect_equal(as.matrix(assay(es, "P")), pCounts)
  es <- createSubset(es, "C", rows = 1:2, cols = 1:2, parentAssay = "P")
  testthat::expect_equal(as.matrix(assay(es, "C")), pCounts[1:2, 1:2])
  testthat::expect_equal(as.matrix(assay(es, "C/counts")), pCounts[1:2, 1:2])
  # each name is listed once, own or linked
  testthat::expect_equal(subsetAssayNames(es, "C"), assayNames(es))
  # the sibling Q still sees the main counts
  testthat::expect_equal(as.matrix(assay(es, "Q")), mainCounts[1:5, 1:4])
})

testthat::test_that("assay<- stores linked names only with unlink = TRUE", {
  es <- ExperimentSubset(sce_chcl)
  assay(es, "logcounts") <- log1p(assay(es, "counts"))
  es <- createSubset(es, "filtered", rows = 1:20, cols = 1:10)
  filteredLog <- as.matrix(assay(es, "filtered/logcounts"))
  testthat::expect_error(assay(es, "filtered/logcounts") <- filteredLog,
                         regexp = "unlink = TRUE")
  assay(es, "filtered/logcounts", unlink = TRUE) <- filteredLog + 1
  testthat::expect_equal(as.matrix(assay(es, "filtered/logcounts")),
                         filteredLog + 1)
  plusTwo <- filteredLog + 2
  suppressMessages(assay(es, "filtered", subsetAssayName = "logcounts",
                         unlink = TRUE) <- plusTwo)
  testthat::expect_equal(as.matrix(assay(es, "filtered/logcounts")),
                         filteredLog + 2)
  testthat::expect_output(show(es), "filtered/logcounts")

  # a subset can be made from an assay its parent only links to
  es <- createSubset(es, "fromLinked", rows = 1:3, cols = 1:3,
                     parentAssay = "filtered/counts")
  mainCounts <- as.matrix(assay(es, "counts"))
  testthat::expect_equal(as.matrix(assay(es, "fromLinked")),
                         mainCounts[1:3, 1:3])
  linkedView <- assay(es, "filtered/counts", withDimnames = FALSE)
  testthat::expect_equal(dim(linkedView), c(20L, 10L))

  # arguments that don't fit the call are errors
  testthat::expect_error(assay(es, "filtered/x", subsetAssayName = "y") <-
                           filteredLog, regexp = "already names")
  testthat::expect_error(assay(es, "counts", unlink = TRUE) <- mainCounts,
                         regexp = "only applies")
  testthat::expect_error(setSubsetAssay(es, "filtered", filteredLog, ""),
                         regexp = "non-empty")
})

testthat::test_that("Saved objects whose subsets record bare parents work", {
  # build a subset the way older versions recorded it, with a bare parent
  .subsets <- ExperimentSubset:::.subsets
  `.subsets<-` <- ExperimentSubset:::`.subsets<-`
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  stored <- as.matrix(assay(es, "s1")) + 1
  es <- setSubsetAssay(es, "s1", stored, "stored")
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2,
                     parentAssay = "s1/stored")
  child <- .subsets(es)[["child"]]
  makeAssaySubset <- ExperimentSubset:::AssaySubset
  childRows <- ExperimentSubset:::.rowIndices(child)
  childCols <- ExperimentSubset:::.colIndices(child)
  childData <- ExperimentSubset:::.internalAssay(child)
  oldChild <- makeAssaySubset(subsetName = "child", rowIndices = childRows,
                              colIndices = childCols, parentAssay = "stored",
                              internalAssay = childData)
  subsetList <- .subsets(es)
  subsetList[["child"]] <- oldChild
  .subsets(es) <- subsetList

  testthat::expect_equal(as.matrix(assay(es, "child")), stored[1:2, 1:2])
  testthat::expect_equal(subsetParent(es, "child"),
                         list("stored", "s1", "counts"))
  # a new subset called "stored" doesn't change what the child points to
  es <- createSubset(es, "stored", rows = 10:20, cols = 10:20)
  testthat::expect_equal(as.matrix(assay(es, "child")), stored[1:2, 1:2])
  testthat::expect_equal(subsetParent(es, "child"),
                         list("s1/stored", "s1", "counts"))
})

testthat::test_that("getSubsetAssay reads stored assays in SE and RSE", {
  counts <- matrix(as.numeric(1:20), nrow = 4,
                   dimnames = list(paste0("g", 1:4), paste0("c", 1:5)))
  se <- SummarizedExperiment(list(counts = counts))
  rse <- as(se, "RangedSummarizedExperiment")
  for (object in list(se, rse)) {
    es <- ExperimentSubset(object)
    testthat::expect_true(is(es, paste0("Subset", class(object))))
    es <- createSubset(es, "s1", rows = 1:2, cols = 1:3)
    doubled <- counts[1:2, 1:3] * 2
    es <- setSubsetAssay(es, "s1", doubled, "doubled")
    s1Doubled <- getSubsetAssay(es, "s1", subsetAssayName = "doubled")
    testthat::expect_equal(as.matrix(s1Doubled), doubled)
  }
})

testthat::test_that("Qualified names that don't exist give clear errors", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  s1Counts <- as.matrix(assay(es, "s1"))
  testthat::expect_error(assay(es, "s1/none"), regexp = "s1/none")
  testthat::expect_error(assay(es, "s1/a/b"), regexp = "s1/a/b")
  testthat::expect_error(assay(es, "nope/counts"), regexp = "nope/counts")
  testthat::expect_error(assay(es, "s1", subsetAssayName = c("a", "b")),
                         regexp = "subsetAssayName")
  # subsetAssayName with an i that isn't a subset name gives the same
  # reason as the setter, not that the subset doesn't exist (#76)
  testthat::expect_error(assay(es, "s1/counts", subsetAssayName = "x"),
                         regexp = "already names")
  testthat::expect_error(assay(es, "counts", subsetAssayName = "x"),
                         regexp = "not a subset")
  # a "/" name that isn't "subset/assay" for an existing subset, as for
  # the setter
  testthat::expect_error(assay(es, "nope/counts", subsetAssayName = "x"),
                         regexp = "existing subset")
  testthat::expect_error(assay(es, "s1/a/b", subsetAssayName = "x"),
                         regexp = "existing subset")
  # getSubsetAssay() names its own argument, subsetName, not assay()'s i
  testthat::expect_error(getSubsetAssay(es, "counts", subsetAssayName = "x"),
                         regexp = "'subsetName' must be")
  testthat::expect_error(assay(es, "nope/extra") <- s1Counts,
                         regexp = "existing subset")
  testthat::expect_error(subsetParent(es, "counts"),
                         regexp = "does not exist")
})

testthat::test_that("Subset and stored assay names can't contain a slash", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  s1Counts <- as.matrix(assay(es, "s1"))
  testthat::expect_error(createSubset(es, "a/b", rows = 1:2), regexp = "/")
  testthat::expect_error(setSubsetAssay(es, "s1", s1Counts, "a/b"),
                         regexp = "/")
})

testthat::test_that("Parent lookups follow .assay() for ambiguous names", {
  sce <- sce_chcl
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 101:200, cols = 501:600)
  es <- createSubset(es, "t1", rows = 1:10, cols = 1:10)
  # two subsets store an assay named "dup", so a plain "dup" is ambiguous
  # and the parent is named with a qualified name (#2)
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "dup")
  es <- setSubsetAssay(es, "t1", assay(es, "t1"), "dup")
  testthat::expect_error(createSubset(es, "s5", 2:3, 2:3, parentAssay = "dup"),
                         regexp = "stored in subset")
  es <- createSubset(es, "s5", rows = 2:3, cols = 2:3, parentAssay = "t1/dup")
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
  testthat::expect_equal(dimnames(assay(es, "s1/unnamed")), list(genes, cells))

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
  roundTrip <- assay(es, "s1/roundTrip")
  testthat::expect_equal(colnames(roundTrip), paste0("cell", 1:4))
})

testthat::test_that("A subset of a subset round-trips through setSubsetAssay", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "p", rows = 5:8, cols = 3:6)
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2, parentAssay = "p")
  childCounts <- assay(es, "child")
  testthat::expect_equal(rownames(childCounts), rownames(sce_chcl)[5:6])
  testthat::expect_equal(colnames(childCounts), colnames(sce_chcl)[3:4])
  es <- setSubsetAssay(es, "child", childCounts, "copy")
  childCopy <- assay(es, "child/copy")
  testthat::expect_equal(dimnames(childCopy), dimnames(childCounts))
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
  testthat::expect_equal(subsetParent(es, "mysubset/stored"),
                         list("mysubset", "counts"))
  # the cleaned name can't overwrite an existing subset
  testthat::expect_error(suppressWarnings(createSubset(es, spaced, 1:2, 1:2)),
                         regexp = "'mysubset' is already a subset name.",
                         fixed = TRUE)
})

testthat::test_that("createSubset rejects invalid rows and cols", {
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_error(createSubset(es, "s1", rows = -1),
                         regexp = "positive whole numbers")
  testthat::expect_error(createSubset(es, "s1", rows = c(0, 1, 2)),
                         regexp = "positive whole numbers")
  testthat::expect_error(createSubset(es, "s1", cols = 1.5),
                         regexp = "positive whole numbers")
  testthat::expect_error(createSubset(es, "s1", rows = c(1, 1)),
                         regexp = "duplicates")
  testthat::expect_error(createSubset(es, "s1", cols = c(2, 3, 2)),
                         regexp = "duplicates")
  sameGene <- rownames(sce_chcl)[c(1, 1)]
  testthat::expect_error(createSubset(es, "s1", rows = sameGene),
                         regexp = "duplicates")
  testthat::expect_error(createSubset(es, "s1", rows = integer(0)),
                         regexp = "at least one")

  # names given for an assay without row or column names
  unnamedSce <- SingleCellExperiment(list(counts = matrix(1, 6, 5)))
  unnamed <- ExperimentSubset(unnamedSce)
  testthat::expect_error(createSubset(unnamed, "s1", rows = c("a", "b")),
                         regexp = "no row names")
  testthat::expect_error(createSubset(unnamed, "s1", cols = "a"),
                         regexp = "no column names")
})

testthat::test_that("createSubset rejects invalid names", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "stored")

  testthat::expect_error(createSubset(es, c("a", "b"), rows = 1:2),
                         regexp = "subsetName")
  testthat::expect_error(createSubset(es, NA_character_, rows = 1:2),
                         regexp = "subsetName")
  testthat::expect_error(createSubset(es, "", rows = 1:2),
                         regexp = "subsetName")
  # spaces are removed first, so a name of only spaces is empty
  testthat::expect_error(suppressWarnings(createSubset(es, "  ", rows = 1:2)),
                         regexp = "subsetName")
  # a subset can't share its name with an assay
  testthat::expect_error(createSubset(es, "counts", rows = 1:2),
                         regexp = "already an assay name")
  # a stored assay is named "s1/stored", so a subset may be called "stored"
  es <- createSubset(es, "stored", rows = 1:2)
  testthat::expect_true("stored" %in% subsetNames(es))
  twoParents <- c("counts", "s1")
  testthat::expect_error(createSubset(es, "s2", 1:2, parentAssay = twoParents),
                         regexp = "parentAssay")
})

testthat::test_that("ExperimentSubset checks the subset argument and lists", {
  withNA <- list(subsetName = "s1", rows = 1:10, cols = NA, parentAssay = NA)
  testthat::expect_error(ExperimentSubset(sce_chcl, subset = withNA),
                         regexp = "'subset'")
  testthat::expect_error(ExperimentSubset(list(matrix(1, 3, 3))),
                         regexp = "unnamed assays")

  # leaving out cols and parentAssay keeps all columns and the first assay
  es <- ExperimentSubset(sce_chcl, subset = list(subsetName = "s1",
                                                 rows = 1:10))
  testthat::expect_equal(subsetDim(es, "s1"), c(10L, ncol(sce_chcl)))
})

testthat::test_that("Unknown subset names give a clear error", {
  es <- ExperimentSubset(sce_chcl)
  # no subsets yet
  testthat::expect_error(assay(es, "nope"), regexp = "nope")
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  testthat::expect_error(assay(es, "nope"), regexp = "nope")
  testthat::expect_error(getSubsetAssay(es, "nope"), regexp = "nope")
  testthat::expect_error(subsetDim(es, "nope"),
                         regexp = "nope subset does not exist")
  testthat::expect_error(
    subsetParent(es, "nope"),
    regexp = "nope does not exist in the subsets slot of the object."
  )
  testthat::expect_error(subsetColnames(es, "nope"),
                         regexp = "nope subset does not exist")
  testthat::expect_error(subsetRownames(es, "nope"),
                         regexp = "nope subset does not exist")
  testthat::expect_error(subsetColnames(es, "nope") <- letters[1:5],
                         regexp = "nope subset does not exist")
  testthat::expect_error(subsetRownames(es, "nope") <- letters[1:5],
                         regexp = "nope subset does not exist")
  testthat::expect_error(reducedDims(es, subsetName = "nope") <- list(),
                         regexp = "nope subset does not exist")
})

testthat::test_that("assay<- stores into a subset chosen by its name", {
  es <- ExperimentSubset(sce_chcl)
  # a subset with the same dimensions as the main assay
  es <- createSubset(es, "all")
  counts <- assay(es, "counts")
  suppressMessages(assay(es, "all", subsetAssayName = "copy") <- counts)
  testthat::expect_true("all/copy" %in% subsetAssayNames(es))
  testthat::expect_false("all" %in% assayNames(es))

  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  testthat::expect_error(assay(es, "s1") <- assay(es, "s1"),
                         regexp = "is needed to store")
  # subsetAssayName with a name that is not a subset
  s1Counts <- assay(es, "s1")
  testthat::expect_error(assay(es, "s1typo", subsetAssayName = "n") <- s1Counts,
                         regexp = "not a subset")
})

testthat::test_that("Testing createSubset with multiple paramter options",{
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_error(es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c(20,21,40,45,90,99,100,123,166,299), parentAssay = "none"), "Input parentAssay does not exist.")
  
  #colnames
  es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
  namedCols <- c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC")
  testthat::expect_equal(subsetNames(es), "subset1")
  expectSubsetOf(es, "subset1",
                 rows = c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200),
                 cols = match(namedCols, colnames(es)))
  
  #null rows and null cols
  es <- createSubset(es, "subset2", parentAssay = "counts")
  testthat::expect_equal(subsetNames(es), c("subset1", "subset2"))
  expectSubsetOf(es, "subset2", rows = seq_len(nrow(es)),
                 cols = seq_len(ncol(es)))
  
  # rows or columns given by names that the parent doesn't have
  testthat::expect_error(
    createSubset(es, "subset22", rows = c("as", "asd", "asd"),
                 cols = namedCols, parentAssay = "counts"),
    regexp = "Selected rows not available in the specified assay."
  )
  testthat::expect_error(
    createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"),
                 parentAssay = "counts"),
    regexp = "Selected columns not available in the specified assay."
  )
})

testthat::test_that("Completing test coverage, move this into relevant places later",{
  es <- ExperimentSubset(sce_chcl)
  testthat::expect_warning(es <- createSubset(es, "subset 1", rows = c(1:5), cols = (1:5)),
                           "Removing spaces from the specified subsetName.")
  testthat::expect_equal(subsetNames(es), "subset1")
  
  # a named list of assays gives a SubsetSingleCellExperiment of them
  es <- ExperimentSubset(list(counts = assay(sce_chcl, "counts")))
  testthat::expect_s4_class(es, "SubsetSingleCellExperiment")
  testthat::expect_equal(assayNames(es), "counts")
  testthat::expect_equal(assay(es, "counts"), assay(sce_chcl, "counts"))
  testthat::expect_equal(subsetCount(es), 0)
})
# Takes a subset out as a standalone experiment and puts it back unchanged
checkExperimentRoundTrip <- function(object) {
  object$batch <- rep(c("a", "b"), 3)
  es <- ExperimentSubset(object)
  es <- createSubset(es, "s1", rows = 2:5, cols = 2:5)
  sub <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(class(sub), class(object))
  testthat::expect_equal(dim(sub), c(4L, 4L))
  testthat::expect_equal(rownames(sub), subsetRownames(es, "s1"))
  testthat::expect_equal(colnames(sub), subsetColnames(es, "s1"))
  testthat::expect_equal(as.matrix(assay(sub, "counts")),
                         as.matrix(assay(es, "s1")))
  testthat::expect_equal(sub$batch, object$batch[2:5])
  # unchanged linked assays and parent columns are not stored
  ownColumns <- colnames(subsetColData(es, "s1", FALSE))
  es <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(subsetAssayNames(es), "s1")
  testthat::expect_equal(colnames(subsetColData(es, "s1", FALSE)),
                         ownColumns)
  invisible(sub)
}

testthat::test_that("Subset experiments round-trip for every class", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("cell", 1:6)))
  se <- SummarizedExperiment(list(counts = counts))
  checkExperimentRoundTrip(se)
  checkExperimentRoundTrip(as(se, "RangedSummarizedExperiment"))
  checkExperimentRoundTrip(SingleCellExperiment(list(counts = counts)))

  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  treeCounts <- counts
  rownames(treeCounts) <- tinyTree$tip.label
  tse <- TreeSummarizedExperiment(assays = list(counts = treeCounts),
                                  rowTree = tinyTree)
  treeSub <- checkExperimentRoundTrip(tse)
  # the tree and links come from the main object
  testthat::expect_false(is.null(rowTree(treeSub)))
  testthat::expect_equal(rowLinks(treeSub)$nodeLab, rowLinks(tse)$nodeLab[2:5])

  coords <- cbind(x = as.numeric(1:6), y = as.numeric(1:6) * 10)
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = coords)
  spatialSub <- checkExperimentRoundTrip(spe)
  testthat::expect_equal(unname(spatialCoords(spatialSub)),
                         unname(coords[2:5, ]))
  testthat::expect_equal(spatialSub$sample_id, spe$sample_id[2:5])
})

testthat::test_that("setSubsetExperiment stores what a tool added", {
  sce <- simSce
  sce$batch <- rep(c("a", "b"), 60)
  rowData(sce)$symbol <- toupper(rownames(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  # parent columns are included and tagged
  testthat::expect_equal(sub$batch, sce$batch[1:10])
  testthat::expect_equal(rowData(sub)$symbol, rowData(sce)$symbol[1:20])

  scaled <- as.matrix(assay(sub, "counts")) / 2
  assay(sub, "scaled") <- scaled
  pca <- matrix(as.numeric(1:20), nrow = 10,
                dimnames = list(colnames(sub), c("PC1", "PC2")))
  reducedDim(sub, "PCA") <- pca
  sub$cluster <- rep(1:2, 5)
  rowData(sub)$variable <- rep(c(TRUE, FALSE), 10)
  metadata(sub)$note <- "clustered"
  es <- setSubsetExperiment(es, "s1", sub)

  testthat::expect_equal(subsetAssayNames(es), c("s1", "s1/scaled"))
  testthat::expect_equal(as.matrix(assay(es, "s1/scaled")), scaled)
  testthat::expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca)
  testthat::expect_equal(colnames(colData(es, subsetName = "s1")), "cluster")
  testthat::expect_equal(colData(es, subsetName = "s1")$cluster, rep(1:2, 5))
  testthat::expect_equal(colnames(rowData(es, subsetName = "s1")),
                         "variable")
  testthat::expect_equal(metadata(es, subsetName = "s1")$note, "clustered")
  # the main object is unchanged
  testthat::expect_equal(colnames(colData(es)), "batch")
  testthat::expect_equal(assayNames(es), "counts")

  # without parent data only the subset's own columns come back
  own <- getSubsetExperiment(es, "s1", parentData = FALSE)
  testthat::expect_equal(colnames(colData(own)), "cluster")
  testthat::expect_equal(colnames(rowData(own)), "variable")
  testthat::expect_equal(assayNames(own), c("counts", "scaled"))
})

testthat::test_that("setSubsetExperiment replaces the subset's own data", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  assay(sub, "scaled") <- as.matrix(assay(sub, "counts")) / 2
  reducedDim(sub, "PCA") <- matrix(0, nrow = 10, ncol = 2)
  es <- setSubsetExperiment(es, "s1", sub)
  # what is removed from value is removed from the subset
  sub <- getSubsetExperiment(es, "s1")
  assay(sub, "scaled") <- NULL
  reducedDim(sub, "PCA") <- NULL
  es <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(subsetAssayNames(es), "s1")
  testthat::expect_equal(reducedDimNames(es, subsetName = "s1"), character(0))
})

testthat::test_that("Changed linked assays are stored only with unlink", {
  sce <- simSce
  logcounts(sce) <- log1p(counts(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(assayNames(sub), c("counts", "logcounts"))
  # unchanged, so not copied
  es <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(subsetAssayNames(es), "s1")

  changed <- as.matrix(logcounts(sub)) * 2
  logcounts(sub) <- changed
  testthat::expect_error(setSubsetExperiment(es, "s1", sub),
                         regexp = "logcounts.*unlink = TRUE")
  es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE)
  testthat::expect_equal(subsetAssayNames(es), c("s1", "s1/logcounts"))
  testthat::expect_equal(as.matrix(assay(es, "s1/logcounts")), changed)
  # the main object's logcounts is unchanged
  testthat::expect_equal(assay(es, "logcounts"), logcounts(sce))
})

testthat::test_that("setSubsetExperiment checks value against the subset", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  testthat::expect_error(setSubsetExperiment(es, "s1", sub[1:5, ]),
                         regexp = "dimensions")
  testthat::expect_error(setSubsetExperiment(es, "s1", sub[, 10:1]),
                         regexp = "colnames")
  testthat::expect_error(setSubsetExperiment(es, "s1", sub[20:1, ]),
                         regexp = "rownames")
  testthat::expect_error(setSubsetExperiment(es, "none", sub),
                         regexp = "none")
  testthat::expect_error(getSubsetExperiment(es, "none"), regexp = "none")
  se <- SummarizedExperiment(list(counts = counts(sub)))
  testthat::expect_error(setSubsetExperiment(es, "s1", se),
                         regexp = "SingleCellExperiment")
  # an SCE for a SummarizedExperiment main object
  seEs <- ExperimentSubset(SummarizedExperiment(list(counts = simCounts)))
  seEs <- createSubset(seEs, "s1", rows = 1:20, cols = 1:10)
  testthat::expect_error(setSubsetExperiment(seEs, "s1", sub),
                         regexp = "must be a SummarizedExperiment")
  bad <- sub
  assay(bad, "a/b") <- counts(sub)
  testthat::expect_error(setSubsetExperiment(es, "s1", bad), regexp = "/")
})

testthat::test_that("Subset experiments of subsets use the subset's assays", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  scaled <- as.matrix(assay(es, "s1")) / 2
  es <- setSubsetAssay(es, "s1", scaled, "scaled")
  es <- createSubset(es, "s2", rows = 1:5, cols = 1:4,
                     parentAssay = "s1/scaled")
  sub <- getSubsetExperiment(es, "s2")
  # the subset's primary assay comes first
  testthat::expect_equal(assayNames(sub), c("scaled", "counts"))
  testthat::expect_equal(as.matrix(assay(sub, "scaled")), scaled[1:5, 1:4])
  testthat::expect_equal(as.matrix(assay(sub)), as.matrix(assay(es, "s2")))
  es <- setSubsetExperiment(es, "s2", sub)
  testthat::expect_equal(subsetAssayNames(es), c("s1", "s2", "s1/scaled"))
  # s1's scaled can't be dropped while s2 is made from it
  s1 <- getSubsetExperiment(es, "s1")
  assay(s1, "scaled") <- NULL
  testthat::expect_error(setSubsetExperiment(es, "s1", s1),
                         regexp = "'scaled', which subsets 's2'")
  testthat::expect_equal(as.matrix(assay(es, "s2")), scaled[1:5, 1:4])
})

testthat::test_that("setSubsetExperiment keeps assays that descendants use", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  scaled <- as.matrix(assay(es, "s1")) / 2
  es <- setSubsetAssay(es, "s1", scaled, "scaled")
  # child reaches s1's scaled through another subset
  es <- createSubset(es, "child", rows = 1:10, cols = 1:8, parentAssay = "s1")
  es <- createSubset(es, "grandchild", rows = 1:5, cols = 1:4,
                     parentAssay = "child/scaled")
  s1 <- getSubsetExperiment(es, "s1")
  assay(s1, "scaled") <- NULL
  testthat::expect_error(setSubsetExperiment(es, "s1", s1),
                         regexp = "'scaled', which subsets 'grandchild'")

  # a subset made from s1 uses s1's own counts, stored with unlink = TRUE
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  ownCounts <- as.matrix(assay(es, "s1")) * 100
  es <- setSubsetAssay(es, "s1", ownCounts, "counts", unlink = TRUE)
  es <- createSubset(es, "child", rows = 1:5, cols = 1:4, parentAssay = "s1")
  s1 <- getSubsetExperiment(es, "s1")
  assay(s1, "counts") <- NULL
  testthat::expect_error(setSubsetExperiment(es, "s1", s1),
                         regexp = "'counts', which subsets 'child'")
  testthat::expect_equal(as.matrix(assay(es, "child")), ownCounts[1:5, 1:4])
})

testthat::test_that("Subset experiments keep names a subset was given", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  es <- createSubset(es, "s1", rows = 1:5, cols = 2:5)
  subsetRownames(es, "s1") <- paste0("row", 1:5)
  subsetColnames(es, "s1") <- paste0("col", 1:4)
  sub <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(rownames(sub), paste0("row", 1:5))
  testthat::expect_equal(colnames(sub), paste0("col", 1:4))
  testthat::expect_equal(rowLinks(sub)$nodeLab, tinyTree$tip.label[1:5])
  es <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(subsetRownames(es, "s1"), paste0("row", 1:5))
})

testthat::test_that("setSubsetExperiment needs the subset's names on value", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  noColnames <- sub
  colnames(noColnames) <- NULL
  testthat::expect_error(setSubsetExperiment(es, "s1", noColnames),
                         regexp = "colnames")
  noRownames <- sub
  rownames(noRownames) <- NULL
  testthat::expect_error(setSubsetExperiment(es, "s1", noRownames),
                         regexp = "rownames")
})

testthat::test_that("Linked assays with the same values count as unchanged", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  # integer counts become double with the same values
  counts(sub) <- counts(sub) + 0
  testthat::expect_type(counts(sub), "double")
  es <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(subsetAssayNames(es), "s1")
})

# A 10 x 20 object with subsets "s1" (rows 1:5, cols 1:8), "s2" (rows 6:10,
# cols 9:20) and "s1a" (rows 2:4, cols 3:6 of s1), where s1 stores
# "doubled", a colData column "tag" and a reducedDim "PCA"
bracketExample <- function() {
  counts <- matrix(as.numeric(seq_len(200)), nrow = 10,
                   dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))
  es <- ExperimentSubset(SingleCellExperiment(list(counts = counts)))
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  es <- createSubset(es, "s2", rows = 6:10, cols = 9:20)
  es <- createSubset(es, "s1a", rows = 2:4, cols = 3:6, parentAssay = "s1")
  es <- setSubsetAssay(es, "s1", counts[1:5, 1:8] * 2, "doubled")
  colData(es, subsetName = "s1")$tag <- paste0("t", 1:8)
  pca <- matrix(as.numeric(1:16), nrow = 8,
                dimnames = list(paste0("c", 1:8), c("PC1", "PC2")))
  reducedDim(es, "PCA", subsetName = "s1") <- pca
  list(es = es, counts = counts, pca = pca)
}

testthat::test_that("[ keeps the rows of subsets that remain (#4)", {
  example <- bracketExample()
  es <- example$es
  counts <- example$counts

  kept <- es[3:10, ]
  testthat::expect_equal(as.matrix(assay(kept, "s1")), counts[3:5, 1:8])
  testthat::expect_equal(as.matrix(assay(kept, "s1/doubled")),
                         counts[3:5, 1:8] * 2)
  testthat::expect_equal(as.matrix(assay(kept, "s1a")), counts[3:4, 3:6])
  testthat::expect_equal(as.matrix(assay(kept, "s2")), counts[6:10, 9:20])
  testthat::expect_equal(subsetRownames(kept, "s1"), paste0("g", 3:5))

  # logical, negative and name subscripts select the same rows
  rowsKept <- rep(c(FALSE, TRUE), c(2, 8))
  testthat::expect_equal(assay(es[rowsKept, ], "s1a"), assay(kept, "s1a"))
  testthat::expect_equal(assay(es[-(1:2), ], "s1a"), assay(kept, "s1a"))
  testthat::expect_equal(assay(es[paste0("g", 3:10), ], "s1a"),
                         assay(kept, "s1a"))

  # with no subscripts nothing changes
  testthat::expect_equal(es[, ], es)

  # a subset keeps its own order when the main object is reordered
  reversed <- es[10:1, ]
  testthat::expect_equal(as.matrix(assay(reversed, "s1")), counts[1:5, 1:8])
  testthat::expect_equal(as.matrix(assay(reversed, "s1a")), counts[2:4, 3:6])
  testthat::expect_equal(as.matrix(assay(reversed[1:3, ], "s1")),
                         counts[character(0), 1:8, drop = FALSE])
})

testthat::test_that("[ keeps the columns of subsets that remain (#4)", {
  example <- bracketExample()
  es <- example$es
  counts <- example$counts

  kept <- es[, c(2, 5, 9, 10)]
  testthat::expect_equal(as.matrix(assay(kept, "s1")), counts[1:5, c(2, 5)])
  testthat::expect_equal(as.matrix(assay(kept, "s1/doubled")),
                         counts[1:5, c(2, 5)] * 2)
  testthat::expect_equal(as.matrix(assay(kept, "s1a")),
                         counts[2:4, 5, drop = FALSE])
  testthat::expect_equal(as.matrix(assay(kept, "s2")), counts[6:10, 9:10])
  testthat::expect_equal(colData(kept, subsetName = "s1")$tag, c("t2", "t5"))
  testthat::expect_equal(reducedDim(kept, "PCA", subsetName = "s1"),
                         example$pca[c(2, 5), ])

  both <- es[3:10, c(2, 5, 9, 10)]
  testthat::expect_equal(as.matrix(assay(both, "s1")), counts[3:5, c(2, 5)])
  testthat::expect_equal(subsetDim(both, "s1a"), c(2L, 1L))
})

testthat::test_that("[ keeps subsets left with no rows or columns (#4)", {
  example <- bracketExample()
  es <- example$es

  firstRows <- es[1:4, ]
  testthat::expect_equal(subsetNames(firstRows), subsetNames(es))
  testthat::expect_equal(dim(assay(firstRows, "s2")), c(0L, 12L))
  testthat::expect_equal(subsetDim(firstRows, "s2"), c(0L, 12L))
  testthat::expect_equal(dim(assay(firstRows, "s1a")), c(3L, 4L))

  lastCols <- es[, 9:20]
  testthat::expect_equal(dim(assay(lastCols, "s1")), c(5L, 0L))
  testthat::expect_equal(dim(assay(lastCols, "s1a")), c(3L, 0L))
  testthat::expect_equal(dim(assay(lastCols, "s2")), c(5L, 12L))
})

testthat::test_that("[ remaps subsets in every Subset class (#4)", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  se <- SummarizedExperiment(list(counts = counts))
  experiments <- list(
    se,
    as(se, "RangedSummarizedExperiment"),
    TreeSummarizedExperiment(assays = list(counts = counts),
                             rowTree = tinyTree),
    SpatialExperiment(assays = list(counts = counts),
                      sample_id = rep(c("A", "B"), each = 3),
                      spatialCoords = cbind(x = 1:6, y = 1:6))
  )
  for (experiment in experiments) {
    es <- ExperimentSubset(experiment)
    es <- createSubset(es, "s1", rows = 2:5, cols = 1:3)
    kept <- es[3:10, 2:6]
    testthat::expect_s4_class(kept, class(es))
    testthat::expect_equal(as.matrix(assay(kept, "s1")), counts[3:5, 2:3])
    # spatial coordinates of the subset follow its columns
    if (is(experiment, "SpatialExperiment"))
      testthat::expect_equal(unname(subsetSpatialCoords(kept, "s1")[, "x"]),
                             c(2, 3))
  }
})

testthat::test_that("[ keeps one copy of repeated rows in subsets (#4)", {
  example <- bracketExample()
  repeated <- example$es[c(3, 3, 4, 1), ]
  testthat::expect_equal(subsetRownames(repeated, "s1"),
                         paste0("g", c(1, 3, 4)))
  testthat::expect_equal(as.matrix(assay(repeated, "s1")),
                         example$counts[c(1, 3, 4), 1:8])
  testthat::expect_equal(as.matrix(assay(repeated, "s1a")),
                         example$counts[3:4, 3:6])
})

testthat::test_that("[ remaps subsets made from a stored assay (#4)", {
  example <- bracketExample()
  es <- createSubset(example$es, "fromDoubled", rows = 2:5, cols = 1:4,
                     parentAssay = "s1/doubled")
  kept <- es[3:10, 2:20]
  testthat::expect_equal(as.matrix(assay(kept, "fromDoubled")),
                         example$counts[3:5, 2:4] * 2)
})

testthat::test_that("createSubset says when the parent is empty (#68)", {
  example <- bracketExample()
  es <- example$es

  noRows <- es[6:10, ]
  testthat::expect_error(createSubset(noRows, "child", parentAssay = "s1"),
                         "Parent 's1' has no rows")
  testthat::expect_error(createSubset(noRows, "child", rows = 1,
                                      parentAssay = "s1"),
                         "Parent 's1' has no rows")
  testthat::expect_error(createSubset(noRows, "child", rows = "g1",
                                      parentAssay = "s1"),
                         "Parent 's1' has no rows")

  noCols <- es[, 9:20]
  testthat::expect_error(createSubset(noCols, "child", parentAssay = "s1"),
                         "Parent 's1' has no columns")
  testthat::expect_error(createSubset(noCols, "child", cols = 1,
                                      parentAssay = "s1"),
                         "Parent 's1' has no columns")
})

# Two objects to combine, each with a subset; the second also has a subset
# of its subset, an assay and a colData column stored in its subset. For
# cbind() they share rows, for rbind() they share columns
combineExample <- function(along) {
  counts1 <- matrix(as.numeric(seq_len(200)), nrow = 10,
                    dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))
  counts2 <- counts1 * 100
  if (along == "cols")
    colnames(counts2) <- paste0("d", 1:20)
  else
    rownames(counts2) <- paste0("h", 1:10)
  es1 <- ExperimentSubset(SingleCellExperiment(list(counts = counts1)))
  es1 <- createSubset(es1, "s1", rows = 1:5, cols = 1:8)
  es2 <- ExperimentSubset(SingleCellExperiment(list(counts = counts2)))
  es2 <- createSubset(es2, "s2", rows = 6:10, cols = 9:20)
  es2 <- createSubset(es2, "s2a", rows = 2:4, cols = 3:6, parentAssay = "s2")
  es2 <- setSubsetAssay(es2, "s2", counts2[6:10, 9:20] * 2, "doubled")
  colData(es2, subsetName = "s2")$tag <- paste0("t", 1:12)
  es2 <- createSubset(es2, "s2d", rows = 1:2, cols = 1:3,
                      parentAssay = "s2/doubled")
  list(es1 = es1, es2 = es2, counts1 = counts1, counts2 = counts2)
}

testthat::test_that("cbind keeps the subsets of every object (#67)", {
  example <- combineExample("cols")
  counts1 <- example$counts1
  counts2 <- example$counts2

  combined <- cbind(example$es1, example$es2)
  testthat::expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  testthat::expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  testthat::expect_equal(as.matrix(assay(combined, "s2")),
                         counts2[6:10, 9:20])
  testthat::expect_equal(subsetColnames(combined, "s2"), paste0("d", 9:20))
  # the subset's own assay and colData, and the subset of the subset
  testthat::expect_equal(as.matrix(assay(combined, "s2/doubled")),
                         counts2[6:10, 9:20] * 2)
  testthat::expect_equal(colData(combined, subsetName = "s2")$tag,
                         paste0("t", 1:12))
  testthat::expect_equal(as.matrix(assay(combined, "s2a")),
                         counts2[7:9, 11:14])
  # a subset made from an assay stored in a subset
  testthat::expect_equal(as.matrix(assay(combined, "s2d")),
                         counts2[6:7, 9:11] * 2)
})

testthat::test_that("rbind keeps the subsets of every object (#67)", {
  example <- combineExample("rows")
  counts1 <- example$counts1
  counts2 <- example$counts2

  combined <- rbind(example$es1, example$es2)
  testthat::expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  testthat::expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  testthat::expect_equal(as.matrix(assay(combined, "s2")),
                         counts2[6:10, 9:20])
  testthat::expect_equal(subsetRownames(combined, "s2"), paste0("h", 6:10))
  testthat::expect_equal(as.matrix(assay(combined, "s2/doubled")),
                         counts2[6:10, 9:20] * 2)
  testthat::expect_equal(colData(combined, subsetName = "s2")$tag,
                         paste0("t", 1:12))
  testthat::expect_equal(as.matrix(assay(combined, "s2a")),
                         counts2[7:9, 11:14])
  testthat::expect_equal(as.matrix(assay(combined, "s2d")),
                         counts2[6:7, 9:11] * 2)
})

testthat::test_that("cbind counts objects without subsets (#67)", {
  example <- combineExample("cols")
  counts3 <- example$counts1 * 10000
  colnames(counts3) <- paste0("e", 1:20)
  noSubsets <- ExperimentSubset(SingleCellExperiment(list(counts = counts3)))

  combined <- cbind(example$es1, noSubsets, example$es2)
  testthat::expect_equal(dim(combined), c(10L, 60L))
  testthat::expect_equal(as.matrix(assay(combined, "s2")),
                         example$counts2[6:10, 9:20])
  testthat::expect_equal(as.matrix(assay(combined, "s2a")),
                         example$counts2[7:9, 11:14])
})

testthat::test_that("cbind and rbind stop on repeated subset names (#67)", {
  example <- combineExample("cols")
  testthat::expect_error(cbind(example$es1, example$es1),
                         "unique.*'s1' is in more than one.*different names")
  rowsExample <- combineExample("rows")
  testthat::expect_error(rbind(rowsExample$es2, rowsExample$es2),
                         "'s2', 's2a', 's2d' are in more than one")
})

testthat::test_that("cbind and rbind stop when names don't line up (#67)", {
  example <- combineExample("cols")
  # the second object's genes in the opposite order
  reversed <- example$es2[10:1, ]
  testthat::expect_error(cbind(example$es1, reversed),
                         "same rows.*object 2")
  rowsExample <- combineExample("rows")
  reversed <- rowsExample$es2[, 20:1]
  # rbind() of the class the objects extend already checks the columns
  testthat::expect_error(rbind(rowsExample$es1, reversed),
                         "same colnames")
})

testthat::test_that("cbind updates the sample_id of spatial subsets (#67)", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  counts2 <- counts * 100
  colnames(counts2) <- paste0("other", 1:6)
  # both objects use samples A and B, which cbind renames
  spe1 <- SpatialExperiment(assays = list(counts = counts),
                            sample_id = rep(c("A", "B"), each = 3),
                            spatialCoords = cbind(x = 1:6, y = 1:6))
  spe2 <- SpatialExperiment(assays = list(counts = counts2),
                            sample_id = rep(c("A", "B"), each = 3),
                            spatialCoords = cbind(x = 1:6, y = 1:6))
  es1 <- createSubset(ExperimentSubset(spe1), "s1", cols = 2:4)
  es2 <- createSubset(ExperimentSubset(spe2), "s2", cols = 3:5)
  es2 <- createSubset(es2, "s2a", cols = 2:3, parentAssay = "s2")
  combined <- suppressMessages(cbind(es1, es2))
  for (name in c("s1", "s2", "s2a")) {
    rootCols <- match(subsetColnames(combined, name), colnames(combined))
    testthat::expect_equal(subsetColData(combined, name)$sample_id,
                           combined$sample_id[rootCols])
  }
  testthat::expect_false(identical(combined$sample_id[9:10],
                                   spe2$sample_id[3:4]))
})

testthat::test_that("All Subset classes keep subsets in cbind, rbind (#67)", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  counts2 <- counts * 100
  colnames(counts2) <- paste0("other", 1:6)
  makeExperiments <- function(counts, samples) {
    se <- SummarizedExperiment(list(counts = counts))
    list(
      se,
      as(se, "RangedSummarizedExperiment"),
      TreeSummarizedExperiment(assays = list(counts = counts),
                               rowTree = tinyTree),
      SpatialExperiment(assays = list(counts = counts),
                        sample_id = rep(samples, each = 3),
                        spatialCoords = cbind(x = 1:6, y = 1:6))
    )
  }
  first <- makeExperiments(counts, c("A", "B"))
  second <- makeExperiments(counts2, c("C", "D"))
  # for rbind, the same columns and samples, and the same rows of the tree
  third <- makeExperiments(counts * 100, c("A", "B"))
  for (k in seq_along(first)) {
    es1 <- createSubset(ExperimentSubset(first[[k]]), "s1",
                        rows = 2:5, cols = 1:3)
    es2 <- createSubset(ExperimentSubset(second[[k]]), "s2",
                        rows = 6:10, cols = 4:6)
    combined <- cbind(es1, es2)
    testthat::expect_s4_class(combined, class(es1))
    testthat::expect_equal(subsetNames(combined), c("s1", "s2"))
    testthat::expect_equal(as.matrix(assay(combined, "s2")),
                           counts2[6:10, 4:6])

    es3 <- createSubset(ExperimentSubset(third[[k]]), "s3",
                        rows = 6:10, cols = 4:6)
    stacked <- rbind(es1, es3)
    testthat::expect_s4_class(stacked, class(es1))
    testthat::expect_equal(subsetNames(stacked), c("s1", "s3"))
    testthat::expect_equal(as.matrix(assay(stacked, "s3")),
                           counts[6:10, 4:6] * 100)
  }
})

# Two objects for combineCols() (along = "cols") or combineRows() ("rows"),
# as in combineExample(), but they share only some rows (or columns), in a
# different order, so matching by name differs from matching by position
unionExample <- function(along) {
  example <- combineExample(along)
  counts2 <- example$counts2
  if (along == "cols")
    rownames(counts2) <- c(paste0("h", 1:4), paste0("g", 10:5))
  else
    colnames(counts2) <- c(paste0("e", 1:8), paste0("c", 20:9))
  es2 <- ExperimentSubset(SingleCellExperiment(list(counts = counts2)))
  es2 <- createSubset(es2, "s2", rows = 1:6, cols = 9:20)
  es2 <- createSubset(es2, "s2a", rows = 2:4, cols = 3:6, parentAssay = "s2")
  es2 <- setSubsetAssay(es2, "s2", counts2[1:6, 9:20] * 2, "doubled")
  colData(es2, subsetName = "s2")$tag <- paste0("t", 1:12)
  es2 <- createSubset(es2, "s2d", rows = 1:2, cols = 1:3,
                      parentAssay = "s2/doubled")
  list(es1 = example$es1, es2 = es2, counts1 = example$counts1,
       counts2 = counts2)
}

# Checks that combined holds the subsets of unionExample()'s objects
expectUnionSubsets <- function(combined, example) {
  counts1 <- example$counts1
  counts2 <- example$counts2
  testthat::expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  testthat::expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  testthat::expect_equal(as.matrix(assay(combined, "s2")),
                         counts2[1:6, 9:20])
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         rownames(counts2)[1:6])
  testthat::expect_equal(subsetColnames(combined, "s2"),
                         colnames(counts2)[9:20])
  # the subset's own assay and colData, and the subset of the subset
  testthat::expect_equal(as.matrix(assay(combined, "s2/doubled")),
                         counts2[1:6, 9:20] * 2)
  testthat::expect_equal(colData(combined, subsetName = "s2")$tag,
                         paste0("t", 1:12))
  testthat::expect_equal(as.matrix(assay(combined, "s2a")),
                         counts2[2:4, 11:14])
  # a subset made from an assay stored in a subset
  testthat::expect_equal(as.matrix(assay(combined, "s2d")),
                         counts2[1:2, 9:11] * 2)
}

testthat::test_that("combineCols keeps the subsets of every object (#73)", {
  example <- unionExample("cols")
  combined <- combineCols(example$es1, example$es2)
  testthat::expect_s4_class(combined, "SubsetSingleCellExperiment")
  testthat::expect_equal(nrow(combined), 14L)
  expectUnionSubsets(combined, example)
})

testthat::test_that("combineRows keeps the subsets of every object (#73)", {
  example <- unionExample("rows")
  combined <- combineRows(example$es1, example$es2)
  # combineRows() of SingleCellExperiment returns the class it extends
  testthat::expect_s4_class(combined, "SubsetRangedSummarizedExperiment")
  testthat::expect_equal(ncol(combined), 28L)
  expectUnionSubsets(combined, example)
})

testthat::test_that("combineCols, combineRows match by position (#73)", {
  example <- unionExample("cols")
  # the rows have different names, so the parent method warns that it
  # keeps the first object's rowRanges
  combined <- suppressWarnings(combineCols(example$es1, example$es2,
                                           use.names = FALSE))
  testthat::expect_equal(nrow(combined), 10L)
  # the rows keep their positions and take the first object's names
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         rownames(example$counts1)[1:6])
  testthat::expect_equal(unname(as.matrix(assay(combined, "s2"))),
                         unname(example$counts2[1:6, 9:20]))

  example <- unionExample("rows")
  combined <- combineRows(example$es1, example$es2, use.names = FALSE)
  testthat::expect_equal(ncol(combined), 20L)
  testthat::expect_equal(subsetColnames(combined, "s2"),
                         colnames(example$counts1)[9:20])
  testthat::expect_equal(unname(as.matrix(assay(combined, "s2"))),
                         unname(example$counts2[1:6, 9:20]))
})

testthat::test_that("combineCols counts objects without subsets (#73)", {
  example <- unionExample("cols")
  counts3 <- example$counts1 * 10000
  rownames(counts3) <- paste0("k", 1:10)
  colnames(counts3) <- paste0("e", 1:20)
  noSubsets <- ExperimentSubset(SingleCellExperiment(list(counts = counts3)))

  combined <- combineCols(example$es1, noSubsets, example$es2)
  testthat::expect_equal(dim(combined), c(24L, 60L))
  expectUnionSubsets(combined, example)
  # a plain object, not converted with ExperimentSubset()
  plain <- SingleCellExperiment(list(counts = counts3))
  combined <- combineCols(example$es1, plain, example$es2)
  testthat::expect_equal(dim(combined), c(24L, 60L))
  expectUnionSubsets(combined, example)
})

testthat::test_that("combineCols, combineRows stop on repeated names (#73)", {
  example <- unionExample("cols")
  testthat::expect_error(combineCols(example$es1, example$es1),
                         "unique.*'s1' is in more than one.*different names")
  rowsExample <- unionExample("rows")
  testthat::expect_error(combineRows(rowsExample$es2, rowsExample$es2),
                         "'s2', 's2a', 's2d' are in more than one")
})

testthat::test_that("All Subset classes keep subsets in combineCols (#73)", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  genes <- tinyTree$tip.label
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(genes, paste0("cell", 1:6)))
  # the second object shares 6 genes, in reverse order, and has 4 others
  counts2 <- counts * 100
  rownames(counts2) <- c(paste0("other", 1:4), rev(genes[5:10]))
  colnames(counts2) <- paste0("more", 1:6)
  # for combineRows, 3 shared cells in reverse order, and 3 others
  counts3 <- counts * 10000
  rownames(counts3) <- paste0("feature", 1:10)
  colnames(counts3) <- c(paste0("cell", 6:4), paste0("new", 1:3))
  makeExperiments <- function(counts, samples) {
    se <- SummarizedExperiment(list(counts = counts))
    list(
      se,
      as(se, "RangedSummarizedExperiment"),
      SingleCellExperiment(list(counts = counts)),
      TreeSummarizedExperiment(assays = list(counts = counts)),
      SpatialExperiment(assays = list(counts = counts),
                        sample_id = samples,
                        spatialCoords = cbind(x = 1:6, y = 1:6))
    )
  }
  first <- makeExperiments(counts, rep(c("A", "B"), each = 3))
  second <- makeExperiments(counts2, rep(c("C", "D"), each = 3))
  third <- makeExperiments(counts3, c("B", "B", "B", "E", "E", "E"))
  colsClass <- c("SubsetSummarizedExperiment",
                 "SubsetRangedSummarizedExperiment",
                 "SubsetSingleCellExperiment",
                 "SubsetSingleCellExperiment",
                 "SubsetSingleCellExperiment")
  rowsClass <- c("SubsetSummarizedExperiment",
                 rep("SubsetRangedSummarizedExperiment", 4))
  for (k in seq_along(first)) {
    es1 <- createSubset(ExperimentSubset(first[[k]]), "s1",
                        rows = 2:5, cols = 1:3)
    es2 <- createSubset(ExperimentSubset(second[[k]]), "s2",
                        rows = 1:6, cols = 4:6)
    combined <- combineCols(es1, es2)
    testthat::expect_s4_class(combined, colsClass[k])
    testthat::expect_equal(subsetNames(combined), c("s1", "s2"))
    testthat::expect_equal(as.matrix(assay(combined, "s1")), counts[2:5, 1:3])
    testthat::expect_equal(as.matrix(assay(combined, "s2")),
                           counts2[1:6, 4:6])
    for (name in subsetNames(combined)) {
      subset <- ExperimentSubset:::.subsets(combined)[[name]]
      internal <- ExperimentSubset:::.internalAssay(subset)
      testthat::expect_equal(class(internal)[1],
                             sub("Subset", "", colsClass[k]))
    }

    es3 <- createSubset(ExperimentSubset(third[[k]]), "s3",
                        rows = 6:10, cols = 1:4)
    stacked <- combineRows(es1, es3)
    testthat::expect_s4_class(stacked, rowsClass[k])
    testthat::expect_equal(subsetNames(stacked), c("s1", "s3"))
    testthat::expect_equal(as.matrix(assay(stacked, "s1")), counts[2:5, 1:3])
    testthat::expect_equal(as.matrix(assay(stacked, "s3")),
                           counts3[6:10, 1:4])
    for (name in subsetNames(stacked)) {
      subset <- ExperimentSubset:::.subsets(stacked)[[name]]
      internal <- ExperimentSubset:::.internalAssay(subset)
      testthat::expect_equal(class(internal)[1],
                             sub("Subset", "", rowsClass[k]))
    }
  }
})

testthat::test_that("Combining keeps names set with subsetRownames<- (#73)", {
  example <- combineExample("cols")
  es2 <- example$es2
  subsetRownames(es2, "s2") <- paste0("renamed", 1:5)
  combined <- cbind(example$es1, es2)
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         paste0("renamed", 1:5))

  example <- unionExample("cols")
  es2 <- example$es2
  ownNames <- c("custom", rownames(example$counts2)[2:6])
  subsetRownames(es2, "s2") <- ownNames
  combined <- combineCols(example$es1, es2)
  testthat::expect_equal(subsetRownames(combined, "s2"), ownNames)
  # matched by position, the renamed row keeps its name and the others
  # take the first object's
  combined <- suppressWarnings(combineCols(example$es1, es2,
                                           use.names = FALSE))
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         c("custom", rownames(example$counts1)[2:6]))

  example <- unionExample("rows")
  es2 <- example$es2
  subsetColnames(es2, "s2") <- c("custom", colnames(example$counts2)[10:20])
  combined <- combineRows(example$es1, es2, use.names = FALSE)
  testthat::expect_equal(subsetColnames(combined, "s2"),
                         c("custom", colnames(example$counts1)[10:20]))
})

testthat::test_that("Combining stops on a subset named as an assay (#73)", {
  example <- unionExample("cols")
  # the second object has a main assay named as the first object's subset
  es2 <- example$es2
  assay(es2, "s1") <- assay(es2, "counts") * 2
  testthat::expect_error(combineCols(example$es1, es2),
                         "'s1' is a subset in one object and an assay")
  example <- unionExample("rows")
  es1 <- example$es1
  assay(es1, "s2") <- assay(es1, "counts") * 2
  testthat::expect_error(combineRows(es1, example$es2),
                         "'s2' is a subset in one object and an assay")
})

testthat::test_that("use.names = FALSE names subsets of unnamed data (#73)", {
  example <- unionExample("cols")
  # the second object's genes have no names
  counts2 <- unname(example$counts2)
  colnames(counts2) <- paste0("d", 1:20)
  es2 <- ExperimentSubset(SingleCellExperiment(list(counts = counts2)))
  es2 <- createSubset(es2, "s2", rows = 1:6, cols = 9:20)
  combined <- suppressWarnings(combineCols(example$es1, es2,
                                           use.names = FALSE))
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         rownames(example$counts1)[1:6])

  # the first object's genes have no names, so the result takes the
  # second object's, and the subset keeps them
  counts1 <- unname(example$counts1)
  colnames(counts1) <- paste0("c", 1:20)
  es1 <- ExperimentSubset(SingleCellExperiment(list(counts = counts1)))
  combined <- suppressWarnings(combineCols(es1, example$es2,
                                           use.names = FALSE))
  testthat::expect_equal(subsetRownames(combined, "s2"),
                         rownames(combined)[1:6])
})

testthat::test_that("Combining allows a subset named as its own assay (#73)", {
  # an object saved before createSubset() checked names, with a subset
  # named as one of its own assays
  example <- unionExample("cols")
  es1 <- example$es1
  assay(es1, "extra") <- assay(es1, "counts") * 2
  s1 <- ExperimentSubset:::.subsets(es1)[["s1"]]
  legacy <- ExperimentSubset:::AssaySubset(
    subsetName = "extra",
    rowIndices = ExperimentSubset:::.rowIndices(s1),
    colIndices = ExperimentSubset:::.colIndices(s1),
    parentAssay = "counts",
    internalAssay = ExperimentSubset:::.internalAssay(s1)
  )
  subsets <- ExperimentSubset:::.subsets(es1)
  subsets[["extra"]] <- legacy
  ExperimentSubset:::.subsets(es1) <- subsets
  combined <- combineCols(es1, example$es2)
  testthat::expect_true("extra" %in% subsetNames(combined))
})

testthat::test_that("A subsetName that isn't one string gives a clear error", {
  # (#83) NULL, NA or several names failed inside .isSubset() with base R
  # errors that didn't name subsetName
  es <- ExperimentSubset(makeReducedDimObjects()$spe)
  es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
  es <- createSubset(es, "s2", rows = 1:5, cols = 2:4)
  msg <- "'subsetName' must be one subset name"
  pca <- matrix(as.numeric(seq_len(6)), nrow = 3)
  badValues <- list(NULL, c("s1", "s2"), NA, NA_character_)
  for (bad in badValues) {
    info <- deparse(bad)
    # methods that take subsetName through ...
    testthat::expect_error(metadata(es, subsetName = bad), msg, info = info)
    testthat::expect_error(metadata(es, subsetName = bad) <- list(),
                           msg, info = info)
    testthat::expect_error(altExp(es, subsetName = bad), msg, info = info)
    testthat::expect_error(altExp(es, "a", subsetName = bad) <- es,
                           msg, info = info)
    testthat::expect_error(altExps(es, subsetName = bad), msg, info = info)
    testthat::expect_error(altExps(es, subsetName = bad) <- list(),
                           msg, info = info)
    testthat::expect_error(altExpNames(es, subsetName = bad), msg,
                           info = info)
    testthat::expect_error(altExpNames(es, subsetName = bad) <- character(0),
                           msg, info = info)
    testthat::expect_error(reducedDim(es, subsetName = bad), msg, info = info)
    testthat::expect_error(reducedDim(es, "PCA", subsetName = bad) <- pca,
                           msg, info = info)
    testthat::expect_error(reducedDims(es, subsetName = bad), msg,
                           info = info)
    testthat::expect_error(reducedDims(es, subsetName = bad) <- list(),
                           msg, info = info)
    testthat::expect_error(reducedDimNames(es, subsetName = bad), msg,
                           info = info)
    testthat::expect_error(reducedDimNames(es, subsetName = bad) <- "X",
                           msg, info = info)
    testthat::expect_error(rowData(es, subsetName = bad), msg, info = info)
    testthat::expect_error(rowData(es, subsetName = bad) <- DataFrame(),
                           msg, info = info)
    testthat::expect_error(colData(es, subsetName = bad), msg, info = info)
    testthat::expect_error(colData(es, subsetName = bad) <- DataFrame(),
                           msg, info = info)
    testthat::expect_error(spatialCoords(es, subsetName = bad), msg,
                           info = info)
    testthat::expect_error(spatialData(es, subsetName = bad), msg,
                           info = info)
  }
  # NULL is subsetAssayNames()'s default, which lists every subset
  testthat::expect_true(all(c("s1", "s2") %in%
                              subsetAssayNames(es, subsetName = NULL)))
  for (bad in badValues[-1]) {
    info <- deparse(bad)
    testthat::expect_error(subsetAssayNames(es, subsetName = bad), msg,
                           info = info)
  }
  # (#97) these dispatched on subsetName = "character", so NULL, NA or a
  # number failed with "unable to find an inherited method"
  for (bad in c(badValues, list(1L))) {
    info <- deparse(bad)
    testthat::expect_error(subsetDim(es, bad), msg, info = info)
    testthat::expect_error(getSubsetAssay(es, bad), msg, info = info)
    testthat::expect_error(getSubsetAssay(es, bad, subsetAssayName = "x"),
                           msg, info = info)
    testthat::expect_error(subsetColData(es, bad) <- DataFrame(),
                           msg, info = info)
    testthat::expect_error(subsetRowData(es, bad) <- DataFrame(),
                           msg, info = info)
    testthat::expect_error(
      suppressWarnings(subsetSpatialData(es, bad)), msg, info = info
    )
    testthat::expect_error(subsetColnames(es, bad), msg, info = info)
    testthat::expect_error(subsetColnames(es, bad) <- letters[1:3],
                           msg, info = info)
    testthat::expect_error(subsetRownames(es, bad), msg, info = info)
    testthat::expect_error(subsetRownames(es, bad) <- letters[1:5],
                           msg, info = info)
    testthat::expect_error(subsetColData(es, bad), msg, info = info)
    testthat::expect_error(subsetRowData(es, bad), msg, info = info)
    testthat::expect_error(getSubsetExperiment(es, bad), msg, info = info)
    testthat::expect_error(setSubsetExperiment(es, bad, es), msg,
                           info = info)
    testthat::expect_error(subsetSpatialCoords(es, bad), msg, info = info)
    testthat::expect_error(setSubsetAssay(es, bad, pca, "a"), msg,
                           info = info)
  }
  # subsetParent() dispatches on x only, so it gets every bad value
  for (bad in badValues) {
    testthat::expect_error(subsetParent(es, bad), msg, info = deparse(bad))
  }
  # the links of a TreeSummarizedExperiment subset (#97)
  tes <- ExperimentSubset(makeReducedDimObjects()$tse)
  tes <- createSubset(tes, "s1", rows = 1:5, cols = 2:4)
  for (bad in c(badValues, list(1L))) {
    info <- deparse(bad)
    testthat::expect_error(subsetRowLinks(tes, bad), msg, info = info)
    testthat::expect_error(subsetColLinks(tes, bad), msg, info = info)
  }
  # on an object without subsets, that comes first, whatever subsetName is
  sce <- makeReducedDimObjects()$sce
  testthat::expect_error(reducedDimNames(sce, subsetName = NULL),
                         "which has no subsets")
  testthat::expect_error(reducedDimNames(sce, subsetName = NULL) <- "X",
                         "which has no subsets")
  testthat::expect_error(reducedDimNames(sce, NA_character_),
                         "which has no subsets")
})

testthat::test_that("An i that isn't one string gives a clear error (#85)", {
  # assay() and assay<- used i in if() without checking it, and
  # getSubsetAssay() without subsetAssayName passed subsetName on to them
  objects <- makeReducedDimObjects()
  se <- SummarizedExperiment(assays = list(counts = assay(objects$sce)))
  objects$se <- se
  objects$rse <- as(se, "RangedSummarizedExperiment")
  iMsg <- "'i' must be one name of an assay, a subset or a \"subset/assay\""
  subsetNameMsg <- "'subsetName' must be one subset name"
  m <- matrix(0, nrow = 5, ncol = 3)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    es <- createSubset(es, "s2", rows = 1:5, cols = 2:4)
    for (bad in list(c("s1", "s2"), c("counts", "s1"), NA_character_)) {
      info <- paste(name, deparse(bad))
      testthat::expect_error(assay(es, bad), iMsg, info = info)
      testthat::expect_error(assay(es, bad) <- m, iMsg, info = info)
      testthat::expect_error(getSubsetAssay(es, bad), subsetNameMsg,
                             info = info)
    }
    # an NA subsetAssayName became "s1/NA" and said that isn't an assay
    testthat::expect_error(
      assay(es, "s1", subsetAssayName = NA_character_),
      "'subsetAssayName' must be a single character string", info = name
    )
    # a numeric i still goes to the main object's assays
    testthat::expect_equal(assay(es, 1), assay(es, "counts"), info = name)
  }
})

testthat::test_that("ExperimentSubset() keeps an ExperimentSubset (#89)", {
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  inputs <- list(
    se = SummarizedExperiment(assays = list(counts = counts)),
    rse = as(SummarizedExperiment(assays = list(counts = counts)),
             "RangedSummarizedExperiment"),
    sce = SingleCellExperiment(assays = list(counts = counts)),
    tse = TreeSummarizedExperiment(assays = list(counts = counts)),
    spe = SpatialExperiment(assays = list(counts = counts))
  )
  for (name in names(inputs)) {
    es <- ExperimentSubset(inputs[[name]])
    testthat::expect_identical(ExperimentSubset(es), es, info = name)
  }

  # subset = still creates a subset on it
  es <- ExperimentSubset(inputs$sce)
  withSubset <- ExperimentSubset(es, subset = list(subsetName = "s1",
                                                   rows = 1:3,
                                                   cols = 1:4))
  testthat::expect_equal(subsetNames(withSubset), "s1")
  testthat::expect_equal(subsetDim(withSubset, subsetName = "s1"), c(3, 4))
})

testthat::test_that("Subset* subclasses store the base class (#89)", {
  setClass("Issue89SubsetSce", contains = "SubsetSingleCellExperiment")
  on.exit(removeClass("Issue89SubsetSce"), add = TRUE)
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  es <- ExperimentSubset(list(counts = counts))
  mine <- new("Issue89SubsetSce", es)

  # ExperimentSubset() returns it as is
  testthat::expect_identical(ExperimentSubset(mine), mine)

  # createSubset() stores a SingleCellExperiment, not the subclass
  mine <- createSubset(mine, "s1", rows = 1:3, cols = 1:4)
  stored <- getSubsetExperiment(mine, "s1")
  testthat::expect_equal(as.character(class(stored)), "SingleCellExperiment")

  # cbind() keeps the subsets, still stored as SingleCellExperiment
  other <- es
  colnames(other) <- paste0("other", 1:10)
  other <- createSubset(other, "s2", cols = 1:2)
  other <- new("Issue89SubsetSce", other)
  combined <- cbind(mine, other)
  testthat::expect_equal(subsetNames(combined), c("s1", "s2"))
  for (name in c("s1", "s2")) {
    stored <- getSubsetExperiment(combined, name)
    testthat::expect_equal(as.character(class(stored)),
                           "SingleCellExperiment", info = name)
  }

  # so do rbind() and combineCols()
  below <- es
  rownames(below) <- paste0("other", 1:6)
  below <- createSubset(below, "s3", rows = 1:2)
  below <- new("Issue89SubsetSce", below)
  combinations <- list(rbind = rbind(mine, below),
                       combineCols = combineCols(mine, other))
  for (how in names(combinations)) {
    combined <- combinations[[how]]
    for (name in subsetNames(combined)) {
      stored <- getSubsetExperiment(combined, name)
      testthat::expect_equal(as.character(class(stored)),
                             "SingleCellExperiment",
                             info = paste(how, name))
    }
  }
  testthat::expect_equal(subsetNames(combinations$rbind), c("s1", "s3"))
  testthat::expect_equal(subsetNames(combinations$combineCols),
                         c("s1", "s2"))
})

testthat::test_that("ExperimentSubset() converts subclasses (#89)", {
  setClass("Issue89Sce", contains = "SingleCellExperiment",
           slots = c(extra = "character"))
  on.exit(removeClass("Issue89Sce"), add = TRUE)
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  sce <- SingleCellExperiment(assays = list(counts = counts))
  mine <- new("Issue89Sce", sce, extra = "dropped")

  # converted as a SingleCellExperiment, with a warning
  testthat::expect_warning(
    es <- ExperimentSubset(mine),
    regexp = "'x' is an object of class Issue89Sce, a subclass of"
  )
  testthat::expect_equal(as.character(class(es)),
                         "SubsetSingleCellExperiment")
  testthat::expect_equal(assay(es, "counts"), counts)

  # subset = still creates a subset
  subset <- list(subsetName = "s1", rows = 1:3)
  testthat::expect_warning(es <- ExperimentSubset(mine, subset = subset),
                           regexp = "Issue89Sce")
  testthat::expect_equal(subsetNames(es), "s1")
})

testthat::test_that("ExperimentSubset() rejects an unsupported x (#89)", {
  counts <- matrix(seq_len(60), nrow = 6)
  testthat::expect_error(ExperimentSubset(counts),
                         regexp = "'x' must be a SummarizedExperiment")
})

testthat::test_that("unlink can name the assays to unlink (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  assay(es, "decontXcounts") <- simCounts * 0.5
  sub <- getSubsetExperiment(es, "s1")
  ownDecontX <- as.matrix(assay(sub, "counts")) * 0.9
  assay(sub, "decontXcounts") <- ownDecontX
  # a changed counts not named in unlink is an error, as with FALSE
  changedCounts <- sub
  assay(changedCounts, "counts") <- as.matrix(assay(sub, "counts")) + 1L
  # and says nothing about decontXcounts, which is not stored either
  messages <- testthat::capture_messages(testthat::expect_error(
    setSubsetExperiment(es, "s1", changedCounts, unlink = "decontXcounts"),
    regexp = "unlink = c(\"decontXcounts\", \"counts\")", fixed = TRUE
  ))
  testthat::expect_length(messages, 0)
  # only the named assay is unlinked, with a message saying so
  messages <- testthat::capture_messages(
    es2 <- setSubsetExperiment(es, "s1", sub, unlink = "decontXcounts")
  )
  testthat::expect_length(messages, 1)
  testthat::expect_match(messages, "'s1' now stores its own 'decontXcounts'")
  testthat::expect_equal(subsetAssayNames(es2, "s1"),
                         c("counts", "decontXcounts"))
  # s1 stores only decontXcounts; counts stays linked
  testthat::expect_equal(subsetAssayNames(es2),
                         c("s1", "s1/decontXcounts"))
  testthat::expect_equal(as.matrix(assay(es2, "s1/decontXcounts")),
                         ownDecontX)
  testthat::expect_equal(as.matrix(assay(es2, "s1/counts")),
                         simCounts[1:5, 1:8])
  testthat::expect_equal(assay(es2, "decontXcounts"), simCounts * 0.5)
  # names that are not linked or not changed are ignored
  messages <- testthat::capture_messages(
    es3 <- setSubsetExperiment(es, "s1", sub,
                               unlink = c("decontXcounts", "counts", "other"))
  )
  testthat::expect_length(messages, 1)
  testthat::expect_equal(subsetAssayNames(es3),
                         c("s1", "s1/decontXcounts"))
  testthat::expect_equal(as.matrix(assay(es3, "s1/counts")),
                         simCounts[1:5, 1:8])
  testthat::expect_equal(as.matrix(assay(es3, "s1/decontXcounts")),
                         ownDecontX)
  # character(0) unlinks nothing, as FALSE
  testthat::expect_error(
    setSubsetExperiment(es, "s1", sub, unlink = character(0)),
    regexp = "'decontXcounts' in subset 's1' is linked"
  )
  # every changed linked assay unlink doesn't allow is named at once
  assay(changedCounts, "decontXcounts") <- ownDecontX + 1
  testthat::expect_error(
    setSubsetExperiment(es, "s1", changedCounts, unlink = c(a = "other")),
    regexp = paste0("'counts' in subset 's1' is linked to the main assay ",
                    "'counts' and 'decontXcounts' to the main assay ",
                    "'decontXcounts'. Set unlink = TRUE or unlink = ",
                    "c(\"other\", \"counts\", \"decontXcounts\") to store ",
                    "separate 'counts' and 'decontXcounts'"),
    fixed = TRUE
  )
  # an error after the unlink check stores nothing and says nothing
  badName <- sub
  assay(badName, "a/b") <- ownDecontX
  messages <- testthat::capture_messages(testthat::expect_error(
    setSubsetExperiment(es, "s1", badName, unlink = TRUE),
    regexp = "can't contain '/'"
  ))
  testthat::expect_length(messages, 0)
})

testthat::test_that("unlink = TRUE says which assays it unlinks (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  assay(es, "decontXcounts") <- simCounts * 0.5
  sub <- getSubsetExperiment(es, "s1")
  assay(sub, "decontXcounts") <- as.matrix(assay(sub, "counts")) * 0.9
  assay(sub, "counts") <- as.matrix(assay(sub, "counts")) + 1L
  messages <- testthat::capture_messages(
    es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE)
  )
  countsMessage <- "'s1' now stores its own 'counts', unlinked from the main"
  testthat::expect_true(any(grepl(countsMessage, messages)))
  decontXMessage <- "'s1' now stores its own 'decontXcounts'"
  testthat::expect_true(any(grepl(decontXMessage, messages)))
  # storing an own assay again unlinks nothing, so says nothing
  testthat::expect_silent(
    es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE)
  )
  # a subset made from s1 names s1's assay as the source
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2, parentAssay = "s1")
  testthat::expect_message(
    es <- setSubsetAssay(es, "child", matrix(0, 2, 2), "counts",
                         unlink = TRUE),
    regexp = "'child' now stores its own 'counts', unlinked from 's1/counts'"
  )
})

testthat::test_that("setSubsetAssay and assay<- take unlink names (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  own <- simCounts[1:5, 1:8] + 1
  testthat::expect_error(
    setSubsetAssay(es, "s1", own, "counts", unlink = "decontXcounts"),
    regexp = "unlink = c(\"decontXcounts\", \"counts\")", fixed = TRUE
  )
  testthat::expect_message(
    es2 <- setSubsetAssay(es, "s1", own, "counts", unlink = "counts"),
    regexp = "'s1' now stores its own 'counts'"
  )
  testthat::expect_equal(as.matrix(assay(es2, "s1/counts")), own)
  testthat::expect_error(
    assay(es, "s1/counts", unlink = "decontXcounts") <- own,
    regexp = "is linked"
  )
  testthat::expect_message(assay(es, "s1/counts", unlink = "counts") <- own,
                           regexp = "now stores its own 'counts'")
  testthat::expect_equal(as.matrix(assay(es, "s1/counts")), own)
  # on a main assay, unlink names are an error, as TRUE is
  testthat::expect_error(
    assay(es, "counts", unlink = "counts") <- simCounts,
    regexp = "only applies"
  )
  # character(0) unlinks nothing, so is allowed there, as FALSE is
  assay(es, "counts", unlink = character(0)) <- simCounts + 1
  testthat::expect_equal(assay(es, "counts"), simCounts + 1)
})

testthat::test_that("An invalid unlink is an error (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  sub <- getSubsetExperiment(es, "s1")
  own <- simCounts[1:5, 1:8] + 1
  msg <- "'unlink' must be TRUE, FALSE or a character vector of assay names"
  for (bad in list(NA, c(TRUE, FALSE), 1, NA_character_, c("counts", ""))) {
    testthat::expect_error(setSubsetExperiment(es, "s1", sub, unlink = bad),
                           regexp = msg, fixed = TRUE)
    testthat::expect_error(setSubsetAssay(es, "s1", own, "x", unlink = bad),
                           regexp = msg, fixed = TRUE)
    testthat::expect_error(assay(es, "s1/x", unlink = bad) <- own,
                           regexp = msg, fixed = TRUE)
  }
})

# An object with 10 genes, 20 cells and an altExp "ADT" of 2 proteins,
# with a subset "s1" of the first 5 genes and 8 cells, as in #58
adtCounts <- matrix(seq_len(40), 2, 20,
                    dimnames = list(c("p1", "p2"), paste0("c", 1:20)))
makeAdtFixture <- function(class = "SingleCellExperiment") {
  counts <- matrix(seq_len(200) %% 7, 10, 20,
                   dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))
  constructor <- match.fun(class)
  object <- constructor(assays = list(counts = counts))
  altExp(object, "ADT") <- SummarizedExperiment(list(counts = adtCounts))
  es <- ExperimentSubset(object)
  createSubset(es, "s1", rows = 1:5, cols = 1:8, parentAssay = "counts")
}
# The part of an altExp a subset stores itself, NULL if none
ownAltExp <- function(es, subsetName, name) {
  internal <- ExperimentSubset:::.internalAssay(
    ExperimentSubset:::.subsets(es)[[subsetName]])
  if (!name %in% altExpNames(internal))
    return(NULL)
  altExp(internal, name, withDimnames = FALSE)
}

testthat::test_that("A subset sees the main object's altExps (#58)", {
  es <- makeAdtFixture()
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(dim(adt), c(2L, 8L))
  testthat::expect_equal(colnames(adt), paste0("c", 1:8))
  testthat::expect_equal(assay(adt, "counts"), adtCounts[, 1:8])
  # by position and in altExps(), in the same order
  testthat::expect_equal(
    assay(altExp(es, 1, subsetName = "s1"), "counts"), adtCounts[, 1:8])
  testthat::expect_equal(
    assay(altExp(es, subsetName = "s1"), "counts"), adtCounts[, 1:8])
  testthat::expect_equal(names(altExps(es, subsetName = "s1")), "ADT")
  # and in getSubsetExperiment()
  sub <- getSubsetExperiment(es, "s1")
  testthat::expect_equal(altExpNames(sub), "ADT")
  testthat::expect_equal(assay(altExp(sub, "ADT"), "counts"),
                         adtCounts[, 1:8])
  # a view: a later change to the main object's altExp shows in the subset
  altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts * 2L))
  testthat::expect_equal(
    assay(altExp(es, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] * 2L)
  # the subset's own altExps come after the linked ones
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es, "spikes", subsetName = "s1") <- spikes
  testthat::expect_equal(altExpNames(es, subsetName = "s1"),
                         c("ADT", "spikes"))
  testthat::expect_equal(nrow(altExp(es, 2, subsetName = "s1")), 3L)
  # an altExp the subset doesn't see gives SingleCellExperiment's error
  testthat::expect_error(altExp(es, "none", subsetName = "s1"),
                         regexp = "'none' not in 'altExpNames")
  testthat::expect_error(altExp(es, 3, subsetName = "s1"))
  none <- createSubset(ExperimentSubset(simSce), "s", cols = 1:3)
  testthat::expect_error(altExp(none, subsetName = "s"),
                         regexp = "no available entries")
})

testthat::test_that("withColData and withDimnames work on a subset's altExps (#58)", {
  es <- makeAdtFixture()
  colData(es, subsetName = "s1") <- S4Vectors::DataFrame(group = letters[1:8])
  adt <- altExp(es, "ADT", subsetName = "s1", withColData = TRUE)
  testthat::expect_equal(adt$group, letters[1:8])
  # the subset's colData columns are dropped again when stored
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1", withColData = TRUE) <- adt
  testthat::expect_equal(ncol(colData(ownAltExp(es, "s1", "ADT"))), 0L)
  # a value with other colnames is accepted only without withDimnames
  renamed <- adt
  colnames(renamed) <- paste0("x", 1:8)
  testthat::expect_warning(altExp(es, "ADT", subsetName = "s1") <- renamed,
                           regexp = "are not the same as")
  altExp(es, "ADT", subsetName = "s1", withDimnames = FALSE) <- renamed
  testthat::expect_equal(colnames(altExp(es, "ADT", subsetName = "s1")),
                         paste0("c", 1:8))
})

testthat::test_that("A subset stores only the altExp assays it adds (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  messages <- testthat::capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- adt)
  testthat::expect_length(messages, 0)
  testthat::expect_equal(
    assayNames(altExp(es, "ADT", subsetName = "s1")),
    c("counts", "decontPro"))
  testthat::expect_equal(assayNames(altExp(es, "ADT")), "counts")
  testthat::expect_equal(assayNames(ownAltExp(es, "s1", "ADT")),
                         "decontPro")
  # counts stays linked
  altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts + 1L))
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(assay(adt, "counts"), adtCounts[, 1:8] + 1L)
  testthat::expect_equal(assay(adt, "decontPro"), adtCounts[, 1:8] / 2)
  # storing what the subset sees stores nothing
  es <- makeAdtFixture()
  altExp(es, "ADT", subsetName = "s1") <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_null(ownAltExp(es, "s1", "ADT"))
})

testthat::test_that("A changed linked altExp assay needs unlink (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "counts") <- assay(adt, "counts") + 100L
  messages <- testthat::capture_messages(testthat::expect_error(
    altExp(es, "ADT", subsetName = "s1") <- adt,
    regexp = paste0("'ADT/counts' in subset 's1' is linked to the main ",
                    "altExp 'ADT'. Set unlink = TRUE or unlink = ",
                    "\"ADT/counts\" to store a separate 'ADT/counts'"),
    fixed = TRUE
  ))
  testthat::expect_length(messages, 0)
  testthat::expect_null(ownAltExp(es, "s1", "ADT"))
  # a bare name is a main assay, not the altExp's
  testthat::expect_error(
    altExp(es, "ADT", subsetName = "s1", unlink = "counts") <- adt,
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  for (unlink in list(TRUE, "ADT/counts")) {
    es2 <- es
    messages <- testthat::capture_messages(
      altExp(es2, "ADT", subsetName = "s1", unlink = unlink) <- adt)
    testthat::expect_equal(
      messages,
      paste0("Subset 's1' now stores its own 'ADT/counts', unlinked from ",
             "the main altExp 'ADT'.\n"))
    testthat::expect_equal(
      assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
      adtCounts[, 1:8] + 100L)
    testthat::expect_equal(assay(altExp(es2, "ADT"), "counts"), adtCounts)
  }
  # also through altExps<-
  testthat::expect_error(
    altExps(es, subsetName = "s1") <- list(ADT = adt),
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  # unlink only applies in a subset
  testthat::expect_error(
    altExp(es, "ADT", unlink = TRUE) <- altExp(es, "ADT"),
    regexp = "'unlink' only applies when storing in a subset.", fixed = TRUE)
  testthat::expect_error(
    altExp(es, "ADT", subsetName = "s1", unlink = NA) <- adt,
    regexp = "'unlink' must be TRUE, FALSE", fixed = TRUE)
})

testthat::test_that("An altExp of another shape replaces the linked one (#58)", {
  es <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  testthat::expect_error(
    altExp(es, "ADT", subsetName = "s1") <- small,
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  suppressMessages(
    altExp(es, "ADT", subsetName = "s1", unlink = "ADT/counts") <- small)
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(rownames(adt), "p9")
  testthat::expect_equal(assayNames(adt), "counts")
  testthat::expect_equal(assay(altExp(es, "ADT"), "counts"), adtCounts)
})

testthat::test_that("A subset's altExp rowData and metadata replace the linked ones (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  rowData(adt)$target <- c("CD3", "CD4")
  metadata(adt)$decontPro <- list(iterations = 10)
  altExp(es, "ADT", subsetName = "s1") <- adt
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(rowData(adt)$target, c("CD3", "CD4"))
  testthat::expect_equal(metadata(adt)$decontPro, list(iterations = 10))
  testthat::expect_equal(ncol(rowData(altExp(es, "ADT"))), 0L)
  own <- ownAltExp(es, "s1", "ADT")
  testthat::expect_length(assayNames(own), 0)
  testthat::expect_equal(colnames(rowData(own)), "target")
})

testthat::test_that("A chained subset sees altExps through its parent (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  spikes <- SummarizedExperiment(
    list(counts = matrix(seq_len(24), 3, 8,
                         dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es, "spikes", subsetName = "s1") <- spikes
  es <- createSubset(es, "s2", cols = 3:5, parentAssay = "s1")
  testthat::expect_equal(altExpNames(es, subsetName = "s2"),
                         c("ADT", "spikes"))
  view <- altExp(es, "ADT", subsetName = "s2")
  testthat::expect_equal(colnames(view), paste0("c", 3:5))
  testthat::expect_equal(assay(view, "counts"), adtCounts[, 3:5])
  testthat::expect_equal(assay(view, "decontPro"), adtCounts[, 3:5] / 2)
  testthat::expect_equal(
    as.vector(assay(altExp(es, "spikes", subsetName = "s2"), "counts")),
    as.vector(matrix(seq_len(24), 3, 8)[, 3:5]))
  # a change to an assay of the parent subset names where it is stored
  assay(view, "decontPro") <- assay(view, "decontPro") + 1
  testthat::expect_error(
    altExp(es, "ADT", subsetName = "s2") <- view,
    regexp = paste0("'ADT/decontPro' in subset 's2' is linked to altExp ",
                    "'ADT' of subset 's1'"),
    fixed = TRUE)
  messages <- testthat::capture_messages(
    altExp(es, "ADT", subsetName = "s2", unlink = TRUE) <- view)
  testthat::expect_equal(
    messages,
    paste0("Subset 's2' now stores its own 'ADT/decontPro', unlinked from ",
           "altExp 'ADT' of subset 's1'.\n"))
})

testthat::test_that("A subset removes and renames only its own altExps (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  # NULL removes the subset's part; the linked altExp shows again
  altExp(es, "ADT", subsetName = "s1") <- NULL
  testthat::expect_null(ownAltExp(es, "s1", "ADT"))
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  testthat::expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
                         "counts")
  # altExps<- can't remove a linked altExp
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExps(es, subsetName = "s1") <- list(spikes = spikes)
  testthat::expect_equal(altExpNames(es, subsetName = "s1"),
                         c("ADT", "spikes"))
  altExps(es, subsetName = "s1") <- list()
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  # renaming: the subset's own altExps only
  altExp(es, "spikes", subsetName = "s1") <- spikes
  altExpNames(es, subsetName = "s1") <- c("ADT", "ercc")
  testthat::expect_equal(altExpNames(es, subsetName = "s1"),
                         c("ADT", "ercc"))
  testthat::expect_equal(altExpNames(es), "ADT")
  testthat::expect_error(
    altExpNames(es, subsetName = "s1") <- c("protein", "ercc"),
    regexp = paste0("Subset 's1' can't rename 'ADT', which is linked to ",
                    "the main altExp 'ADT'."),
    fixed = TRUE)
  # names must differ, so an own altExp can't take a linked name, which
  # would hide the linked altExp without unlink
  testthat::expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", "ADT"),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  testthat::expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", NA),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  testthat::expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", "ercc", "x"),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  testthat::expect_equal(altExpNames(es, subsetName = "s1"),
                         c("ADT", "ercc"))
})

testthat::test_that("altExpNames<- rejects values that aren't character", {
  # (#105) other values failed at S4 dispatch, before the package's check
  valueMsg <- "'value' must be a character vector of altExp names."
  classes <- c("SingleCellExperiment", "TreeSummarizedExperiment",
               "SpatialExperiment")
  for (class in classes) {
    es <- makeAdtFixture(class)
    # a factor of the right length, NULL and a number
    testthat::expect_error(
      altExpNames(es, subsetName = "s1") <- factor("protein"),
      valueMsg, fixed = TRUE, info = class
    )
    testthat::expect_error(altExpNames(es, subsetName = "s1") <- NULL,
                           valueMsg, fixed = TRUE, info = class)
    testthat::expect_error(altExpNames(es, subsetName = "s1") <- 1L,
                           valueMsg, fixed = TRUE, info = class)
    testthat::expect_error(altExpNames(es) <- factor("protein"),
                           valueMsg, fixed = TRUE, info = class)
    testthat::expect_error(altExpNames(es) <- NULL,
                           valueMsg, fixed = TRUE, info = class)
    testthat::expect_error(altExpNames(es) <- 1L,
                           valueMsg, fixed = TRUE, info = class)
    # nothing was renamed
    testthat::expect_equal(altExpNames(es), "ADT", info = class)
    testthat::expect_equal(altExpNames(es, subsetName = "s1"), "ADT",
                           info = class)
    # a character value still renames the main altExp
    altExpNames(es) <- "protein"
    testthat::expect_equal(altExpNames(es), "protein", info = class)
    testthat::expect_equal(altExpNames(es, subsetName = "s1"), "protein",
                           info = class)
  }
})

testthat::test_that("Writes after an unlink don't need unlink again (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "counts") <- assay(adt, "counts") + 100L
  suppressMessages(altExp(es, "ADT", subsetName = "s1", unlink = TRUE) <- adt)
  # s1's own counts is no longer linked
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  messages <- testthat::capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- adt)
  testthat::expect_length(messages, 0)
  testthat::expect_equal(assayNames(ownAltExp(es, "s1", "ADT")),
                         c("counts", "decontPro"))
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  messages <- testthat::capture_messages(
    altExp(es, "spikes", subsetName = "s1") <- spikes)
  testthat::expect_length(messages, 0)
  messages <- testthat::capture_messages(
    es2 <- setSubsetExperiment(es, "s1", getSubsetExperiment(es, "s1")))
  testthat::expect_length(messages, 0)
  testthat::expect_equal(
    assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] + 100L)
  # also after storing an altExp of another shape
  es <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  suppressMessages(altExp(es, "ADT", subsetName = "s1", unlink = TRUE) <- small)
  assay(small, "counts") <- assay(small, "counts") + 1L
  messages <- testthat::capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- small)
  testthat::expect_length(messages, 0)
  testthat::expect_equal(
    as.vector(assay(altExp(es, "ADT", subsetName = "s1"), "counts")),
    rep(2L, 8))
})

testthat::test_that("A subset's altExp reducedDims replace the linked ones (#58)", {
  es <- makeAdtFixture()
  adtSce <- SingleCellExperiment(list(counts = adtCounts))
  reducedDim(adtSce, "PCA") <- matrix(seq_len(40), 20, 2)
  altExp(es, "ADT") <- adtSce
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(unname(reducedDim(adt, "PCA")),
                         matrix(seq_len(40), 20, 2)[1:8, ])
  # storing the same reducedDims stores nothing
  altExp(es, "ADT", subsetName = "s1") <- adt
  testthat::expect_null(ownAltExp(es, "s1", "ADT"))
  reducedDim(adt, "PCA") <- matrix(0, 8, 2)
  metadata(adt)$note <- "own"
  altExp(es, "ADT", subsetName = "s1") <- adt
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(unname(reducedDim(adt, "PCA")), matrix(0, 8, 2))
  testthat::expect_equal(metadata(adt)$note, "own")
  testthat::expect_equal(unname(reducedDim(altExp(es, "ADT"), "PCA")),
                         matrix(seq_len(40), 20, 2))
  testthat::expect_length(assayNames(ownAltExp(es, "s1", "ADT")), 0)
  # an altExp without reducedDims shows the linked ones and stores nothing
  es <- makeAdtFixture()
  altExp(es, "ADT") <- adtSce
  adt <- altExp(es, "ADT", subsetName = "s1")
  reducedDims(adt) <- list()
  altExp(es, "ADT", subsetName = "s1") <- adt
  testthat::expect_null(ownAltExp(es, "s1", "ADT"))
})

testthat::test_that("A rename carries what subsets store along (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  es <- createSubset(es, "s2", cols = 1:3, parentAssay = "s1")
  adt <- altExp(es, "ADT", subsetName = "s2")
  assay(adt, "zero") <- assay(adt, "counts") * 0L
  altExp(es, "ADT", subsetName = "s2") <- adt
  altExpNames(es) <- "protein"
  testthat::expect_equal(altExpNames(es, subsetName = "s1"), "protein")
  testthat::expect_equal(
    assayNames(altExp(es, "protein", subsetName = "s2")),
    c("counts", "decontPro", "zero"))
  # a parent subset's own altExp: its children's parts follow it
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es, "spikes", subsetName = "s1") <- spikes
  child <- altExp(es, "spikes", subsetName = "s2")
  assay(child, "logcounts") <- log1p(assay(child, "counts"))
  altExp(es, "spikes", subsetName = "s2") <- child
  altExpNames(es, subsetName = "s1") <- c("protein", "ercc")
  testthat::expect_equal(altExpNames(es, subsetName = "s2"),
                         c("protein", "ercc"))
  testthat::expect_equal(assayNames(altExp(es, "ercc", subsetName = "s2")),
                         c("counts", "logcounts"))
  # a rename onto a name a child already stores is an error
  other <- SummarizedExperiment(
    list(counts = matrix(1L, 1, 3, dimnames = list(NULL, paste0("c", 1:3)))))
  altExp(es, "other", subsetName = "s2") <- other
  testthat::expect_error(
    altExpNames(es, subsetName = "s1") <- c("protein", "other"),
    regexp = "Subset 's2' already stores an altExp called 'other'.",
    fixed = TRUE)
  # also when the subset stores only the new name, an altExp of its own
  es <- makeAdtFixture()
  own <- SummarizedExperiment(
    list(counts = matrix(1L, 2, 8, dimnames = list(c("p1", "p2"),
                                                   paste0("c", 1:8)))))
  altExp(es, "X", subsetName = "s1") <- own
  testthat::expect_error(
    altExpNames(es) <- "X",
    regexp = "Subset 's1' already stores an altExp called 'X'.",
    fixed = TRUE)
  testthat::expect_equal(altExpNames(es), "ADT")
})

testthat::test_that("A subset's whole altExp stays whole (#58)", {
  es <- makeAdtFixture()
  twoRows <- function(value, assays = "counts") {
    m <- matrix(value, 2, 20, dimnames = list(c("r1", "r2"),
                                             paste0("c", 1:20)))
    SummarizedExperiment(stats::setNames(rep(list(m), length(assays)),
                                         assays))
  }
  own <- twoRows(5L)[, 1:8]
  suppressMessages(altExp(es, "ADT", subsetName = "s1", unlink = TRUE) <- own)
  # the main ADT takes the same rows: s1's ADT still replaces it
  testthat::expect_warning(
    altExp(es, "ADT") <- twoRows(1L, c("counts", "other")), regexp = NA)
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(assayNames(adt), "counts")
  testthat::expect_equal(as.vector(assay(adt, "counts")), rep(5L, 16))
  # and is kept when the main ADT is reshaped
  testthat::expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts)),
    regexp = NA)
  testthat::expect_equal(rownames(altExp(es, "ADT", subsetName = "s1")),
                         c("r1", "r2"))
  # and stays whole through later writes of other altExps and round trips
  altExp(es, "ADT") <- twoRows(1L, c("counts", "other"))
  altExp(es, "B", subsetName = "s1") <- own
  es <- setSubsetExperiment(es, "s1", getSubsetExperiment(es, "s1"))
  testthat::expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
                         "counts")
  testthat::expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts)),
    regexp = NA)
  testthat::expect_equal(rownames(altExp(es, "ADT", subsetName = "s1")),
                         c("r1", "r2"))
  # an altExp the main object adds under a name s1 already stores is
  # hidden by s1's
  altExp(es, "X", subsetName = "s1") <- own
  altExp(es, "X") <- twoRows(1L, c("counts", "other"))
  testthat::expect_equal(assayNames(altExp(es, "X", subsetName = "s1")),
                         "counts")
})

testthat::test_that("Objects saved before pieces were recorded read as before (#58)", {
  es <- makeAdtFixture()
  # a full copy stored as earlier versions did, with no record of pieces
  subsets <- ExperimentSubset:::.subsets(es)
  internal <- ExperimentSubset:::.internalAssay(subsets[["s1"]])
  copy <- SummarizedExperiment(list(counts = adtCounts[, 1:8] + 5L))
  altExp(internal, "ADT") <- copy
  subsets[["s1"]] <- ExperimentSubset:::`.internalAssay<-`(
    subsets[["s1"]], value = internal)
  es <- ExperimentSubset:::`.subsets<-`(es, value = subsets)
  altExp(es, "ADT") <- SummarizedExperiment(
    list(counts = adtCounts, other = adtCounts * 2L))
  # the copy's counts hide the main ADT's, whose other assay shows too
  adt <- altExp(es, "ADT", subsetName = "s1")
  testthat::expect_equal(assayNames(adt), c("counts", "other"))
  testthat::expect_equal(assay(adt, "counts"), adtCounts[, 1:8] + 5L)
})

testthat::test_that("combineCols carries what subsets store along (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  # a plain object of 3 cells whose ADT has the given proteins
  makeOther <- function(proteins) {
    counts <- matrix(1L, 10, 3, dimnames = list(paste0("g", 1:10),
                                               paste0("n", 1:3)))
    other <- SingleCellExperiment(list(counts = counts))
    altExp(other, "ADT") <- SummarizedExperiment(list(
      counts = matrix(1L, length(proteins), 3,
                      dimnames = list(proteins, paste0("n", 1:3)))))
    other
  }
  # the same proteins: nothing changes
  testthat::expect_warning(
    out <- combineCols(es, makeOther(c("p1", "p2"))), regexp = NA)
  testthat::expect_equal(assayNames(altExp(out, "ADT", subsetName = "s1")),
                         c("counts", "decontPro"))
  # a third protein: ADT gets another row, and s1's part is dropped
  testthat::expect_warning(
    out <- combineCols(es, makeOther(c("p1", "p2", "p3"))),
    regexp = "Dropped what subset 's1' stored in altExp 'ADT'", fixed = TRUE)
  testthat::expect_equal(dim(altExp(out, "ADT", subsetName = "s1")),
                         c(3L, 8L))
  testthat::expect_equal(assayNames(altExp(out, "ADT", subsetName = "s1")),
                         "counts")
})

testthat::test_that("Removing or reshaping an altExp drops what subsets added to it (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  # new values keep the subset's part
  testthat::expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts + 1L)),
    regexp = NA)
  testthat::expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
                         c("counts", "decontPro"))
  dropped <- paste0("Dropped what subset 's1' stored in altExp 'ADT': the ",
                    "altExp it added to was removed, or now has other rows ",
                    "or another class.")
  # other rows
  es2 <- es
  testthat::expect_warning(altExp(es2, "ADT") <- altExp(es2, "ADT")[1, ],
                           regexp = dropped, fixed = TRUE)
  testthat::expect_null(ownAltExp(es2, "s1", "ADT"))
  testthat::expect_equal(dim(altExp(es2, "ADT", subsetName = "s1")),
                         c(1L, 8L))
  # removed, one or all
  es2 <- es
  testthat::expect_warning(altExp(es2, "ADT") <- NULL,
                           regexp = dropped, fixed = TRUE)
  testthat::expect_length(altExpNames(es2, subsetName = "s1"), 0)
  es2 <- es
  testthat::expect_warning(altExps(es2) <- list(),
                           regexp = dropped, fixed = TRUE)
  testthat::expect_length(altExpNames(es2, subsetName = "s1"), 0)
  # a whole altExp of the subset's own is kept
  es3 <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  suppressMessages(altExp(es3, "ADT", subsetName = "s1", unlink = TRUE) <- small)
  testthat::expect_warning(altExp(es3, "ADT") <- NULL, regexp = NA)
  testthat::expect_equal(rownames(altExp(es3, "ADT", subsetName = "s1")),
                         "p9")
  # a parent subset's change reaches its children, also through
  # setSubsetExperiment()
  es4 <- makeAdtFixture()
  es4 <- createSubset(es4, "s2", cols = 1:3, parentAssay = "s1")
  adt <- altExp(es4, "ADT", subsetName = "s2")
  assay(adt, "zero") <- assay(adt, "counts") * 0L
  altExp(es4, "ADT", subsetName = "s2") <- adt
  testthat::expect_warning(
    suppressMessages(
      altExp(es4, "ADT", subsetName = "s1", unlink = TRUE) <- small),
    regexp = "Dropped what subset 's2' stored in altExp 'ADT'", fixed = TRUE)
  es5 <- makeAdtFixture()
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es5, "spikes", subsetName = "s1") <- spikes
  es5 <- createSubset(es5, "s2", cols = 1:3, parentAssay = "s1")
  child <- altExp(es5, "spikes", subsetName = "s2")
  assay(child, "logcounts") <- log1p(assay(child, "counts"))
  altExp(es5, "spikes", subsetName = "s2") <- child
  sub <- getSubsetExperiment(es5, "s1")
  altExp(sub, "spikes") <- NULL
  testthat::expect_warning(
    es5 <- setSubsetExperiment(es5, "s1", sub),
    regexp = "Dropped what subset 's2' stored in altExp 'spikes'",
    fixed = TRUE)
  testthat::expect_equal(altExpNames(es5, subsetName = "s2"), "ADT")
})

testthat::test_that("setSubsetExperiment stores only changed altExp assays (#58)", {
  es <- makeAdtFixture()
  sub <- getSubsetExperiment(es, "s1")
  # a round trip stores no copy of ADT
  es2 <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_null(ownAltExp(es2, "s1", "ADT"))
  testthat::expect_equal(altExpNames(es2, subsetName = "s1"), "ADT")
  # a new assay in ADT is stored without unlink
  assay(altExp(sub, "ADT"), "decontPro") <- adtCounts[, 1:8] / 2
  es2 <- setSubsetExperiment(es, "s1", sub)
  testthat::expect_equal(assayNames(ownAltExp(es2, "s1", "ADT")),
                         "decontPro")
  # changed linked assays, main and altExp, are named in one error
  changed <- sub
  assay(changed, "counts") <- as.matrix(assay(changed, "counts")) + 1
  assay(altExp(changed, "ADT"), "counts") <- adtCounts[, 1:8] + 1L
  messages <- testthat::capture_messages(testthat::expect_error(
    setSubsetExperiment(es, "s1", changed),
    regexp = paste0("'counts' in subset 's1' is linked to the main assay ",
                    "'counts' and 'ADT/counts' to the main altExp 'ADT'. ",
                    "Set unlink = TRUE or unlink = c(\"counts\", ",
                    "\"ADT/counts\")"),
    fixed = TRUE))
  testthat::expect_length(messages, 0)
  testthat::expect_error(
    setSubsetExperiment(es, "s1", changed, unlink = "counts"),
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  messages <- testthat::capture_messages(
    es2 <- setSubsetExperiment(es, "s1", changed, unlink = TRUE))
  testthat::expect_length(messages, 2)
  testthat::expect_match(messages[2], "its own 'ADT/counts'", fixed = TRUE)
  testthat::expect_equal(
    assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] + 1L)
})

testthat::test_that("TreeSE and SpatialExperiment subsets see altExps (#58)", {
  for (class in c("TreeSummarizedExperiment", "SpatialExperiment")) {
    es <- makeAdtFixture(class)
    testthat::expect_equal(altExpNames(es, subsetName = "s1"), "ADT",
                           info = class)
    testthat::expect_equal(
      assay(altExp(es, "ADT", subsetName = "s1"), "counts"),
      adtCounts[, 1:8], info = class)
    sub <- getSubsetExperiment(es, "s1")
    testthat::expect_equal(altExpNames(sub), "ADT", info = class)
  }
})
