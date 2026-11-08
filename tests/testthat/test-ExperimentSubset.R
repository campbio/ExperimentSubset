library(testthat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ExperimentSubset)
# most tests use sce_chcl from scds (Suggests); skip the file without it
skip_if_not_installed("scds")

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
  expect_equal(subsetRownames(es, subsetName), rownames(es)[rows])
  expect_equal(subsetColnames(es, subsetName), colnames(es)[cols])
  expect_equal(subsetParent(es, subsetName), list("counts"))
  expect_equal(assay(es, subsetName),
               assay(es, "counts")[rows, cols, drop = FALSE])
}

test_that("Testing ExperimentSubset constructor by implicitly providing a SingleCellExperiment object",{
  es <- ExperimentSubset(sce_chcl)
  expect_true(validObject(es))
  expect_equal(class(es)[1], "SubsetSingleCellExperiment")
})

test_that("Testing ExperimentSubset constructor by implicitly providing a SummarizedExperiment object",{
  se <- SummarizedExperiment::SummarizedExperiment(sce_chcl)
  expect_error(es <- ExperimentSubset(se), "Experiment objects with unnamed assays are not supported. Use assayNames<- setter method to set assay names before creating ES object.")
  se <- SummarizedExperiment::SummarizedExperiment(list(counts = assay(sce_chcl, "counts"), logcounts = assay(sce_chcl, "logcounts")))
  expect_s4_class(es <- ExperimentSubset(se), "SubsetSummarizedExperiment")
  expect_true(validObject(es))
  expect_equal(length(metadata(es)), 0)
  expect_s4_class(assay(es, "counts"), "dgCMatrix")
  es <- createSubset(es, "s1", rows = c(1:10))
  expect_true("s1" %in% subsetNames(es))
  expect_true("s1" %in% subsetAssayNames(es))
  expectSubsetOf(es, "s1", rows = 1:10, cols = seq_len(ncol(se)))
  expect_output(show(es), "subsets\\(1\\): s1")
  s1Plus <- assay(es, "s1") + 1
  es <- setSubsetAssay(es, "s1", s1Plus, "newS1")
  expect_true("s1/newS1" %in% subsetAssayNames(es))
  expect_false("newS1" %in% subsetNames(es))
  expect_equal(getSubsetAssay(es, "s1", subsetAssayName = "newS1"),
               s1Plus)
  expect_output(show(es), "subsetAssays\\(2\\): s1 s1/newS1")
  assay(es, "assay2") <- assay(es, "counts")
  expect_true("assay2" %in% assayNames(es))
  expect_equal(assay(es, "assay2"), assay(es, "counts"))
  expect_equal(subsetParent(es, "s1"), list("counts"))
  expect_equal(subsetParent(es, "s1/newS1"), list("s1", "counts"))
  expect_output(subsetSummary(es), "s1 +10, 2000 +counts +newS1")
  expect_equal(subsetDim(es, "s1")[1], 10)
  expect_equal(subsetDim(es, "s1")[2], 2000)
  expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  expect_equal(length(metadata(es)), 1)
  expect_s4_class(getSubsetAssay(es, "s1"), "dgCMatrix")
  expect_equal(getSubsetAssay(es, "s1"),
               assay(se, "counts")[1:10, ])
  subsetRownames(es, "s1") <- subsetRownames(es, "s1")
  expect_equal(subsetRownames(es, "s1"), rownames(se)[1:10])
  subsetColnames(es, "s1") <- subsetColnames(es, "s1")
  expect_equal(subsetColnames(es, "s1"), colnames(se))
  subsetRownames(es, "s1/newS1") <- subsetRownames(es, "s1/newS1")
  expect_equal(subsetRownames(es, "s1/newS1"), rownames(se)[1:10])
  subsetColnames(es, "s1/newS1") <- subsetColnames(es, "s1/newS1")
  expect_equal(subsetColnames(es, "s1/newS1"), colnames(se))
  colData(es, subsetName = "s1") <- cbind(colData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[2]))))
  expect_equal(ncol(colData(es, subsetName = "s1")), 1)
  expect_equal(colData(es, subsetName = "s1")$col1, seq_len(2000))
  rowData(es, subsetName = "s1") <- cbind(rowData(es, subsetName = "s1"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s1")[1]))))
  expect_equal(ncol(rowData(es, subsetName = "s1")), 1)
  expect_equal(rowData(es, subsetName = "s1")$col1, seq_len(10))
  
  
  skip_if_not_installed("airway")
  data(airway, package="airway")
  rse <- airway
  expect_s4_class(es <- ExperimentSubset(rse), "SubsetRangedSummarizedExperiment")
  expect_output(subsetSummary(es),
                "Main assay\\(s\\):\n counts.*Subset\\(s\\):\nNULL")
  expect_equal(length(metadata(es)), 1)

  es <- createSubset(es, "s2", rows = c(1:15))
  expect_true("s2" %in% subsetNames(es))
  expect_true("s2" %in% subsetAssayNames(es))
  expectSubsetOf(es, "s2", rows = 1:15, cols = seq_len(ncol(rse)))
  expect_output(show(es), "subsets\\(1\\): s2")
  s2Plus <- assay(es, "s2") + 1
  es <- setSubsetAssay(es, "s2", s2Plus, "newS2")
  expect_true("s2/newS2" %in% subsetAssayNames(es))
  expect_false("newS2" %in% subsetNames(es))
  expect_equal(getSubsetAssay(es, "s2", subsetAssayName = "newS2"),
               s2Plus)
  assay(es, "assay2") <- assay(es, "counts")
  expect_true("assay2" %in% assayNames(es))
  expect_equal(assay(es, "assay2"), assay(es, "counts"))
  expect_equal(subsetParent(es, "s2"), list("counts"))
  expect_output(subsetSummary(es), "s2 +15, 8 +counts +newS2")
  expect_equal(subsetDim(es, "s2")[1], 15)
  expect_equal(subsetDim(es, "s2")[2], 8)
  expect_equal(subsetCount(es), 1)
  metadata(es) <- list(a="meta")
  expect_equal(length(metadata(es)), 1)

  subsetRownames(es, "s2") <- subsetRownames(es, "s2")
  expect_equal(subsetRownames(es, "s2"), rownames(rse)[1:15])
  subsetColnames(es, "s2") <- subsetColnames(es, "s2")
  expect_equal(subsetColnames(es, "s2"), colnames(rse))
  subsetRownames(es, "s2/newS2") <- subsetRownames(es, "s2/newS2")
  expect_equal(subsetRownames(es, "s2/newS2"),
               rownames(rse)[1:15])
  subsetColnames(es, "s2/newS2") <- subsetColnames(es, "s2/newS2")
  expect_equal(subsetColnames(es, "s2/newS2"), colnames(rse))
  colData(es, subsetName = "s2") <- cbind(colData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[2]))))
  expect_equal(ncol(colData(es, subsetName = "s2", parentColData = TRUE)), 10)
  expect_equal(colData(es, subsetName = "s2")$col1, seq_len(8))
  rowData(es, subsetName = "s2") <- cbind(rowData(es, subsetName = "s2"), 
                                          DataFrame(col1 = c(seq(subsetDim(es, subsetName = "s2")[1]))))
  expect_equal(ncol(rowData(es, subsetName = "s2")), 1)
  expect_equal(rowData(es, subsetName = "s2")$col1, seq_len(15))
  
  # compare the whole message, so it is one line with no runs of spaces
  expect_equal(
    tryCatch(createSubset(es, "s2", rows = c(1:15)),
             error = conditionMessage),
    paste0("'subsetName' must differ from the names of existing subsets; ",
           "'s2' is already a subset name.")
  )
  expect_error(es <- createSubset(es, subsetName = 12, rows = c(1:15)),
               "'subsetName' must be a single, non-empty")
  expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", rows = list()),
             error = conditionMessage),
    paste0("'rows' must be NULL, to include all rows, or a numeric or ",
           "character vector of the rows to include.")
  )
  expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", cols = list()),
             error = conditionMessage),
    paste0("'cols' must be NULL, to include all columns, or a numeric or ",
           "character vector of the columns to include.")
  )
  expect_equal(
    tryCatch(createSubset(es, subsetName = "s3", parentAssay = 1),
             error = conditionMessage),
    paste0("'parentAssay' must be NULL, to use the default assay, or a ",
           "single character string naming the assay or subset to make ",
           "the subset from.")
  )
  expect_error(es <- createSubset(es, subsetName = "s3", rows = c(64103:64103)),
               "Selected rows not available in the specified assay.")
  expect_error(es <- createSubset(es, subsetName = "s3", rows = c(seq(1:64102), 64102)),
               "Selected rows not available in the specified assay.")
  expect_error(es <- createSubset(es, subsetName = "s3", rows = c(1,3, NA)),
               "Selected rows not available in the specified assay.")
  expect_error(es <- createSubset(es, subsetName = "s3", cols = c(10:12)),
               "Selected columns not available in the specified assay.")
  expect_error(es <- createSubset(es, subsetName = "s3", cols = c(seq(1:8), 8)),
               "Selected columns not available in the specified assay.")
  expect_error(es <- createSubset(es, subsetName = "s3", cols = c(1,3, NA)),
               "Selected columns not available in the specified assay.")
  
})

test_that("Testing ExperimentSubset constructor with subset & createSubset without parentAssay",{
  es <- ExperimentSubset(sce_chcl, subset = list(subsetName = "subset10", rows = c(1:10), cols = c(1:2), parentAssay = "counts"))
  expect_equal(subsetNames(es), "subset10")
  expectSubsetOf(es, "subset10", rows = 1:10, cols = 1:2)
  
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1)
  expect_equal(subsetNames(es), c("subset10", "subset1"))
  expectSubsetOf(es, "subset1", rows = rows1, cols = cols1)
  expectSubsetOf(es, "subset10", rows = 1:10, cols = 1:2)
})

test_that("Testing createSubset and setSubsetAssay",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  #createSubset
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  expectSubsetOf(es, "subset1", rows = rows1, cols = cols1)
  
  subsetRownames(es, "subset1") <- subsetRownames(es, "subset1")
  subsetColnames(es, "subset1") <- subsetColnames(es, "subset1")
  expect_equal(subsetRownames(es, "subset1"), rownames(es)[rows1])
  expect_equal(subsetColnames(es, "subset1"), colnames(es)[cols1])
  
  # the same rows and columns, given by name
  es <- createSubset(es, "subset2",
                                       rows = c("OPCML", "RP11-338K17.8", "RP11-770J1.3", "NPSR1", "AC060835.1", "GABRA2", "PRDM16", "CAB39", "LYRM1", "RP1-178F10.3"),
                                       cols = c("TTAGTTCTCGTAGGAG", "GGCTGGTGTCTCGTTC", "TCGCGTTGTTAAGTAG", "GGATGTTGTAGGCATG", "AATCGGTTCTGATACG", "GTTCTCGGTATGAAAC", "TCCACACTCTTTAGTC", "GCGCAGTAGATGCGAC", "CGGACACTCAAGAAGT", "TGACTTTGTCGCGAAA"),
                                       parentAssay = "counts")
  expectSubsetOf(es, "subset2", rows = rows1, cols = cols1)
  
  # numbers given as names are not row names
  expect_error(
    createSubset(es, "subset3", rows = c("OPCML", "RP11-338K17.8", "50"),
                 cols = 1:2, parentAssay = "counts"),
    regexp = "Selected rows not available in the specified assay."
  )
  
  #setSubsetAssay
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  # the same message as other functions give for an unknown subset
  expect_equal(
    tryCatch(setSubsetAssay(es, "subset4", counts1p, "scaledSubset1"),
             error = conditionMessage),
    paste0("subset4 subset does not exist. 'subsetName' must be the name ",
           "of a subset.")
  )
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  expect_true("subset1/scaledSubset1" %in% subsetAssayNames(es))
  expect_equal(
    getSubsetAssay(es, "subset1", subsetAssayName = "scaledSubset1"), counts1p
  )
  # the subset's own data is unchanged
  expect_equal(getSubsetAssay(es, "subset1"),
               assay(es, "counts")[rows1, cols1])
  # an assay stored in a subset needs a name
  # NULL, NA and 1L failed at dispatch before reaching this check (#97)
  for (badName in list(NA_character_, "", c("a", "b"), NULL, NA, 1L)) {
    expect_error(
      setSubsetAssay(es, "subset1", counts1p, subsetAssayName = badName),
      regexp = "'subsetAssayName' must be a single, non-empty character",
      info = deparse(badName)
    )
  }
})

test_that("Testing rownames, colnames, rowData and colData",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  counts1p <- assay(es, "subset1")
  counts1p[,] <- counts1p[,] + 1
  es <- setSubsetAssay(es, "subset1", counts1p, "scaledSubset1")
  
  expect_equal(subsetRownames(es, "subset1/scaledSubset1"),
               rownames(es)[rows1])
  expect_equal(subsetColnames(es, "subset1/scaledSubset1"),
               colnames(es)[cols1])
  expect_equal(subsetRownames(es, "subset1"), rownames(es)[rows1])
  expect_equal(subsetColnames(es, "subset1"), colnames(es)[cols1])
  expect_error(colData(es, subsetName = "subset1Internal"),
               "subset1Internal subset does not exist.")
  expect_error(rowData(es, subsetName = "subset1Internal"),
               "subset1Internal subset does not exist.")
  rowData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[1])))
  colData(es, subsetName = "subset1") <- DataFrame(col1 = c(seq(subsetDim(es, subsetName = "subset1")[2])))
  expect_equal(ncol(rowData(es, subsetName = "subset1", parentRowData = TRUE)), 2)
  expect_equal(nrow(rowData(es, subsetName = "subset1", parentRowData = TRUE)), subsetDim(es, subsetName = "subset1")[1])
  expect_equal(ncol(colData(es, subsetName = "subset1", parentColData = TRUE)), 12)
  expect_equal(nrow(colData(es, subsetName = "subset1", parentColData = TRUE)),
               subsetDim(es, subsetName = "subset1")[2])
  expect_equal(rowData(es, subsetName = "subset1")$col1,
               seq_along(rows1))
  expect_equal(colData(es, subsetName = "subset1")$col1,
               seq_along(cols1))
  # the subset is square, so check names to tell its rows from its columns
  expect_equal(
    rownames(rowData(es, subsetName = "subset1", parentRowData = TRUE)),
    rownames(es)[rows1]
  )
  expect_equal(
    rownames(colData(es, subsetName = "subset1", parentColData = TRUE)),
    colnames(es)[cols1]
  )
})

test_that("Testing subset helper/supplementary functions",{
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
  expect_output(subsetSummary(es),
                "subset1 +10, 10 +counts +scaledSubset1")
  expect_equal(subsetParent(es, "subset1/scaledSubset1"),
               list("subset1", "counts"))
  expect_output(show(es),
                "subsetAssays\\(2\\): subset1 subset1/scaledSubset1")
})

test_that("Testing assay function",{
  es <- ExperimentSubset(sce_chcl)
  rows1 <- c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200)
  cols1 <- c(20, 21, 40, 45, 90, 99, 100, 123, 166, 299)
  es <- createSubset(es, "subset1", rows = rows1, cols = cols1,
                     parentAssay = "counts")
  subset1Assay <- assay(es, "subset1")
  expect_equal(subset1Assay, assay(es, "counts")[rows1, cols1])
  assay(es, "counts2") <- assay(es, "counts")
  expect_equal(assay(es, "counts2"), assay(es, "counts"))
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
  # a fixture, so its "Storing assay inside subset" message isn't tested
  suppressMessages(
    assay(
      es, "filteredAssay", withDimnames = FALSE,
      subsetAssayName = "filteredAssayNormalized"
    ) <- scrapper::normalizeCounts(
      filteredCounts, size.factors = filteredSizeFactors
    )
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

test_that("Testing reducedDims",{
  skip_if_not_installed("scrapper")
  skip_if_not_installed("scater")
  es <- makeHvgFixture()
  set.seed(20)
  colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), cluster = kmeans(t(assay(es, "hvg1000")), 5)$cluster)
  expect_true("cluster" %in% colnames(colData(es, subsetName = "hvg1000")))

  hvgPca <- scater::calculatePCA(assay(es, "hvg1000"))
  reducedDim(es, "PCA", subsetName = "hvg1000") <- hvgPca
  expect_true("PCA" %in% reducedDimNames(es, subsetName = "hvg1000"))
  pcaList <- list(a = hvgPca, b = hvgPca)
  reducedDims(es, subsetName = "hvg1000") <- pcaList
  expect_equal(reducedDimNames(es, subsetName = "hvg1000"),
               c("a", "b"))
  expect_equal(reducedDim(es, "b", subsetName = "hvg1000"), hvgPca)
  
  reducedDimNames(es, subsetName = "hvg1000") <- c("PCA_1", "PCA_2")
  expect_equal(reducedDimNames(es, subsetName = "hvg1000"),
               c("PCA_1", "PCA_2"))
  expect_equal(reducedDims(es, subsetName = "hvg1000"),
               SimpleList(PCA_1 = hvgPca, PCA_2 = hvgPca))

  expect_error(reducedDimNames(es, subsetName = "hvg123"),
               "hvg123 subset does not exist.")
  
  expect_error(reducedDimNames(es, subsetName = "hvg123") <- c("PCA_1"),
               "hvg123 subset does not exist.")
  
  countsPca <- scater::calculatePCA(assay(es, "counts"))
  reducedDim(es, type = "PCA") <- countsPca
  expect_equal(reducedDim(es, type = "PCA"), countsPca)

  reducedDims(es) <- list(a = countsPca, b = countsPca * 2)
  expect_equal(reducedDims(es),
               SimpleList(a = countsPca, b = countsPca * 2))
  expect_equal(reducedDimNames(es), c("a", "b"))

  reducedDimNames(es) <- c("PCA_Counts_1", "PCA_Counts_2")
  expect_equal(reducedDimNames(es), c("PCA_Counts_1", "PCA_Counts_2"))
  expect_equal(reducedDim(es, "PCA_Counts_2"), countsPca * 2)

})

test_that("Testing metadata",{
  skip_if_not_installed("scrapper")
  es <- makeHvgFixture()
  metadata(es, subsetName = "hvg1000") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(metadata(es, subsetName = "hvg1000")$meta1, "This is testing meta1")
  expect_equal(metadata(es, subsetName = "hvg1000")$meta2, "This is testing meta2")
  
  expect_error(metadata(es, subsetName = "hvg") <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2"),
               "hvg subset does not exist.")
  
  expect_error(metadata(es, subsetName = "hvg"),
               "hvg subset does not exist.")
  
  metadata(es) <- list(meta1 = "This is testing meta1", meta2 = "This is testing meta2")
  expect_equal(metadata(es)$meta1, "This is testing meta1")
  expect_equal(metadata(es)$meta2, "This is testing meta2")
})

test_that("Testing altExperiments",{
  skip_if_not_installed("scrapper")
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
  expect_equal(altExpNames(es, subsetName = "hvg1000"),
               c("alt1", "alt2"))

  expect_error(altExp(es, "alt1", subsetName = "hvg") <- hvgAlt,
               "hvg subset does not exist.")
  
  # without a name, altExp<- replaces the first altExp
  altExp(es, subsetName = "hvg1000") <- makeHvgAlt(2)
  expect_equal(length(altExpNames(es, subsetName = "hvg1000")), 2)
  
  expect_equal(
    assay(altExp(es, subsetName = "hvg1000"), "counts"),
    assay(makeHvgAlt(2), "counts")
  )
  expect_equal(
    assay(altExp(es, "alt2", subsetName = "hvg1000"), "counts"),
    assay(hvgAlt, "counts")
  )
  
  altExp(es) <- simSce
  expect_equal(assay(altExp(es), "counts"), simCounts)

  altExp(es, "a1") <- simSce
  expect_true("a1" %in% altExpNames(es))
  expect_equal(assay(altExp(es, "a1"), "counts"), simCounts)
  
  expect_error(altExpNames(es, subsetName = "na"),
               "na subset does not exist.")
  
  expect_error(altExpNames(es, subsetName = "na") <- c("na2"),
               "na subset does not exist.")
  
  altExpNames(es) <- c("a1", "a2")
  expect_equal(altExpNames(es), c("a1", "a2"))
  
  expect_error(altExp(es, "alt1", subsetName = "hvg10002"),
               "hvg10002 subset does not exist.")
  
  # the subset sees the main object's altExps a1 and a2, then its own
  # alt1 and alt2, and renames only its own (#58)
  altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2", "b1", "b2")
  expect_equal(altExpNames(es, subsetName = "hvg1000"),
               c("a1", "a2", "b1", "b2"))

  expect_error(altExpNames(es, subsetName = "hvg1000x"),
               "hvg1000x subset does not exist.")
  
  altExps(es, subsetName = "hvg1000") <- list(a = hvgAlt, b = makeHvgAlt(3))
  # its own altExps are replaced; the main object's stay visible (#58)
  expect_equal(altExpNames(es, subsetName = "hvg1000"),
               c("a1", "a2", "a", "b"))
  expect_equal(
    assay(altExp(es, "b", subsetName = "hvg1000"), "counts"),
    assay(makeHvgAlt(3), "counts")
  )

  expect_error(altExps(es, subsetName = "hvg") <- list(a = hvgAlt),
               "hvg subset does not exist.")
  
  expect_error(altExps(es, subsetName = "hvg10002"),
               "hvg10002 subset does not exist.")
  
  altExps(es) <- list(a = simSce, b = simSce)
  expect_equal(altExpNames(es), c("a", "b"))
} 
)

test_that("altExp gets and sets the altExp chosen by name or index", {
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

  expect_equal(nrow(altExp(es, "alt2", subsetName = "small")), 3L)
  expect_equal(nrow(altExp(es, 2, subsetName = "small")), 3L)
  expect_equal(nrow(altExp(es, subsetName = "small")), 2L)

  altExp(es, 2, subsetName = "small") <- alt3
  expect_equal(nrow(altExp(es, "alt1", subsetName = "small")), 2L)
  expect_equal(nrow(altExp(es, "alt2", subsetName = "small")), 4L)

  # same calls on the main object, without subsetName
  altExp(es, "a1") <- makeAlt(1, 2, colnames(es))
  altExp(es, "a2") <- makeAlt(2, 3, colnames(es))
  altExp(es, 2) <- makeAlt(3, 4, colnames(es))
  expect_equal(nrow(altExp(es, "a1")), 2L)
  expect_equal(nrow(altExp(es, 2)), 4L)
})

test_that("Accessors on the main object pass extra arguments on", {
  # withDimnames is passed through the generic's ... for these getters
  sce <- sce_chcl
  reducedDim(sce, "PCA") <- matrix(seq_len(ncol(sce) * 2), ncol = 2)
  # stored altExp colnames differ, so withDimnames changes the result
  altCells <- paste0("alt_", seq_len(ncol(sce)))
  counts <- matrix(1, 2, ncol(sce), dimnames = list(NULL, altCells))
  alt <- SingleCellExperiment(list(counts = counts))
  altExps(sce, withDimnames = FALSE) <- list(alt1 = alt)
  es <- ExperimentSubset(sce)

  expect_equal(reducedDim(es, "PCA", withDimnames = FALSE),
               reducedDim(sce, "PCA", withDimnames = FALSE))
  expect_equal(reducedDims(es, withDimnames = FALSE),
               reducedDims(sce, withDimnames = FALSE))
  expect_equal(altExps(es, withDimnames = FALSE),
               altExps(sce, withDimnames = FALSE))
  expect_equal(colnames(altExps(es, withDimnames = FALSE)[[1]]),
               altCells)
})

test_that("Subsets of subsets use the right rows and columns", {
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

  expect_equal(rownames(assay(es, "s2")), rownames(sce)[101:110])
  expect_equal(subsetRownames(es, "s2"), rownames(sce)[101:110])
  expect_equal(subsetColnames(es, "s2"), colnames(sce)[501:510])
  expect_equal(subsetRownames(es, "s3"), rownames(sce)[102:103])
  expect_equal(subsetColnames(es, "s3"), colnames(sce)[504:505])
  expect_equal(subsetRownames(es, "s4"), rownames(sce)[105:106])
  expect_equal(subsetColnames(es, "s4"), colnames(sce)[502:503])

  s2RowData <- rowData(es, subsetName = "s2", parentRowData = TRUE)
  s2ColData <- colData(es, subsetName = "s2", parentColData = TRUE)
  expect_equal(s2RowData$geneIndex, 101:110)
  expect_equal(s2ColData$cellIndex, 501:510)
  s3ColData <- colData(es, subsetName = "s3", parentColData = TRUE)
  expect_equal(s3ColData$cellIndex, 504:505)
})

test_that("Subset rowData and colData keep the subset's names", {
  sce <- sce_chcl
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  genes <- rownames(sce)[1:5]
  cells <- colnames(sce)[1:8]

  # getters return the subset's row and column names
  expect_equal(rownames(rowData(es, subsetName = "s1")), genes)
  expect_equal(rownames(colData(es, subsetName = "s1")), cells)
  s1ColData <- subsetColData(es, "s1", parentColData = TRUE)
  expect_equal(rownames(s1ColData), cells)
  expect_equal(s1ColData$cellIndex, 1:8)

  # get, modify and set back keeps the cell names
  s1ColData <- subsetColData(es, "s1", parentColData = FALSE)
  s1ColData$score <- seq_len(8)
  subsetColData(es, "s1") <- s1ColData
  expect_equal(subsetColnames(es, "s1"), cells)
  expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  # so a reducedDim with the cell names can be stored without a warning
  umap <- matrix(0, 8, 2, dimnames = list(cells, c("U1", "U2")))
  expect_warning(reducedDim(es, "umap", subsetName = "s1") <- umap,
                 NA)

  # a value without rownames takes the subset's names
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:8])
  expect_equal(subsetColnames(es, "s1"), cells)
  rowData(es, subsetName = "s1") <- DataFrame(rank = 5:1)
  expect_equal(subsetRownames(es, "s1"), genes)

  # colData with other rownames, or of the wrong size, is an error
  otherCells <- DataFrame(group = letters[1:8], row.names = LETTERS[1:8])
  expect_error(subsetColData(es, "s1") <- otherCells,
               regexp = "rownames of 'value'")
  expect_error(colData(es, subsetName = "s1") <- otherCells,
               regexp = "rownames of 'value'")
  expect_error(subsetColData(es, "s1") <- DataFrame(a = 1:5),
               regexp = "nrow")
  expect_equal(subsetColnames(es, "s1"), cells)
  # rowData with other rownames is an error too, so values meant for other
  # rows aren't stored by position (#100)
  rowData(es, subsetName = "s1") <- DataFrame(rank = 11:15, row.names = genes)
  expect_equal(rowData(es, subsetName = "s1")$rank, 11:15)
  rowMsg <- paste0("The rownames of 'value' must be NULL or the row names ",
                   "of subset 's1', in the same order.")
  otherGenes <- DataFrame(rank = 1:5, row.names = LETTERS[1:5])
  expect_error(subsetRowData(es, "s1") <- otherGenes,
               regexp = rowMsg, fixed = TRUE)
  expect_error(rowData(es, subsetName = "s1") <- otherGenes,
               regexp = rowMsg, fixed = TRUE)
  reordered <- DataFrame(rank = 1:5, row.names = rev(genes))
  expect_error(subsetRowData(es, "s1") <- reordered,
               regexp = rowMsg, fixed = TRUE)
  otherFrame <- data.frame(rank = 1:5, row.names = paste0("x", 1:5))
  expect_error(subsetRowData(es, "s1") <- otherFrame,
               regexp = rowMsg, fixed = TRUE)
  expect_equal(subsetRownames(es, "s1"), genes)
  expect_equal(rowData(es, subsetName = "s1")$rank, 11:15)
  # a value of the wrong size gets the size error, even with rownames
  tooFewGenes <- DataFrame(rank = 1:3, row.names = genes[1:3])
  tooFewCells <- DataFrame(group = letters[1:5], row.names = cells[1:5])
  expect_error(subsetRowData(es, "s1") <- tooFewGenes,
               regexp = "3 rows in value to replace 5")
  expect_error(subsetColData(es, "s1") <- tooFewCells,
               regexp = "nrow")
  # the colData message names the subset's column names
  expect_error(subsetColData(es, "s1") <- otherCells,
               regexp = paste0("The rownames of 'value' must be ",
                               "NULL or the column names of ",
                               "subset 's1', in the same order."),
               fixed = TRUE)
  # a subset without row names has none to check against
  unnamed <- ExperimentSubset(list(counts = matrix(1:20, nrow = 5)))
  unnamed <- createSubset(unnamed, "s1", rows = 1:3, cols = 1:2)
  rowData(unnamed, subsetName = "s1") <- otherGenes[1:3, , drop = FALSE]
  expect_equal(rowData(unnamed, subsetName = "s1")$rank, 1:3)
  expect_null(subsetRownames(unnamed, "s1"))

  # with the parent's data, the result has the subset's renamed names
  subsetColnames(es, "s1") <- paste0("cell", 1:8)
  subsetRownames(es, "s1") <- paste0("gene", 1:5)
  s1ColData <- subsetColData(es, "s1", parentColData = TRUE)
  s1RowData <- subsetRowData(es, "s1", parentRowData = TRUE)
  expect_equal(rownames(s1ColData), paste0("cell", 1:8))
  expect_equal(rownames(s1RowData), paste0("gene", 1:5))
})

test_that("Subsets of subsets use the right spatial coordinates", {
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

  expect_equal(subsetColnames(es, "s2"), spots[14:18])
  expect_equal(unname(subsetSpatialCoords(es, "s2")[, "x"]),
               as.numeric(14:18))
  expect_equal(colData(es, subsetName = "s2")$sample_id,
               c("A", "A", "B", "B", "B"))
})

test_that("Subsets of subsets use the right rowLinks", {
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
  expect_equal(links$nodeLab, rowLinks(tse)$nodeLab[4:6])
})

test_that("Every Subset* class is an ExperimentSubset (#5)", {
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
    expect_equal(class(es)[1], expectedClasses[[name]],
                 info = name)
    expect_true(is(es, "ExperimentSubset"), info = name)
  }

  # the issue's case: a SubsetSpatialExperiment isn't a
  # SubsetSingleCellExperiment, but it is an ExperimentSubset
  expect_false(is(ExperimentSubset(spe),
                  "SubsetSingleCellExperiment"))
  expect_false(is(sce, "ExperimentSubset"))
  expect_false(is(se, "ExperimentSubset"))
  expect_false(is(spe, "ExperimentSubset"))

  # a new Subset* class must be added to the union in allClasses.R
  packageClasses <- getClasses(asNamespace("ExperimentSubset"))
  subsetClasses <- grep("^Subset", packageClasses, value = TRUE)
  expect_length(subsetClasses, 5)
  for (subsetClass in subsetClasses) {
    expect_true(extends(subsetClass, "ExperimentSubset"),
                info = subsetClass)
  }
})

# Runs the methods every Subset*Experiment class shares on es, whose main
# assay "counts" has row and column names and at least 5 rows and 6
# columns. Returns es with a subset "s1" and a stored assay "doubled".
checkSharedMethods <- function(es) {
  counts <- as.matrix(assay(es, "counts"))
  es <- createSubset(es, "s1", rows = 2:5, cols = 3:6)
  expect_equal(subsetNames(es), "s1")
  expect_equal(subsetCount(es), 1L)
  expect_equal(subsetDim(es, "s1"), c(4L, 4L))
  expect_equal(subsetRownames(es, "s1"), rownames(counts)[2:5])
  expect_equal(subsetColnames(es, "s1"), colnames(counts)[3:6])
  expect_equal(as.matrix(assay(es, "s1")), counts[2:5, 3:6])

  doubled <- counts[2:5, 3:6] * 2
  es <- setSubsetAssay(es, "s1", doubled, "doubled")
  storedDoubled <- getSubsetAssay(es, "s1", subsetAssayName = "doubled")
  expect_equal(as.matrix(storedDoubled), doubled)
  expect_equal(subsetAssayNames(es), c("s1", "s1/doubled"))
  expect_equal(subsetAssayCount(es), 2L)
  expect_equal(subsetParent(es, "s1/doubled"), list("s1", "counts"))
  expect_output(show(es), "subsets\\(1\\): s1")
  expect_output(subsetSummary(es), "doubled")

  metadata(es, subsetName = "s1") <- list(note = "kept")
  expect_equal(metadata(es, subsetName = "s1")$note, "kept")

  genes <- subsetRownames(es, "s1")
  cells <- subsetColnames(es, "s1")
  rowData(es, subsetName = "s1") <- DataFrame(score = 1:4, row.names = genes)
  expect_equal(rowData(es, subsetName = "s1")$score, 1:4)
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4],
                                              row.names = cells)
  expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  subsetRowData(es, "s1") <- DataFrame(rank = 4:1, row.names = genes)
  expect_equal(subsetRowData(es, "s1", FALSE)$rank, 4:1)
  subsetColData(es, "s1") <- DataFrame(batch = 1:4, row.names = cells)
  expect_equal(subsetColData(es, "s1", FALSE)$batch, 1:4)
  parentCols <- subsetColData(es, "s1", parentColData = TRUE)
  expect_equal(rownames(parentCols), colnames(counts)[3:6])

  # reducedDims and altExps, which these classes get from SCE
  pca <- matrix(as.numeric(1:8), nrow = 4, dimnames = list(cells, NULL))
  reducedDim(es, "PCA", subsetName = "s1") <- pca
  expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca)
  expect_equal(reducedDimNames(es, subsetName = "s1"), "PCA")
  expect_equal(names(reducedDims(es, subsetName = "s1")), "PCA")
  altCounts <- matrix(1, 2, 4, dimnames = list(NULL, cells))
  alt <- SingleCellExperiment(list(counts = altCounts))
  altExp(es, "alt", subsetName = "s1") <- alt
  expect_equal(altExpNames(es, subsetName = "s1"), "alt")
  expect_equal(dim(altExp(es, "alt", subsetName = "s1")), c(2L, 4L))
  expect_equal(names(altExps(es, subsetName = "s1")), "alt")

  # access by position or default, and the remaining setters
  expect_equal(reducedDim(es, 1, subsetName = "s1"), pca)
  reducedDim(es, 1, subsetName = "s1") <- pca * 2
  expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca * 2)
  expect_equal(dim(altExp(es, 1, subsetName = "s1")), c(2L, 4L))
  expect_equal(dim(altExp(es, subsetName = "s1")), c(2L, 4L))
  altExp(es, subsetName = "s1") <- alt
  altExp(es, 1, subsetName = "s1") <- alt
  altExps(es, subsetName = "s1") <- list(first = alt, second = alt)
  expect_equal(altExpNames(es, subsetName = "s1"),
               c("first", "second"))
  altExpNames(es, subsetName = "s1") <- c("x1", "x2")
  expect_equal(altExpNames(es, subsetName = "s1"), c("x1", "x2"))
  tripled <- counts[2:5, 3:6] * 3
  suppressMessages(assay(es, "s1", subsetAssayName = "tripled") <- tripled)
  expect_equal(as.matrix(assay(es, "s1/tripled")), tripled)

  # renaming comes last, since later lookups would see the new names
  subsetRownames(es, "s1") <- paste0("g", 1:4)
  subsetColnames(es, "s1") <- paste0("c", 1:4)
  expect_equal(subsetRownames(es, "s1"), paste0("g", 1:4))
  expect_equal(subsetColnames(es, "s1"), paste0("c", 1:4))
  es
}

test_that("TreeSummarizedExperiment subsets support all methods", {
  data("tinyTree", package = "TreeSummarizedExperiment",
       envir = environment())
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(tinyTree$tip.label, paste0("cell", 1:6)))
  tse <- TreeSummarizedExperiment(assays = list(counts = counts),
                                  rowTree = tinyTree)
  es <- ExperimentSubset(tse)
  expect_s4_class(es, "SubsetTreeSummarizedExperiment")
  es <- checkSharedMethods(es)

  links <- subsetRowLinks(es, "s1", parentRowLinkData = TRUE)
  expect_equal(links$nodeLab, rowLinks(tse)$nodeLab[2:5])
  # the subset itself stores no links
  expect_null(subsetRowLinks(es, "s1", parentRowLinkData = FALSE))
  expect_null(subsetColLinks(es, "s1", parentColLinkData = FALSE))
  # parentColLinkData = TRUE errors without a column tree (#44)
})

test_that("SpatialExperiment subsets support all methods", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  coords <- cbind(x = as.numeric(1:6), y = as.numeric(1:6) * 10)
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = coords)
  es <- ExperimentSubset(spe)
  expect_s4_class(es, "SubsetSpatialExperiment")
  es <- checkSharedMethods(es)

  spotCoords <- spatialCoords(es, subsetName = "s1")
  expect_equal(unname(spotCoords[, "y"]), c(30, 40, 50, 60))
  expect_equal(subsetSpatialCoords(es, "s1"), spotCoords)
  # the main object's coordinates are unchanged
  expect_equal(spatialCoords(es), spatialCoords(spe))
})

# A SpatialExperiment with six spots in two samples, A and B, and an
# in-memory image for each sample
makeImageSpe <- function() {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  image <- SpatialImage(as.raster(matrix(c(0, 0.5, 1, 0.5), nrow = 2)))
  images <- S4Vectors::DataFrame(sample_id = c("A", "B"),
                                 image_id = "lowres",
                                 data = I(list(image, image)),
                                 scaleFactor = c(1, 0.5))
  SpatialExperiment(assays = list(counts = counts),
                    sample_id = rep(c("A", "B"), each = 3),
                    spatialCoords = cbind(x = as.numeric(1:6),
                                          y = as.numeric(1:6) * 10),
                    imgData = images)
}

test_that("subsetImgData returns the images of a subset's samples", {
  es <- ExperimentSubset(makeImageSpe())
  es <- createSubset(es, "sampleB", cols = 4:6)
  es <- createSubset(es, "both", cols = 3:4)
  expect_equal(subsetImgData(es, "sampleB"), imgData(es)[2, ])
  expect_equal(subsetImgData(es, "both"), imgData(es))
  # a subset of a subset
  es <- createSubset(es, "childA", cols = 1, parentAssay = "both")
  expect_equal(subsetImgData(es, "childA"), imgData(es)[1, ])
  # the same images as getSubsetExperiment() carries
  expect_equal(subsetImgData(es, "sampleB"),
               imgData(getSubsetExperiment(es, "sampleB")))
  expect_error(subsetImgData(es, "nope"),
               regexp = "nope subset does not exist")
  # an object without images gives none
  noImages <- makeImageSpe()
  imgData(noImages) <- NULL
  es <- createSubset(ExperimentSubset(noImages), "s1", cols = 1:2)
  expect_equal(nrow(subsetImgData(es, "s1")), 0L)
})

test_that("A subset can store its own spatialCoords", {
  es <- ExperimentSubset(makeImageSpe())
  es <- createSubset(es, "s1", cols = 3:6)
  mainCoords <- spatialCoords(es)
  shifted <- subsetSpatialCoords(es, "s1") + 100
  rownames(shifted) <- NULL

  # stored on the subset; the main object is unchanged
  subsetSpatialCoords(es, "s1") <- shifted
  expect_equal(unname(subsetSpatialCoords(es, "s1")), unname(shifted))
  expect_equal(spatialCoords(es, subsetName = "s1"),
               subsetSpatialCoords(es, "s1"))
  expect_equal(colnames(subsetSpatialCoords(es, "s1")), c("x", "y"))
  expect_equal(spatialCoords(es), mainCoords)

  # a subset made from it sees its coordinates, cut to its spots, also
  # when made from one of its assays
  es <- createSubset(es, "child", cols = 2:3, parentAssay = "s1")
  expect_equal(unname(subsetSpatialCoords(es, "child")),
               unname(shifted[2:3, ]))
  es <- setSubsetAssay(es, "s1", assay(es, "s1") * 2, "doubled")
  es <- createSubset(es, "fromAssay", cols = 4, parentAssay = "s1/doubled")
  expect_equal(unname(subsetSpatialCoords(es, "fromAssay")),
               unname(shifted[4, , drop = FALSE]))

  # the subset's own coordinates follow [ on the main object
  cut <- es[, 2:6]
  expect_equal(unname(subsetSpatialCoords(cut, "s1")), unname(shifted))
  withoutSpot3 <- es[, c(1, 2, 4, 5, 6)]
  expect_equal(unname(subsetSpatialCoords(withoutSpot3, "s1")),
               unname(shifted[2:4, ]))
  # and cbind(), where the first object's subsets keep their spots
  other <- ExperimentSubset(makeImageSpe())
  colnames(other) <- paste0("other", 1:6)
  # SpatialExperiment notes the shared sample_ids
  combined <- suppressMessages(cbind(es, other))
  expect_equal(unname(subsetSpatialCoords(combined, "s1")), unname(shifted))

  # renamed spots rename the coordinates
  renamed <- es
  subsetColnames(renamed, "s1") <- paste0("new", 1:4)
  expect_equal(rownames(subsetSpatialCoords(renamed, "s1")),
               paste0("new", 1:4))
  expect_equal(rownames(spatialCoords(getSubsetExperiment(renamed, "s1"))),
               paste0("new", 1:4))

  # getSubsetExperiment() carries them
  standalone <- getSubsetExperiment(es, "s1")
  expect_equal(unname(spatialCoords(standalone)), unname(shifted))

  # NULL goes back to the main object's coordinates
  subsetSpatialCoords(es, "s1") <- NULL
  expect_equal(subsetSpatialCoords(es, "s1"),
               mainCoords[3:6, , drop = FALSE])
  expect_equal(unname(subsetSpatialCoords(es, "child")),
               unname(mainCoords[4:5, ]))

  # rownames must be NULL or the subset's spot names, in order
  named <- shifted
  rownames(named) <- paste0("spot", 3:6)
  subsetSpatialCoords(es, "s1") <- named
  expect_equal(rownames(subsetSpatialCoords(es, "s1")),
               paste0("spot", 3:6))
  rownames(named) <- paste0("spot", c(4, 3, 5, 6))
  expect_error(subsetSpatialCoords(es, "s1") <- named,
               regexp = paste0("The rownames of 'value' must be NULL or ",
                               "the column names of subset 's1', in the ",
                               "same order."), fixed = TRUE)
})

test_that("subsetSpatialCoords<- checks value", {
  es <- ExperimentSubset(makeImageSpe())
  es <- createSubset(es, "s1", cols = 3:6)
  expect_error(subsetSpatialCoords(es, "s1") <- data.frame(x = 1:4),
               regexp = "'value' must be a numeric matrix or NULL.",
               fixed = TRUE)
  expect_error(subsetSpatialCoords(es, "s1") <- matrix("a", 4, 2),
               regexp = "'value' must be a numeric matrix or NULL.",
               fixed = TRUE)
  expect_error(subsetSpatialCoords(es, "s1") <- matrix(1, 4, 0),
               regexp = paste0("'value' must have one or more columns; use ",
                               "NULL to remove the subset's own ",
                               "coordinates."), fixed = TRUE)
  expect_error(subsetSpatialCoords(es, "s1") <- matrix(1, 3, 2),
               regexp = paste0("'value' must have one row for each of the ",
                               "4 spots in subset 's1'."), fixed = TRUE)
  expect_error(subsetSpatialCoords(es, "nope") <- matrix(1, 4, 2),
               regexp = "nope subset does not exist")
})

test_that("setSubsetExperiment stores changed spatialCoords only", {
  es <- ExperimentSubset(makeImageSpe())
  es <- createSubset(es, "s1", cols = 3:6)
  mainCoords <- spatialCoords(es)

  # unchanged coordinates are not stored: the subset still follows the
  # main object's
  standalone <- getSubsetExperiment(es, "s1")
  expect_no_message(es <- setSubsetExperiment(es, "s1", standalone))
  expect_equal(subsetSpatialCoords(es, "s1"),
               mainCoords[3:6, , drop = FALSE])
  doubled <- es
  spatialCoords(doubled) <- spatialCoords(doubled) * 2
  expect_equal(unname(subsetSpatialCoords(doubled, "s1")),
               unname(mainCoords[3:6, ] * 2))

  # changed coordinates are stored as the subset's own
  spatialCoords(standalone) <- spatialCoords(standalone) + 1
  expect_message(es <- setSubsetExperiment(es, "s1", standalone),
                 regexp = "Subset 's1' now stores its own spatialCoords.",
                 fixed = TRUE)
  expect_equal(unname(subsetSpatialCoords(es, "s1")),
               unname(mainCoords[3:6, ] + 1))
  expect_equal(spatialCoords(es), mainCoords)

  # also when the main object's coordinates have other rownames than its
  # spot names
  es <- createSubset(es, "s2", cols = 1:2)
  namedCoords <- mainCoords
  rownames(namedCoords) <- paste0("other", 1:6)
  spatialCoords(es) <- namedCoords
  standalone <- getSubsetExperiment(es, "s2")
  spatialCoords(standalone) <- spatialCoords(standalone) * 3
  expect_message(es <- setSubsetExperiment(es, "s2", standalone),
                 regexp = "Subset 's2' now stores its own spatialCoords.",
                 fixed = TRUE)
  expect_equal(unname(subsetSpatialCoords(es, "s2")),
               unname(mainCoords[1:2, ] * 3))
  expect_equal(rownames(subsetSpatialCoords(es, "s2")), c("spot1", "spot2"))
})

test_that("Subset spatialData is deprecated, its setter defunct", {
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
  expect_warning(
    s1Data <- suppressMessages(
      subsetSpatialData(es, "s1", parentSpatialData = TRUE)
    ),
    "deprecated"
  )
  expect_equal(s1Data$region, letters[2:5])
  expect_warning(
    s1Data <- suppressMessages(
      spatialData(es, subsetName = "s1", parentSpatialData = TRUE)
    ),
    "deprecated"
  )
  expect_equal(s1Data$region, letters[2:5])
  # the main object's spatialData is SpatialExperiment's, with no warning
  expect_warning(suppressMessages(spatialData(es)), NA)
  # one warning, without SpatialExperiment's own notes as well
  expect_message(
    suppressWarnings(subsetSpatialData(es, "s1", parentSpatialData = TRUE)),
    NA
  )
  # a missing subset is an error, before any deprecation warning
  expect_warning(
    expect_error(subsetSpatialData(es, "s9"), "does not exist"),
    NA
  )

  # the setter fails clearly, for both kinds of data frame
  expect_error(
    subsetSpatialData(es, "s1") <- data.frame(region = letters[1:4]),
    "defunct"
  )
  expect_error(
    subsetSpatialData(es, "s1") <- DataFrame(region = letters[1:4]),
    "defunct"
  )
  # and for a subset given by number
  expect_error(
    subsetSpatialData(es, 1) <- data.frame(region = letters[1:4]),
    "defunct"
  )
})

test_that("reducedDim with subsetName works without a type", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  pca <- matrix(1, nrow = 4, ncol = 2)
  reducedDim(es, "PCA", subsetName = "s1") <- pca
  expect_equal(dim(reducedDim(es, subsetName = "s1")), c(4L, 2L))
  # without a type, the setter replaces the first reducedDim
  reducedDim(es, subsetName = "s1") <- pca * 2
  expect_equal(unname(reducedDim(es, "PCA", subsetName = "s1")),
               pca * 2)
})

test_that("reducedDim on the main object works without a type", {
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
    expect_s4_class(es, paste0("Subset", class(objects[[name]])))
    reducedDim(es, "PCA") <- pca
    reducedDim(es, "UMAP") <- pca * 10
    # without a type, the getter returns the first reducedDim
    expect_equal(unname(reducedDim(es)), pca, info = name)
    # without a type, the setter replaces the first reducedDim
    reducedDim(es) <- pca * 2
    expect_equal(reducedDimNames(es), c("PCA", "UMAP"),
                 info = name)
    expect_equal(unname(reducedDim(es, "PCA")), pca * 2,
                 info = name)
    expect_equal(unname(reducedDim(es, "UMAP")), pca * 10,
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

test_that("reducedDims<- is SingleCellExperiment's generic", {
  # (#26) a package generic of the same name masked SCE's
  expect_equal(attr(ExperimentSubset::`reducedDims<-`, "package"),
               "SingleCellExperiment")
})

test_that("reducedDims<- takes withDimnames", {
  # (#26) the package generic had no withDimnames and no ...
  objects <- makeReducedDimObjects()
  pca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  subsetPca <- pca[2:4, ]
  for (name in names(objects)) {
    plain <- objects[[name]]
    expect_no_error(
      reducedDims(plain, withDimnames = FALSE) <- list(PCA = pca)
    )
    expect_equal(unname(reducedDim(plain, "PCA")), pca,
                 info = name)

    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    expect_no_error(
      reducedDims(es, withDimnames = FALSE) <- list(PCA = pca)
    )
    expect_equal(unname(reducedDim(es, "PCA")), pca, info = name)
    expect_no_error(
      reducedDims(es, withDimnames = FALSE, subsetName = "s1") <-
        list(UMAP = subsetPca)
    )
    expect_equal(reducedDimNames(es, subsetName = "s1"), "UMAP",
                 info = name)
    expect_equal(
      unname(reducedDim(es, "UMAP", subsetName = "s1")), subsetPca,
      info = name
    )
    # the main object's reducedDims are unchanged
    expect_equal(reducedDimNames(es), "PCA", info = name)
  }
})

test_that("reducedDims<- needs subsetName by name", {
  # (#26) the old generic took subsetName second; that is now withDimnames
  objects <- makeReducedDimObjects()
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    expect_error(
      reducedDims(es, "s1") <- list(PCA = subsetPca),
      "name the subset with subsetName =", info = name
    )
    expect_equal(reducedDimNames(es, subsetName = "s1"),
                 character(0), info = name)
  }
})

test_that("reducedDims<- on a plain SCE ignores subsetName", {
  # (#26) as for altExps<-: SingleCellExperiment's method drops it
  sce <- makeReducedDimObjects()$sce
  pca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  reducedDims(sce, subsetName = "s1") <- list(PCA = pca)
  expect_equal(reducedDimNames(sce), "PCA")
  expect_equal(unname(reducedDim(sce, "PCA")), pca)
})

test_that("reducedDimNames() takes only subsetName, by name", {
  # (#82) other arguments, or a subset name by position, were dropped
  objects <- makeReducedDimObjects()
  mainPca <- matrix(as.numeric(seq_len(12)), nrow = 6, ncol = 2)
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    reducedDim(es, "MAIN") <- mainPca
    reducedDim(es, "PCA", subsetName = "s1") <- subsetPca
    expect_error(reducedDimNames(es, "s1"),
                 "name the subset with subsetName =", info = name)
    expect_error(reducedDimNames(es, foo = 1),
                 "takes only 'x' and 'subsetName'", info = name)
    expect_error(reducedDimNames(es, subsetName = "s1", foo = 1),
                 "takes only 'x' and 'subsetName'", info = name)
    expect_error(
      reducedDimNames(es, subsetName = "s1", subsetName = "s1"),
      "takes one 'subsetName'", info = name
    )
    # an empty name can't be given as subsetName =
    expect_error(reducedDimNames(es, character(0)),
                 "takes only 'x' and 'subsetName'", info = name)
    # calls without other arguments are unchanged
    expect_equal(reducedDimNames(es), "MAIN", info = name)
    expect_equal(reducedDimNames(es, subsetName = "s1"), "PCA",
                 info = name)
  }
  # a plain SCE has no subsets, so say that rather than how to name one
  sce <- objects$sce
  expect_error(reducedDimNames(sce, "s1"), "which has no subsets")
  expect_equal(reducedDimNames(sce), character(0))
})

test_that("subsetName on a plain SCE says it has no subsets", {
  # (#82) .isSubset() failed in subsetNames() dispatch for these objects
  sce <- makeReducedDimObjects()$sce
  reducedDim(sce, "PCA") <- matrix(as.numeric(seq_len(12)), nrow = 6)
  expect_error(reducedDimNames(sce, subsetName = "s1"),
               "which has no subsets")
  expect_error(reducedDimNames(sce, subsetName = "s1") <- "X",
               "which has no subsets")
  expect_error(reducedDimNames(sce, "s1") <- "X",
               "which has no subsets")
  # without subsetName, the main object is unchanged
  expect_equal(reducedDimNames(sce), "PCA")
})

test_that("reducedDims and their names check the subset name", {
  # (#26) each class's methods, not only the last file in Collate
  objects <- makeReducedDimObjects()
  subsetPca <- matrix(as.numeric(seq_len(6)), nrow = 3, ncol = 2)
  for (name in names(objects)) {
    es <- ExperimentSubset(objects[[name]])
    es <- createSubset(es, "s1", rows = 1:5, cols = 2:4)
    reducedDims(es, subsetName = "s1") <- list(PCA = subsetPca)
    reducedDimNames(es, subsetName = "s1") <- "renamed"
    expect_equal(reducedDimNames(es, subsetName = "s1"), "renamed",
                 info = name)
    expect_equal(reducedDimNames(es), character(0), info = name)
    expect_error(
      reducedDims(es, subsetName = "nope") <- list(PCA = subsetPca),
      "nope subset does not exist", info = name
    )
    expect_error(
      reducedDimNames(es, subsetName = "nope"),
      "nope subset does not exist", info = name
    )
    expect_error(
      reducedDimNames(es, subsetName = "nope") <- "renamed",
      "nope subset does not exist", info = name
    )
  }
})

test_that("Links and spatialData of a subset keep its names", {
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
  expect_equal(rownames(links), tinyTree$tip.label[2:5])
  expect_null(subsetColLinks(es, "s1", parentColLinkData = TRUE))

  spots <- paste0("spot", 1:6)
  spe <- SpatialExperiment(assays = list(counts = counts[, 1:6]),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  colnames(spe) <- spots
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = 2:4)
  # subsetSpatialData() is deprecated (#34)
  expect_warning(
    s1Data <- subsetSpatialData(es, "s1", parentSpatialData = FALSE),
    "deprecated"
  )
  expect_equal(rownames(s1Data), spots[2:4])
})

test_that("withDimnames, withColData and use.names reach SCE", {
  # names set on the object, not on the stored matrix, so the arguments
  # change what comes back
  sce <- SingleCellExperiment(list(counts = matrix(1, 4, 3)))
  rownames(sce) <- paste0("gene", 1:4)
  colnames(sce) <- paste0("cell", 1:3)
  small <- ExperimentSubset(sce)
  expect_null(rownames(assay(small, "counts", withDimnames = FALSE)))
  # SingleCellExperiment keeps row names here; SummarizedExperiment drops them
  se <- SummarizedExperiment(list(counts = matrix(1, 4, 3)))
  rownames(se) <- paste0("gene", 1:4)
  seRowData <- rowData(ExperimentSubset(se), use.names = FALSE)
  expect_null(rownames(seRowData))

  es <- ExperimentSubset(sce_chcl)
  counts <- assay(es, "counts")
  # values with other dimnames are accepted when withDimnames = FALSE
  renamed <- counts
  rownames(renamed) <- paste0("other", seq_len(nrow(renamed)))
  expect_error(assay(es, "renamed", withDimnames = FALSE) <- renamed,
               NA)
  otherCells <- paste0("other", seq_len(ncol(es)))
  pca <- matrix(1, nrow = ncol(es), ncol = 2,
                dimnames = list(otherCells, NULL))
  expect_warning(reducedDim(es, "PCA", withDimnames = FALSE) <- pca,
                 NA)
  altCounts <- matrix(1, nrow = 2, ncol = ncol(es),
                      dimnames = list(NULL, otherCells))
  alt <- SingleCellExperiment(list(counts = altCounts))
  expect_warning(altExp(es, "alt", withDimnames = FALSE) <- alt, NA)
  altList <- list(alt = alt)
  expect_warning(altExps(es, withDimnames = FALSE) <- altList, NA)

  # and on a subset
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  subsetPca <- pca[1:4, ]
  expect_warning(reducedDim(es, "PCA", subsetName = "s1",
                            withDimnames = FALSE) <- subsetPca,
                 NA)
  storedPca <- reducedDim(es, "PCA", subsetName = "s1", withDimnames = FALSE)
  expect_equal(rownames(storedPca), otherCells[1:4])
})

test_that("Writing back parent rowData or colData stores no copies", {
  sce <- sce_chcl
  colData(sce)$batch <- rep(c("A", "B"), length.out = ncol(sce))
  rowData(sce)$symbol <- rownames(sce)
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)

  s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE)
  s1ColData$score <- 1:4
  expect_message(colData(es, subsetName = "s1") <- s1ColData,
                 regexp = "batch")
  expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  withParent <- colData(es, subsetName = "s1", parentColData = TRUE)
  expect_equal(sum(colnames(withParent) == "batch"), 1L)

  # a column the user adds with a parent column's name is stored
  own <- colData(es, subsetName = "s1")
  own$batch <- letters[1:4]
  expect_message(colData(es, subsetName = "s1") <- own, NA)
  expect_equal(colData(es, subsetName = "s1")$batch, letters[1:4])

  s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE)
  s1RowData$rank <- 5:1
  expect_message(rowData(es, subsetName = "s1") <- s1RowData,
                 regexp = "symbol")
  expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")
})

test_that("A subset's own column hides the main object's (#57)", {
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
  expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score.*colData\\(x\\)"
  )
  expect_equal(colnames(s1ColData), c("batch", "score"))
  expect_equal(s1ColData$score, 1:8)
  expect_equal(s1ColData$batch, colData(sce)$batch[1:8])
  expect_message(
    s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE),
    regexp = "rank.*rowData\\(x\\)"
  )
  expect_equal(colnames(s1RowData), c("symbol", "rank"))
  expect_equal(s1RowData$rank, 5:1)

  # writing back keeps the subset's values and leaves out the main object's
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  expect_equal(colData(es)$score, rep(0, ncol(sce)))
  suppressMessages(rowData(es, subsetName = "s1") <- s1RowData)
  expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")
  expect_equal(rowData(es, subsetName = "s1")$rank, 5:1)
  expect_equal(rowData(es)$rank, rep(0, nrow(sce)))

  # the extracted object has the same columns, and stores back the same way
  s1 <- suppressMessages(getSubsetExperiment(es, "s1"))
  expect_equal(colnames(colData(s1)), c("batch", "score"))
  expect_equal(s1$score, 1:8)
  expect_equal(colnames(rowData(s1)), c("symbol", "rank"))
  expect_equal(rowData(s1)$rank, 5:1)
  es <- suppressMessages(setSubsetExperiment(es, "s1", s1))
  expect_equal(colnames(colData(es, subsetName = "s1")), "score")
  expect_equal(colData(es, subsetName = "s1")$score, 1:8)
  expect_equal(colnames(rowData(es, subsetName = "s1")), "rank")

  # a subset of a subset: its own column hides the main object's too
  es <- createSubset(es, "s2", cols = 1:4, parentAssay = "s1")
  colData(es, subsetName = "s2") <- DataFrame(score = 11:14)
  s2ColData <- suppressMessages(
    colData(es, subsetName = "s2", parentColData = TRUE)
  )
  expect_equal(colnames(s2ColData), c("batch", "score"))
  expect_equal(s2ColData$score, 11:14)

  # the same values: one column, and no message
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0, 8))
  expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  expect_equal(colnames(s1ColData), c("batch", "score"))
  rowData(es, subsetName = "s1") <- DataFrame(rank = rep(0, 5))
  expect_message(
    s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE),
    NA
  )
  expect_equal(colnames(s1RowData), c("symbol", "rank"))

  # equal values of another type count as different: integer, not double
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0L, 8))
  expect_message(
    colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score"
  )
})

test_that("A main object's repeated column name is checked (#57)", {
  sce <- simSce
  sameScores <- DataFrame(score = rep(0, ncol(sce)))
  otherScores <- DataFrame(score = rep(9, ncol(sce)))
  colData(sce) <- cbind(sameScores, otherScores)
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", cols = 1:8)
  colData(es, subsetName = "s1") <- DataFrame(score = rep(0, 8))
  # the second copy differs from the subset's column
  expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    regexp = "score"
  )
  expect_equal(colnames(s1ColData), "score")
})

test_that("SpatialExperiment subsets have one sample_id (#57)", {
  counts <- matrix(as.numeric(seq_len(60)), nrow = 10,
                   dimnames = list(paste0("gene", 1:10),
                                   paste0("spot", 1:6)))
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = 2:5)
  expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  expect_equal(sum(colnames(s1ColData) == "sample_id"), 1L)
  expect_equal(s1ColData$sample_id, c("A", "A", "B", "B"))
  s1 <- getSubsetExperiment(es, "s1")
  expect_equal(sum(colnames(colData(s1)) == "sample_id"), 1L)

  # the main object's current sample_id, not the copy createSubset() made
  colData(es)$sample_id <- rep(c("C", "D"), each = 3)
  expect_message(
    s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE),
    NA
  )
  expect_equal(s1ColData$sample_id, c("C", "C", "D", "D"))
  s1 <- getSubsetExperiment(es, "s1")
  expect_equal(s1$sample_id, c("C", "C", "D", "D"))

  # writing back still works
  s1ColData$group <- letters[1:4]
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  s1$group <- letters[5:8]
  es <- suppressMessages(setSubsetExperiment(es, "s1", s1))
  expect_equal(colData(es, subsetName = "s1")$group, letters[5:8])
})

test_that("Subset data setters take NULL or a data.frame", {
  sce <- sce_chcl
  rowData(sce)$symbol <- rownames(sce)
  colData(sce)$batch <- rep(c("A", "B"), length.out = ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)

  rowData(es, subsetName = "s1") <- data.frame(rank = 5:1)
  expect_equal(rowData(es, subsetName = "s1")$rank, 5:1)
  rowData(es, subsetName = "s1") <- NULL
  expect_equal(ncol(rowData(es, subsetName = "s1")), 0L)
  expect_equal(rowData(es)$symbol, rownames(sce))

  colData(es, subsetName = "s1") <- data.frame(group = letters[1:4])
  expect_equal(colData(es, subsetName = "s1")$group, letters[1:4])
  colData(es, subsetName = "s1") <- NULL
  expect_equal(ncol(colData(es, subsetName = "s1")), 0L)
  expect_equal(subsetColnames(es, "s1"), colnames(sce)[1:4])
  expect_equal(colData(es)$batch, colData(sce)$batch)
  # the main object accepts only a DataFrame or NULL, as SummarizedExperiment
  expect_error(colData(es) <- data.frame(a = seq_len(ncol(es))),
               regexp = "DataFrame or NULL")

  # subsetRowData<- and subsetColData<- take the same values (#97)
  subsetRowData(es, "s1") <- DataFrame(rank = 5:1)
  fromDataFrame <- subsetRowData(es, "s1")
  subsetRowData(es, "s1") <- data.frame(rank = 5:1)
  expect_equal(subsetRowData(es, "s1"), fromDataFrame)
  expect_equal(subsetRowData(es, "s1")$rank, 5:1)
  expect_equal(subsetRownames(es, "s1"), rownames(sce)[1:5])
  subsetRowData(es, "s1") <- NULL
  expect_equal(ncol(subsetRowData(es, "s1")), 0L)
  subsetColData(es, "s1") <- DataFrame(group = letters[1:4])
  fromDataFrame <- subsetColData(es, "s1")
  subsetColData(es, "s1") <- data.frame(group = letters[1:4])
  expect_equal(subsetColData(es, "s1"), fromDataFrame)
  expect_equal(subsetColData(es, "s1")$group, letters[1:4])
  expect_equal(subsetColnames(es, "s1"), colnames(sce)[1:4])
  subsetColData(es, "s1") <- NULL
  expect_equal(ncol(subsetColData(es, "s1")), 0L)
  valueMsg <- "'value' must be a DataFrame, a data.frame or NULL."
  expect_error(subsetRowData(es, "s1") <- 1:5, valueMsg,
               fixed = TRUE)
  expect_error(subsetColData(es, "s1") <- 1:4, valueMsg,
               fixed = TRUE)
  expect_error(rowData(es, subsetName = "s1") <- 1:5, valueMsg,
               fixed = TRUE)

  # the parent-column mark is never stored
  s1RowData <- rowData(es, subsetName = "s1", parentRowData = TRUE)
  rowData(es) <- rowData(es)
  expect_null(S4Vectors::mcols(rowData(es))$fromParentData)
  s1ColData <- colData(es, subsetName = "s1", parentColData = TRUE)
  s1ColData$score <- 1:4
  suppressMessages(colData(es, subsetName = "s1") <- s1ColData)
  storedInfo <- S4Vectors::mcols(colData(es, subsetName = "s1"))
  expect_null(storedInfo$fromParentData)
})

test_that("colData<- rejects other values in every Subset class", {
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
    expect_error(colData(es, subsetName = "s1") <- 1:4,
                 subsetMsg, fixed = TRUE, info = name)
    expect_error(colData(es, subsetName = "s1") <- "a",
                 subsetMsg, fixed = TRUE, info = name)
    expect_error(colData(es) <- seq_len(ncol(es)),
                 mainMsg, fixed = TRUE, info = name)
    # a DataFrame or NULL with subsetName still sets the subset's colData,
    # not the main object's
    colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4])
    expect_equal(colData(es, subsetName = "s1")$group,
                 letters[1:4], info = name)
    expect_equal(colData(es), mainColData, info = name)
    colData(es, subsetName = "s1") <- NULL
    # a SpatialExperiment subset keeps its sample_id column
    keptColumns <- if (name == "spe") "sample_id" else character(0)
    expect_equal(colnames(colData(es, subsetName = "s1")),
                 keptColumns, info = name)
    expect_equal(colData(es), mainColData, info = name)
    # a data.frame with subsetName sets the subset's colData too
    colData(es, subsetName = "s1") <- data.frame(group = letters[5:8])
    expect_equal(colData(es, subsetName = "s1")$group,
                 letters[5:8], info = name)
    expect_equal(colData(es), mainColData, info = name)
  }
})

test_that("altExps and reducedDims of a subset pass arguments on", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  cells <- subsetColnames(es, "s1")
  colData(es, subsetName = "s1") <- DataFrame(group = letters[1:4],
                                              row.names = cells)
  altCounts <- matrix(1, nrow = 2, ncol = 4, dimnames = list(NULL, cells))
  altExp(es, "alt", subsetName = "s1") <-
    SingleCellExperiment(list(counts = altCounts))
  withGroup <- altExps(es, subsetName = "s1", withColData = TRUE)
  expect_equal(colData(withGroup$alt)$group, letters[1:4])
  withoutGroup <- altExps(es, subsetName = "s1")
  expect_null(colData(withoutGroup$alt)$group)

  pca <- matrix(1, nrow = 4, ncol = 2,
                dimnames = list(paste0("other", 1:4), NULL))
  reducedDim(es, "PCA", subsetName = "s1", withDimnames = FALSE) <- pca
  storedPcas <- reducedDims(es, subsetName = "s1", withDimnames = FALSE)
  expect_equal(rownames(storedPcas$PCA), paste0("other", 1:4))
})

test_that("Subset getters use renamed and repeated names", {
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
  expect_equal(rownames(links), paste0("g", 1:4))

  # spot barcodes repeat across the two samples
  barcodes <- rep(paste0("b", 1:3), 2)
  spe <- SpatialExperiment(assays = list(counts = counts[, 1:6]),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = cbind(x = 1:6, y = 1:6))
  colnames(spe) <- barcodes
  es <- ExperimentSubset(spe)
  es <- createSubset(es, "s1", cols = c(1, 2, 4, 5))
  # subsetSpatialData() is deprecated (#34)
  expect_warning(s1Data <- subsetSpatialData(es, "s1", FALSE),
                 "deprecated")
  expect_equal(rownames(s1Data), c("b1", "b2", "b1", "b2"))
  subsetColnames(es, "s1") <- paste0("spot", 1:4)
  expect_warning(s1Data <- subsetSpatialData(es, "s1", TRUE),
                 "deprecated")
  expect_equal(rownames(s1Data), paste0("spot", 1:4))
})

test_that("subsetParent lists each parent up to the main assay", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:10, cols = 1:10)
  es <- createSubset(es, "s2", rows = 1:5, cols = 1:5, parentAssay = "s1")
  expect_equal(subsetParent(es, "s2"), list("s1", "counts"))
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "stored")
  expect_equal(subsetParent(es, "s1/stored"), list("s1", "counts"))
  es <- createSubset(es, "s3", rows = 1:2, cols = 1:2,
                     parentAssay = "s1/stored")
  expect_equal(subsetParent(es, "s3"),
               list("s1/stored", "s1", "counts"))
})

test_that("Subsets see their parent's assays as linked views", {
  es <- ExperimentSubset(sce_chcl)
  assay(es, "logcounts") <- log1p(assay(es, "counts"))
  es <- createSubset(es, "filtered", rows = 1:20, cols = 1:10)
  es <- createSubset(es, "tcells", rows = 1:5, cols = 2:4,
                     parentAssay = "filtered")
  mainLog <- as.matrix(assay(es, "logcounts"))

  # a main assay is linked into every subset under the same name
  filteredLog <- as.matrix(assay(es, "filtered/logcounts"))
  expect_equal(filteredLog, mainLog[1:20, 1:10])
  tcellLog <- as.matrix(assay(es, "tcells", subsetAssayName = "logcounts"))
  expect_equal(tcellLog, mainLog[1:5, 2:4])
  expect_equal(subsetAssayNames(es, "tcells"),
               c("counts", "logcounts"))

  # an assay stored in filtered is linked into tcells
  filteredScaled <- filteredLog * 2
  es <- setSubsetAssay(es, "filtered", filteredScaled, "scaled")
  tcellView <- as.matrix(assay(es, "tcells/scaled"))
  expect_equal(tcellView, filteredScaled[1:5, 2:4])

  # storing over a linked name needs unlink = TRUE
  tcellScaled <- tcellLog + 1
  expect_error(setSubsetAssay(es, "tcells", tcellScaled, "scaled"),
               regexp = "linked to 'filtered/scaled'")
  expect_message(
    es <- setSubsetAssay(es, "tcells", tcellScaled, "scaled", unlink = TRUE),
    regexp = paste0("Subset 'tcells' now stores its own 'scaled', ",
                    "unlinked from 'filtered/scaled'."),
    fixed = TRUE
  )
  expect_equal(as.matrix(assay(es, "tcells/scaled")), tcellScaled)
  # the parent's assay is unchanged
  expect_equal(as.matrix(assay(es, "filtered/scaled")),
               filteredScaled)
  # a subset of tcells sees tcells' own scaled
  es <- createSubset(es, "cd4", rows = 1:2, cols = 1:2,
                     parentAssay = "tcells")
  expect_equal(as.matrix(assay(es, "cd4/scaled")),
               tcellScaled[1:2, 1:2])
  # main assays are linked too
  filteredCounts <- as.matrix(assay(es, "filtered"))
  counts <- filteredCounts
  expect_error(setSubsetAssay(es, "filtered", counts, "logcounts"),
               regexp = "linked to the main assay 'logcounts'")
})

test_that("Sibling subsets can store the same assay name", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  es <- createSubset(es, "s2", rows = 6:10, cols = 5:8)
  s1Norm <- as.matrix(assay(es, "s1")) + 1
  s2Norm <- as.matrix(assay(es, "s2")) + 2
  es <- setSubsetAssay(es, "s1", s1Norm, "norm")
  es <- setSubsetAssay(es, "s2", s2Norm, "norm")

  # a qualified name or subsetAssayName picks one
  expect_equal(as.matrix(assay(es, "s1/norm")), s1Norm)
  s2Stored <- assay(es, "s2", subsetAssayName = "norm")
  expect_equal(as.matrix(s2Stored), s2Norm)
  s1Stored <- getSubsetAssay(es, "s1", subsetAssayName = "norm")
  expect_equal(as.matrix(s1Stored), s1Norm)
  expect_equal(as.matrix(getSubsetAssay(es, "s2/norm")), s2Norm)
  # a bare name means a main assay or a subset
  expect_error(assay(es, "norm"),
               regexp = "stored in subset 's1', 's2'")
  es <- setSubsetAssay(es, "s1", s1Norm, "onlyInS1")
  expect_error(assay(es, "onlyInS1"), regexp = "s1/onlyInS1")
  expect_error(createSubset(es, "s3", 1:2, parentAssay = "norm"),
               regexp = "stored in subset")

  # a qualified parentAssay makes a subset of that stored assay
  es <- createSubset(es, "s3", rows = 1:2, cols = 1:2,
                     parentAssay = "s2/norm")
  expect_equal(subsetRownames(es, "s3"), rownames(sce_chcl)[6:7])
  expect_equal(as.matrix(assay(es, "s3")), s2Norm[1:2, 1:2])
  expect_equal(subsetParent(es, "s3"),
               list("s2/norm", "s2", "counts"))
  expect_true("s2/norm" %in% subsetAssayNames(es))

  # a later sibling storing "norm" doesn't change what s3 points to
  es <- createSubset(es, "s4", rows = 11:15, cols = 9:12)
  es <- setSubsetAssay(es, "s4", as.matrix(assay(es, "s4")), "norm")
  expect_equal(as.matrix(assay(es, "s3")), s2Norm[1:2, 1:2])

  # assay<- with a qualified name stores into that subset
  assay(es, "s1/extra") <- s1Norm
  expect_equal(as.matrix(assay(es, "s1", subsetAssayName = "extra")),
               s1Norm)
})

test_that("A subset's own assay is used by it and its descendants", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "P", rows = 1:5, cols = 1:4)
  es <- createSubset(es, "Q", rows = 1:5, cols = 1:4)
  mainCounts <- as.matrix(assay(es, "counts"))
  pCounts <- as.matrix(assay(es, "P")) + 100
  expect_message(
    es <- setSubsetAssay(es, "P", pCounts, "counts", unlink = TRUE),
    regexp = paste0("Subset 'P' now stores its own 'counts', unlinked ",
                    "from the main assay 'counts'."),
    fixed = TRUE
  )
  # P's primary assay is now its own counts
  expect_equal(as.matrix(assay(es, "P")), pCounts)
  es <- createSubset(es, "C", rows = 1:2, cols = 1:2, parentAssay = "P")
  expect_equal(as.matrix(assay(es, "C")), pCounts[1:2, 1:2])
  expect_equal(as.matrix(assay(es, "C/counts")), pCounts[1:2, 1:2])
  # each name is listed once, own or linked
  expect_equal(subsetAssayNames(es, "C"), assayNames(es))
  # the sibling Q still sees the main counts
  expect_equal(as.matrix(assay(es, "Q")), mainCounts[1:5, 1:4])
})

test_that("assay<- stores linked names only with unlink = TRUE", {
  es <- ExperimentSubset(sce_chcl)
  assay(es, "logcounts") <- log1p(assay(es, "counts"))
  es <- createSubset(es, "filtered", rows = 1:20, cols = 1:10)
  filteredLog <- as.matrix(assay(es, "filtered/logcounts"))
  expect_error(assay(es, "filtered/logcounts") <- filteredLog,
               regexp = "unlink = TRUE")
  expect_no_message(
    expect_error(
      assay(es, "filtered", subsetAssayName = "logcounts") <- filteredLog,
      regexp = "unlink = TRUE"
    )
  )
  expect_message(
    assay(es, "filtered/logcounts", unlink = TRUE) <- filteredLog + 1,
    regexp = paste0("Subset 'filtered' now stores its own 'logcounts', ",
                    "unlinked from the main assay 'logcounts'."),
    fixed = TRUE
  )
  expect_equal(as.matrix(assay(es, "filtered/logcounts")),
               filteredLog + 1)
  plusTwo <- filteredLog + 2
  suppressMessages(assay(es, "filtered", subsetAssayName = "logcounts",
                         unlink = TRUE) <- plusTwo)
  expect_equal(as.matrix(assay(es, "filtered/logcounts")),
               filteredLog + 2)
  expect_output(show(es), "filtered/logcounts")

  # a subset can be made from an assay its parent only links to
  es <- createSubset(es, "fromLinked", rows = 1:3, cols = 1:3,
                     parentAssay = "filtered/counts")
  mainCounts <- as.matrix(assay(es, "counts"))
  expect_equal(as.matrix(assay(es, "fromLinked")),
               mainCounts[1:3, 1:3])
  linkedView <- assay(es, "filtered/counts", withDimnames = FALSE)
  expect_equal(dim(linkedView), c(20L, 10L))

  # arguments that don't fit the call are errors
  expect_error(assay(es, "filtered/x", subsetAssayName = "y") <-
                 filteredLog, regexp = "already names")
  expect_error(assay(es, "counts", unlink = TRUE) <- mainCounts,
               regexp = "only applies")
  expect_error(setSubsetAssay(es, "filtered", filteredLog, ""),
               regexp = "non-empty")
})

test_that("Saved objects whose subsets record bare parents work", {
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

  expect_equal(as.matrix(assay(es, "child")), stored[1:2, 1:2])
  expect_equal(subsetParent(es, "child"),
               list("stored", "s1", "counts"))
  # a new subset called "stored" doesn't change what the child points to
  es <- createSubset(es, "stored", rows = 10:20, cols = 10:20)
  expect_equal(as.matrix(assay(es, "child")), stored[1:2, 1:2])
  expect_equal(subsetParent(es, "child"),
               list("s1/stored", "s1", "counts"))
})

test_that("getSubsetAssay reads stored assays in SE and RSE", {
  counts <- matrix(as.numeric(1:20), nrow = 4,
                   dimnames = list(paste0("g", 1:4), paste0("c", 1:5)))
  se <- SummarizedExperiment(list(counts = counts))
  rse <- as(se, "RangedSummarizedExperiment")
  for (object in list(se, rse)) {
    es <- ExperimentSubset(object)
    expect_true(is(es, paste0("Subset", class(object))))
    es <- createSubset(es, "s1", rows = 1:2, cols = 1:3)
    doubled <- counts[1:2, 1:3] * 2
    es <- setSubsetAssay(es, "s1", doubled, "doubled")
    s1Doubled <- getSubsetAssay(es, "s1", subsetAssayName = "doubled")
    expect_equal(as.matrix(s1Doubled), doubled)
  }
})

test_that("Qualified names that don't exist give clear errors", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  s1Counts <- as.matrix(assay(es, "s1"))
  expect_error(assay(es, "s1/none"), regexp = "s1/none")
  expect_error(assay(es, "s1/a/b"), regexp = "s1/a/b")
  expect_error(assay(es, "nope/counts"), regexp = "nope/counts")
  expect_error(assay(es, "s1", subsetAssayName = c("a", "b")),
               regexp = "subsetAssayName")
  # subsetAssayName with an i that isn't a subset name gives the same
  # reason as the setter, not that the subset doesn't exist (#76)
  expect_error(assay(es, "s1/counts", subsetAssayName = "x"),
               regexp = "already names")
  expect_error(assay(es, "counts", subsetAssayName = "x"),
               regexp = "not a subset")
  # a "/" name that isn't "subset/assay" for an existing subset, as for
  # the setter
  expect_error(assay(es, "nope/counts", subsetAssayName = "x"),
               regexp = "existing subset")
  expect_error(assay(es, "s1/a/b", subsetAssayName = "x"),
               regexp = "existing subset")
  # getSubsetAssay() names its own argument, subsetName, not assay()'s i
  expect_error(getSubsetAssay(es, "counts", subsetAssayName = "x"),
               regexp = "'subsetName' must be")
  expect_error(assay(es, "nope/extra") <- s1Counts,
               regexp = "existing subset")
  expect_error(subsetParent(es, "counts"),
               regexp = "does not exist")
})

test_that("Subset and stored assay names can't contain a slash", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  s1Counts <- as.matrix(assay(es, "s1"))
  expect_error(createSubset(es, "a/b", rows = 1:2), regexp = "/")
  expect_error(setSubsetAssay(es, "s1", s1Counts, "a/b"),
               regexp = "/")
})

test_that("Parent lookups follow .assay() for ambiguous names", {
  sce <- sce_chcl
  colData(sce)$cellIndex <- seq_len(ncol(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 101:200, cols = 501:600)
  es <- createSubset(es, "t1", rows = 1:10, cols = 1:10)
  # two subsets store an assay named "dup", so a plain "dup" is ambiguous
  # and the parent is named with a qualified name (#2)
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "dup")
  es <- setSubsetAssay(es, "t1", assay(es, "t1"), "dup")
  expect_error(createSubset(es, "s5", 2:3, 2:3, parentAssay = "dup"),
               regexp = "stored in subset")
  es <- createSubset(es, "s5", rows = 2:3, cols = 2:3, parentAssay = "t1/dup")
  expect_equal(subsetRownames(es, "s5"), rownames(assay(es, "s5")))
  expect_equal(subsetRownames(es, "s5"), rownames(sce)[2:3])

  # renaming the main assay must not break subsets made from it
  assayNames(es)[1] <- "raw"
  t1ColData <- colData(es, subsetName = "t1", parentColData = TRUE)
  expect_equal(t1ColData$cellIndex, 1:10)
})

test_that("Subsets with one row or column keep their dimensions", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "oneRow", rows = 5, cols = 1:10)
  es <- createSubset(es, "oneCol", rows = 1:10, cols = 5)
  expect_equal(dim(assay(es, "oneRow")), c(1L, 10L))
  expect_equal(dim(assay(es, "oneCol")), c(10L, 1L))
  # a subset made from a one-row subset
  es <- createSubset(es, "child", rows = 1, cols = 1:3, parentAssay = "oneRow")
  expect_equal(dim(assay(es, "child")), c(1L, 3L))
})

test_that("setSubsetAssay keeps the subset's row and column names", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:4)
  genes <- subsetRownames(es, "s1")
  cells <- subsetColnames(es, "s1")
  counts <- as.matrix(assay(es, "s1"))

  # input without names: the subset keeps its names, and so does the assay
  es <- setSubsetAssay(es, "s1", unname(counts), "unnamed")
  expect_equal(subsetRownames(es, "s1"), genes)
  expect_equal(subsetColnames(es, "s1"), cells)
  expect_equal(dimnames(assay(es, "s1/unnamed")), list(genes, cells))

  # input with the same names is accepted
  es <- setSubsetAssay(es, "s1", counts, "named")
  expect_equal(subsetColnames(es, "s1"), cells)

  # input with reordered or different names is rejected
  reordered <- counts[, c(2, 1, 3, 4)]
  expect_error(setSubsetAssay(es, "s1", reordered, "reordered"),
               regexp = "inputMatrix")
  renamed <- counts
  rownames(renamed) <- paste0("other", 1:5)
  expect_error(setSubsetAssay(es, "s1", renamed, "renamed"),
               regexp = "inputMatrix")
  # assay<- says it stored the assay only after setSubsetAssay() succeeds
  expect_no_message(
    expect_error(assay(es, "s1", subsetAssayName = "x") <- reordered,
                 regexp = "inputMatrix")
  )
  expect_equal(subsetColnames(es, "s1"), cells)
  expect_message(assay(es, "s1", subsetAssayName = "viaAssay") <- counts,
                 regexp = "Storing assay inside subset 's1'.", fixed = TRUE)
  expect_equal(assay(es, "s1/viaAssay"), counts)

  # a matrix of the wrong size gets a dimension error, not a names error
  expect_error(setSubsetAssay(es, "s1", counts[1:3, ], "small"),
               regexp = "same dimensions")

  # after renaming the subset's cells, its assay round-trips
  subsetColnames(es, "s1") <- paste0("cell", 1:4)
  s1Counts <- assay(es, "s1")
  expect_equal(colnames(s1Counts), paste0("cell", 1:4))
  es <- setSubsetAssay(es, "s1", s1Counts, "roundTrip")
  roundTrip <- assay(es, "s1/roundTrip")
  expect_equal(colnames(roundTrip), paste0("cell", 1:4))

  # assay<- says it stored the assay after saying it unlinked it
  storeMessages <- capture_messages(
    assay(es, "s1", subsetAssayName = "counts", unlink = TRUE) <-
      s1Counts + 1
  )
  expect_equal(storeMessages, c(
    paste0("Subset 's1' now stores its own 'counts', unlinked from the ",
           "main assay 'counts'.\n"),
    "Storing assay inside subset 's1'.\n"
  ))
})

test_that("A subset of a subset round-trips through setSubsetAssay", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "p", rows = 5:8, cols = 3:6)
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2, parentAssay = "p")
  childCounts <- assay(es, "child")
  expect_equal(rownames(childCounts), rownames(sce_chcl)[5:6])
  expect_equal(colnames(childCounts), colnames(sce_chcl)[3:4])
  es <- setSubsetAssay(es, "child", childCounts, "copy")
  childCopy <- assay(es, "child/copy")
  expect_equal(dimnames(childCopy), dimnames(childCounts))
})

test_that("A subset of an object without names takes the input's", {
  sce <- SingleCellExperiment(list(counts = matrix(1, 6, 5)))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:3, cols = 1:2)
  counts <- matrix(2, 3, 2, dimnames = list(paste0("g", 1:3), c("a", "b")))
  es <- setSubsetAssay(es, "s1", counts, "named")
  expect_equal(subsetRownames(es, "s1"), paste0("g", 1:3))
  expect_equal(subsetColnames(es, "s1"), c("a", "b"))
})

test_that("Spaces are removed from a new subset's name", {
  es <- ExperimentSubset(sce_chcl)
  spaced <- "my subset"
  expect_warning(es <- createSubset(es, spaced, 1:5, 1:5),
                 regexp = "Removing spaces")
  expect_equal(subsetNames(es), "mysubset")
  expect_equal(subsetDim(es, "mysubset"), c(5L, 5L))
  es <- setSubsetAssay(es, "mysubset", assay(es, "mysubset"), "stored")
  expect_equal(subsetParent(es, "mysubset/stored"),
               list("mysubset", "counts"))
  # the cleaned name can't overwrite an existing subset
  expect_error(suppressWarnings(createSubset(es, spaced, 1:2, 1:2)),
               regexp = "'mysubset' is already a subset name.",
               fixed = TRUE)
})

test_that("createSubset rejects invalid rows and cols", {
  es <- ExperimentSubset(sce_chcl)
  expect_error(createSubset(es, "s1", rows = -1),
               regexp = "positive whole numbers")
  expect_error(createSubset(es, "s1", rows = c(0, 1, 2)),
               regexp = "positive whole numbers")
  expect_error(createSubset(es, "s1", cols = 1.5),
               regexp = "positive whole numbers")
  expect_error(createSubset(es, "s1", rows = c(1, 1)),
               regexp = "duplicates")
  expect_error(createSubset(es, "s1", cols = c(2, 3, 2)),
               regexp = "duplicates")
  sameGene <- rownames(sce_chcl)[c(1, 1)]
  expect_error(createSubset(es, "s1", rows = sameGene),
               regexp = "duplicates")
  expect_error(createSubset(es, "s1", rows = integer(0)),
               regexp = "at least one")

  # names given for an assay without row or column names
  unnamedSce <- SingleCellExperiment(list(counts = matrix(1, 6, 5)))
  unnamed <- ExperimentSubset(unnamedSce)
  expect_error(createSubset(unnamed, "s1", rows = c("a", "b")),
               regexp = "no row names")
  expect_error(createSubset(unnamed, "s1", cols = "a"),
               regexp = "no column names")
})

test_that("createSubset rejects invalid names", {
  es <- ExperimentSubset(sce_chcl)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  es <- setSubsetAssay(es, "s1", assay(es, "s1"), "stored")

  expect_error(createSubset(es, c("a", "b"), rows = 1:2),
               regexp = "subsetName")
  expect_error(createSubset(es, NA_character_, rows = 1:2),
               regexp = "subsetName")
  expect_error(createSubset(es, "", rows = 1:2),
               regexp = "subsetName")
  # spaces are removed first, so a name of only spaces is empty
  expect_error(suppressWarnings(createSubset(es, "  ", rows = 1:2)),
               regexp = "subsetName")
  # a subset can't share its name with an assay
  expect_error(createSubset(es, "counts", rows = 1:2),
               regexp = "already an assay name")
  # a stored assay is named "s1/stored", so a subset may be called "stored"
  es <- createSubset(es, "stored", rows = 1:2)
  expect_true("stored" %in% subsetNames(es))
  twoParents <- c("counts", "s1")
  expect_error(createSubset(es, "s2", 1:2, parentAssay = twoParents),
               regexp = "parentAssay")
})

test_that("ExperimentSubset checks the subset argument and lists", {
  withNA <- list(subsetName = "s1", rows = 1:10, cols = NA, parentAssay = NA)
  expect_error(ExperimentSubset(sce_chcl, subset = withNA),
               regexp = "'subset'")
  expect_error(ExperimentSubset(list(matrix(1, 3, 3))),
               regexp = "unnamed assays")

  # leaving out cols and parentAssay keeps all columns and the first assay
  es <- ExperimentSubset(sce_chcl, subset = list(subsetName = "s1",
                                                 rows = 1:10))
  expect_equal(subsetDim(es, "s1"), c(10L, ncol(sce_chcl)))
})

test_that("Unknown subset names give a clear error", {
  es <- ExperimentSubset(sce_chcl)
  # no subsets yet
  expect_error(assay(es, "nope"), regexp = "nope")
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  expect_error(assay(es, "nope"), regexp = "nope")
  expect_error(getSubsetAssay(es, "nope"), regexp = "nope")
  expect_error(subsetDim(es, "nope"),
               regexp = "nope subset does not exist")
  expect_equal(
    tryCatch(subsetParent(es, "nope"), error = conditionMessage),
    paste0("nope subset does not exist. 'subsetName' must be the name of ",
           "a subset or a \"subset/assay\".")
  )
  expect_error(subsetColnames(es, "nope"),
               regexp = "nope subset does not exist")
  expect_error(subsetRownames(es, "nope"),
               regexp = "nope subset does not exist")
  expect_error(subsetColnames(es, "nope") <- letters[1:5],
               regexp = "nope subset does not exist")
  expect_error(subsetRownames(es, "nope") <- letters[1:5],
               regexp = "nope subset does not exist")
  expect_error(reducedDims(es, subsetName = "nope") <- list(),
               regexp = "nope subset does not exist")
})

test_that("assay<- stores into a subset chosen by its name", {
  es <- ExperimentSubset(sce_chcl)
  # a subset with the same dimensions as the main assay
  es <- createSubset(es, "all")
  counts <- assay(es, "counts")
  suppressMessages(assay(es, "all", subsetAssayName = "copy") <- counts)
  expect_true("all/copy" %in% subsetAssayNames(es))
  expect_false("all" %in% assayNames(es))

  es <- createSubset(es, "s1", rows = 1:5, cols = 1:5)
  expect_error(assay(es, "s1") <- assay(es, "s1"),
               regexp = "is needed to store")
  # subsetAssayName with a name that is not a subset
  s1Counts <- assay(es, "s1")
  expect_error(assay(es, "s1typo", subsetAssayName = "n") <- s1Counts,
               regexp = "not a subset")
})

test_that("Testing createSubset with multiple paramter options",{
  es <- ExperimentSubset(sce_chcl)
  expect_equal(
    tryCatch(createSubset(es, "subset1", rows = 1:5, cols = 1:5,
                          parentAssay = "none"),
             error = conditionMessage),
    paste0("'parentAssay' must be NULL or name an assay, a subset or a ",
           "\"subset/assay\"; 'none' is none of these.")
  )
  
  #colnames
  es <- createSubset(es, "subset1", rows = c(10,11,50,56,98,99,102,105,109, 200), cols = c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC"), parentAssay = "counts")
  namedCols <- c("CTGCTGTCAGGGTATG", "CAGTCCTTCGGTTAAC")
  expect_equal(subsetNames(es), "subset1")
  expectSubsetOf(es, "subset1",
                 rows = c(10, 11, 50, 56, 98, 99, 102, 105, 109, 200),
                 cols = match(namedCols, colnames(es)))
  
  #null rows and null cols
  es <- createSubset(es, "subset2", parentAssay = "counts")
  expect_equal(subsetNames(es), c("subset1", "subset2"))
  expectSubsetOf(es, "subset2", rows = seq_len(nrow(es)),
                 cols = seq_len(ncol(es)))
  
  # rows or columns given by names that the parent doesn't have
  expect_error(
    createSubset(es, "subset22", rows = c("as", "asd", "asd"),
                 cols = namedCols, parentAssay = "counts"),
    regexp = "Selected rows not available in the specified assay."
  )
  expect_error(
    createSubset(es, "subset22", cols = c("CTGTG", "CAGTCC"),
                 parentAssay = "counts"),
    regexp = "Selected columns not available in the specified assay."
  )
})

test_that("Completing test coverage, move this into relevant places later",{
  es <- ExperimentSubset(sce_chcl)
  expect_warning(es <- createSubset(es, "subset 1", rows = c(1:5), cols = (1:5)),
                 "Removing spaces from the specified subsetName.")
  expect_equal(subsetNames(es), "subset1")
  
  # a named list of assays gives a SubsetSingleCellExperiment of them
  es <- ExperimentSubset(list(counts = assay(sce_chcl, "counts")))
  expect_s4_class(es, "SubsetSingleCellExperiment")
  expect_equal(assayNames(es), "counts")
  expect_equal(assay(es, "counts"), assay(sce_chcl, "counts"))
  expect_equal(subsetCount(es), 0)
})
# Takes a subset out as a standalone experiment and puts it back unchanged
checkExperimentRoundTrip <- function(object) {
  object$batch <- rep(c("a", "b"), 3)
  es <- ExperimentSubset(object)
  es <- createSubset(es, "s1", rows = 2:5, cols = 2:5)
  sub <- getSubsetExperiment(es, "s1")
  expect_equal(class(sub), class(object))
  expect_equal(dim(sub), c(4L, 4L))
  expect_equal(rownames(sub), subsetRownames(es, "s1"))
  expect_equal(colnames(sub), subsetColnames(es, "s1"))
  expect_equal(as.matrix(assay(sub, "counts")),
               as.matrix(assay(es, "s1")))
  expect_equal(sub$batch, object$batch[2:5])
  # unchanged linked assays and parent columns are not stored
  ownColumns <- colnames(subsetColData(es, "s1", FALSE))
  es <- setSubsetExperiment(es, "s1", sub)
  expect_equal(subsetAssayNames(es), "s1")
  expect_equal(colnames(subsetColData(es, "s1", FALSE)),
               ownColumns)
  invisible(sub)
}

test_that("Subset experiments round-trip for every class", {
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
  expect_false(is.null(rowTree(treeSub)))
  expect_equal(rowLinks(treeSub)$nodeLab, rowLinks(tse)$nodeLab[2:5])

  coords <- cbind(x = as.numeric(1:6), y = as.numeric(1:6) * 10)
  spe <- SpatialExperiment(assays = list(counts = counts),
                           sample_id = rep(c("A", "B"), each = 3),
                           spatialCoords = coords)
  spatialSub <- checkExperimentRoundTrip(spe)
  expect_equal(unname(spatialCoords(spatialSub)),
               unname(coords[2:5, ]))
  expect_equal(spatialSub$sample_id, spe$sample_id[2:5])
})

test_that("setSubsetExperiment stores what a tool added", {
  sce <- simSce
  sce$batch <- rep(c("a", "b"), 60)
  rowData(sce)$symbol <- toupper(rownames(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  # parent columns are included and tagged
  expect_equal(sub$batch, sce$batch[1:10])
  expect_equal(rowData(sub)$symbol, rowData(sce)$symbol[1:20])

  scaled <- as.matrix(assay(sub, "counts")) / 2
  assay(sub, "scaled") <- scaled
  pca <- matrix(as.numeric(1:20), nrow = 10,
                dimnames = list(colnames(sub), c("PC1", "PC2")))
  reducedDim(sub, "PCA") <- pca
  sub$cluster <- rep(1:2, 5)
  rowData(sub)$variable <- rep(c(TRUE, FALSE), 10)
  metadata(sub)$note <- "clustered"
  es <- setSubsetExperiment(es, "s1", sub)

  expect_equal(subsetAssayNames(es), c("s1", "s1/scaled"))
  expect_equal(as.matrix(assay(es, "s1/scaled")), scaled)
  expect_equal(reducedDim(es, "PCA", subsetName = "s1"), pca)
  expect_equal(colnames(colData(es, subsetName = "s1")), "cluster")
  expect_equal(colData(es, subsetName = "s1")$cluster, rep(1:2, 5))
  expect_equal(colnames(rowData(es, subsetName = "s1")),
               "variable")
  expect_equal(metadata(es, subsetName = "s1")$note, "clustered")
  # the main object is unchanged
  expect_equal(colnames(colData(es)), "batch")
  expect_equal(assayNames(es), "counts")

  # without parent data only the subset's own columns come back
  own <- getSubsetExperiment(es, "s1", parentData = FALSE)
  expect_equal(colnames(colData(own)), "cluster")
  expect_equal(colnames(rowData(own)), "variable")
  expect_equal(assayNames(own), c("counts", "scaled"))
})

test_that("setSubsetExperiment replaces the subset's own data", {
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
  expect_equal(subsetAssayNames(es), "s1")
  expect_equal(reducedDimNames(es, subsetName = "s1"), character(0))
})

test_that("Changed linked assays are stored only with unlink", {
  sce <- simSce
  logcounts(sce) <- log1p(counts(sce))
  es <- ExperimentSubset(sce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  expect_equal(assayNames(sub), c("counts", "logcounts"))
  # unchanged, so not copied
  es <- setSubsetExperiment(es, "s1", sub)
  expect_equal(subsetAssayNames(es), "s1")

  changed <- as.matrix(logcounts(sub)) * 2
  logcounts(sub) <- changed
  expect_error(setSubsetExperiment(es, "s1", sub),
               regexp = "logcounts.*unlink = TRUE")
  expect_message(
    es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE),
    regexp = paste0("Subset 's1' now stores its own 'logcounts', ",
                    "unlinked from the main assay 'logcounts'."),
    fixed = TRUE
  )
  expect_equal(subsetAssayNames(es), c("s1", "s1/logcounts"))
  expect_equal(as.matrix(assay(es, "s1/logcounts")), changed)
  # the main object's logcounts is unchanged
  expect_equal(assay(es, "logcounts"), logcounts(sce))
})

test_that("setSubsetExperiment checks value against the subset", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  expect_error(setSubsetExperiment(es, "s1", sub[1:5, ]),
               regexp = "dimensions")
  expect_error(setSubsetExperiment(es, "s1", sub[, 10:1]),
               regexp = "colnames")
  expect_error(setSubsetExperiment(es, "s1", sub[20:1, ]),
               regexp = "rownames")
  expect_error(setSubsetExperiment(es, "none", sub),
               regexp = "none")
  expect_error(getSubsetExperiment(es, "none"), regexp = "none")
  se <- SummarizedExperiment(list(counts = counts(sub)))
  expect_error(setSubsetExperiment(es, "s1", se),
               regexp = "SingleCellExperiment")
  # an SCE for a SummarizedExperiment main object
  seEs <- ExperimentSubset(SummarizedExperiment(list(counts = simCounts)))
  seEs <- createSubset(seEs, "s1", rows = 1:20, cols = 1:10)
  expect_error(setSubsetExperiment(seEs, "s1", sub),
               regexp = "must be a SummarizedExperiment")
  bad <- sub
  assay(bad, "a/b") <- counts(sub)
  expect_error(setSubsetExperiment(es, "s1", bad), regexp = "/")
})

test_that("Subset experiments of subsets use the subset's assays", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  scaled <- as.matrix(assay(es, "s1")) / 2
  es <- setSubsetAssay(es, "s1", scaled, "scaled")
  es <- createSubset(es, "s2", rows = 1:5, cols = 1:4,
                     parentAssay = "s1/scaled")
  sub <- getSubsetExperiment(es, "s2")
  # the subset's primary assay comes first
  expect_equal(assayNames(sub), c("scaled", "counts"))
  expect_equal(as.matrix(assay(sub, "scaled")), scaled[1:5, 1:4])
  expect_equal(as.matrix(assay(sub)), as.matrix(assay(es, "s2")))
  es <- setSubsetExperiment(es, "s2", sub)
  expect_equal(subsetAssayNames(es), c("s1", "s2", "s1/scaled"))
  # s1's scaled can't be dropped while s2 is made from it
  s1 <- getSubsetExperiment(es, "s1")
  assay(s1, "scaled") <- NULL
  expect_error(setSubsetExperiment(es, "s1", s1),
               regexp = "'scaled', which subsets 's2'")
  expect_equal(as.matrix(assay(es, "s2")), scaled[1:5, 1:4])
})

test_that("setSubsetExperiment keeps assays that descendants use", {
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
  expect_error(setSubsetExperiment(es, "s1", s1),
               regexp = "'scaled', which subsets 'grandchild'")

  # a subset made from s1 uses s1's own counts, stored with unlink = TRUE
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  ownCounts <- as.matrix(assay(es, "s1")) * 100
  expect_message(
    es <- setSubsetAssay(es, "s1", ownCounts, "counts", unlink = TRUE),
    regexp = paste0("Subset 's1' now stores its own 'counts', unlinked ",
                    "from the main assay 'counts'."),
    fixed = TRUE
  )
  es <- createSubset(es, "child", rows = 1:5, cols = 1:4, parentAssay = "s1")
  s1 <- getSubsetExperiment(es, "s1")
  assay(s1, "counts") <- NULL
  expect_error(setSubsetExperiment(es, "s1", s1),
               regexp = "'counts', which subsets 'child'")
  expect_equal(as.matrix(assay(es, "child")), ownCounts[1:5, 1:4])
})

test_that("Subset experiments keep names a subset was given", {
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
  expect_equal(rownames(sub), paste0("row", 1:5))
  expect_equal(colnames(sub), paste0("col", 1:4))
  expect_equal(rowLinks(sub)$nodeLab, tinyTree$tip.label[1:5])
  es <- setSubsetExperiment(es, "s1", sub)
  expect_equal(subsetRownames(es, "s1"), paste0("row", 1:5))
})

test_that("setSubsetExperiment needs the subset's names on value", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  noColnames <- sub
  colnames(noColnames) <- NULL
  expect_error(setSubsetExperiment(es, "s1", noColnames),
               regexp = "colnames")
  noRownames <- sub
  rownames(noRownames) <- NULL
  expect_error(setSubsetExperiment(es, "s1", noRownames),
               regexp = "rownames")
})

test_that("Linked assays with the same values count as unchanged", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:20, cols = 1:10)
  sub <- getSubsetExperiment(es, "s1")
  # integer counts become double with the same values
  counts(sub) <- counts(sub) + 0
  expect_type(counts(sub), "double")
  es <- setSubsetExperiment(es, "s1", sub)
  expect_equal(subsetAssayNames(es), "s1")
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

test_that("[ keeps the rows of subsets that remain (#4)", {
  example <- bracketExample()
  es <- example$es
  counts <- example$counts

  kept <- es[3:10, ]
  expect_equal(as.matrix(assay(kept, "s1")), counts[3:5, 1:8])
  expect_equal(as.matrix(assay(kept, "s1/doubled")),
               counts[3:5, 1:8] * 2)
  expect_equal(as.matrix(assay(kept, "s1a")), counts[3:4, 3:6])
  expect_equal(as.matrix(assay(kept, "s2")), counts[6:10, 9:20])
  expect_equal(subsetRownames(kept, "s1"), paste0("g", 3:5))

  # logical, negative and name subscripts select the same rows
  rowsKept <- rep(c(FALSE, TRUE), c(2, 8))
  expect_equal(assay(es[rowsKept, ], "s1a"), assay(kept, "s1a"))
  expect_equal(assay(es[-(1:2), ], "s1a"), assay(kept, "s1a"))
  expect_equal(assay(es[paste0("g", 3:10), ], "s1a"),
               assay(kept, "s1a"))

  # with no subscripts nothing changes
  expect_equal(es[, ], es)

  # a subset keeps its own order when the main object is reordered
  reversed <- es[10:1, ]
  expect_equal(as.matrix(assay(reversed, "s1")), counts[1:5, 1:8])
  expect_equal(as.matrix(assay(reversed, "s1a")), counts[2:4, 3:6])
  expect_equal(as.matrix(assay(reversed[1:3, ], "s1")),
               counts[character(0), 1:8, drop = FALSE])
})

test_that("[ keeps the columns of subsets that remain (#4)", {
  example <- bracketExample()
  es <- example$es
  counts <- example$counts

  kept <- es[, c(2, 5, 9, 10)]
  expect_equal(as.matrix(assay(kept, "s1")), counts[1:5, c(2, 5)])
  expect_equal(as.matrix(assay(kept, "s1/doubled")),
               counts[1:5, c(2, 5)] * 2)
  expect_equal(as.matrix(assay(kept, "s1a")),
               counts[2:4, 5, drop = FALSE])
  expect_equal(as.matrix(assay(kept, "s2")), counts[6:10, 9:10])
  expect_equal(colData(kept, subsetName = "s1")$tag, c("t2", "t5"))
  expect_equal(reducedDim(kept, "PCA", subsetName = "s1"),
               example$pca[c(2, 5), ])

  both <- es[3:10, c(2, 5, 9, 10)]
  expect_equal(as.matrix(assay(both, "s1")), counts[3:5, c(2, 5)])
  expect_equal(subsetDim(both, "s1a"), c(2L, 1L))
})

test_that("[ keeps subsets left with no rows or columns (#4)", {
  example <- bracketExample()
  es <- example$es

  firstRows <- es[1:4, ]
  expect_equal(subsetNames(firstRows), subsetNames(es))
  expect_equal(dim(assay(firstRows, "s2")), c(0L, 12L))
  expect_equal(subsetDim(firstRows, "s2"), c(0L, 12L))
  expect_equal(dim(assay(firstRows, "s1a")), c(3L, 4L))

  lastCols <- es[, 9:20]
  expect_equal(dim(assay(lastCols, "s1")), c(5L, 0L))
  expect_equal(dim(assay(lastCols, "s1a")), c(3L, 0L))
  expect_equal(dim(assay(lastCols, "s2")), c(5L, 12L))
})

test_that("[ remaps subsets in every Subset class (#4)", {
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
    expect_s4_class(kept, class(es))
    expect_equal(as.matrix(assay(kept, "s1")), counts[3:5, 2:3])
    # spatial coordinates of the subset follow its columns
    if (is(experiment, "SpatialExperiment"))
      expect_equal(unname(subsetSpatialCoords(kept, "s1")[, "x"]),
                   c(2, 3))
  }
})

test_that("[ keeps one copy of repeated rows in subsets (#4)", {
  example <- bracketExample()
  repeated <- example$es[c(3, 3, 4, 1), ]
  expect_equal(subsetRownames(repeated, "s1"),
               paste0("g", c(1, 3, 4)))
  expect_equal(as.matrix(assay(repeated, "s1")),
               example$counts[c(1, 3, 4), 1:8])
  expect_equal(as.matrix(assay(repeated, "s1a")),
               example$counts[3:4, 3:6])
})

test_that("[ remaps subsets made from a stored assay (#4)", {
  example <- bracketExample()
  es <- createSubset(example$es, "fromDoubled", rows = 2:5, cols = 1:4,
                     parentAssay = "s1/doubled")
  kept <- es[3:10, 2:20]
  expect_equal(as.matrix(assay(kept, "fromDoubled")),
               example$counts[3:5, 2:4] * 2)
})

test_that("createSubset says when the parent is empty (#68)", {
  example <- bracketExample()
  es <- example$es

  noRows <- es[6:10, ]
  expect_error(createSubset(noRows, "child", parentAssay = "s1"),
               "Parent 's1' has no rows")
  expect_error(createSubset(noRows, "child", rows = 1,
                            parentAssay = "s1"),
               "Parent 's1' has no rows")
  expect_error(createSubset(noRows, "child", rows = "g1",
                            parentAssay = "s1"),
               "Parent 's1' has no rows")

  noCols <- es[, 9:20]
  expect_error(createSubset(noCols, "child", parentAssay = "s1"),
               "Parent 's1' has no columns")
  expect_error(createSubset(noCols, "child", cols = 1,
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

test_that("cbind keeps the subsets of every object (#67)", {
  example <- combineExample("cols")
  counts1 <- example$counts1
  counts2 <- example$counts2

  combined <- cbind(example$es1, example$es2)
  expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  expect_equal(as.matrix(assay(combined, "s2")),
               counts2[6:10, 9:20])
  expect_equal(subsetColnames(combined, "s2"), paste0("d", 9:20))
  # the subset's own assay and colData, and the subset of the subset
  expect_equal(as.matrix(assay(combined, "s2/doubled")),
               counts2[6:10, 9:20] * 2)
  expect_equal(colData(combined, subsetName = "s2")$tag,
               paste0("t", 1:12))
  expect_equal(as.matrix(assay(combined, "s2a")),
               counts2[7:9, 11:14])
  # a subset made from an assay stored in a subset
  expect_equal(as.matrix(assay(combined, "s2d")),
               counts2[6:7, 9:11] * 2)
})

test_that("rbind keeps the subsets of every object (#67)", {
  example <- combineExample("rows")
  counts1 <- example$counts1
  counts2 <- example$counts2

  combined <- rbind(example$es1, example$es2)
  expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  expect_equal(as.matrix(assay(combined, "s2")),
               counts2[6:10, 9:20])
  expect_equal(subsetRownames(combined, "s2"), paste0("h", 6:10))
  expect_equal(as.matrix(assay(combined, "s2/doubled")),
               counts2[6:10, 9:20] * 2)
  expect_equal(colData(combined, subsetName = "s2")$tag,
               paste0("t", 1:12))
  expect_equal(as.matrix(assay(combined, "s2a")),
               counts2[7:9, 11:14])
  expect_equal(as.matrix(assay(combined, "s2d")),
               counts2[6:7, 9:11] * 2)
})

test_that("cbind counts objects without subsets (#67)", {
  example <- combineExample("cols")
  counts3 <- example$counts1 * 10000
  colnames(counts3) <- paste0("e", 1:20)
  noSubsets <- ExperimentSubset(SingleCellExperiment(list(counts = counts3)))

  combined <- cbind(example$es1, noSubsets, example$es2)
  expect_equal(dim(combined), c(10L, 60L))
  expect_equal(as.matrix(assay(combined, "s2")),
               example$counts2[6:10, 9:20])
  expect_equal(as.matrix(assay(combined, "s2a")),
               example$counts2[7:9, 11:14])
})

test_that("cbind and rbind stop on repeated subset names (#67)", {
  example <- combineExample("cols")
  expect_error(cbind(example$es1, example$es1),
               "unique.*'s1' is in more than one.*different names")
  rowsExample <- combineExample("rows")
  expect_error(rbind(rowsExample$es2, rowsExample$es2),
               "'s2', 's2a', 's2d' are in more than one")
})

test_that("cbind and rbind stop when names don't line up (#67)", {
  example <- combineExample("cols")
  # the second object's genes in the opposite order
  reversed <- example$es2[10:1, ]
  expect_error(cbind(example$es1, reversed),
               "same rows.*object 2")
  rowsExample <- combineExample("rows")
  reversed <- rowsExample$es2[, 20:1]
  # rbind() of the class the objects extend already checks the columns
  expect_error(rbind(rowsExample$es1, reversed),
               "same colnames")
})

test_that("cbind updates the sample_id of spatial subsets (#67)", {
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
    expect_equal(subsetColData(combined, name)$sample_id,
                 combined$sample_id[rootCols])
  }
  expect_false(identical(combined$sample_id[9:10],
                         spe2$sample_id[3:4]))
})

test_that("All Subset classes keep subsets in cbind, rbind (#67)", {
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
    expect_s4_class(combined, class(es1))
    expect_equal(subsetNames(combined), c("s1", "s2"))
    expect_equal(as.matrix(assay(combined, "s2")),
                 counts2[6:10, 4:6])

    es3 <- createSubset(ExperimentSubset(third[[k]]), "s3",
                        rows = 6:10, cols = 4:6)
    stacked <- rbind(es1, es3)
    expect_s4_class(stacked, class(es1))
    expect_equal(subsetNames(stacked), c("s1", "s3"))
    expect_equal(as.matrix(assay(stacked, "s3")),
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
  expect_equal(subsetNames(combined), c("s1", "s2", "s2a", "s2d"))
  expect_equal(as.matrix(assay(combined, "s1")), counts1[1:5, 1:8])
  expect_equal(as.matrix(assay(combined, "s2")),
               counts2[1:6, 9:20])
  expect_equal(subsetRownames(combined, "s2"),
               rownames(counts2)[1:6])
  expect_equal(subsetColnames(combined, "s2"),
               colnames(counts2)[9:20])
  # the subset's own assay and colData, and the subset of the subset
  expect_equal(as.matrix(assay(combined, "s2/doubled")),
               counts2[1:6, 9:20] * 2)
  expect_equal(colData(combined, subsetName = "s2")$tag,
               paste0("t", 1:12))
  expect_equal(as.matrix(assay(combined, "s2a")),
               counts2[2:4, 11:14])
  # a subset made from an assay stored in a subset
  expect_equal(as.matrix(assay(combined, "s2d")),
               counts2[1:2, 9:11] * 2)
}

test_that("combineCols keeps the subsets of every object (#73)", {
  example <- unionExample("cols")
  combined <- combineCols(example$es1, example$es2)
  expect_s4_class(combined, "SubsetSingleCellExperiment")
  expect_equal(nrow(combined), 14L)
  expectUnionSubsets(combined, example)
})

test_that("combineRows keeps the subsets of every object (#73)", {
  example <- unionExample("rows")
  combined <- combineRows(example$es1, example$es2)
  # combineRows() of SingleCellExperiment returns the class it extends
  expect_s4_class(combined, "SubsetRangedSummarizedExperiment")
  expect_equal(ncol(combined), 28L)
  expectUnionSubsets(combined, example)
})

test_that("combineCols, combineRows match by position (#73)", {
  example <- unionExample("cols")
  # the rows have different names, so the parent method warns that it
  # keeps the first object's rowRanges
  combined <- suppressWarnings(combineCols(example$es1, example$es2,
                                           use.names = FALSE))
  expect_equal(nrow(combined), 10L)
  # the rows keep their positions and take the first object's names
  expect_equal(subsetRownames(combined, "s2"),
               rownames(example$counts1)[1:6])
  expect_equal(unname(as.matrix(assay(combined, "s2"))),
               unname(example$counts2[1:6, 9:20]))

  example <- unionExample("rows")
  combined <- combineRows(example$es1, example$es2, use.names = FALSE)
  expect_equal(ncol(combined), 20L)
  expect_equal(subsetColnames(combined, "s2"),
               colnames(example$counts1)[9:20])
  expect_equal(unname(as.matrix(assay(combined, "s2"))),
               unname(example$counts2[1:6, 9:20]))
})

test_that("combineCols counts objects without subsets (#73)", {
  example <- unionExample("cols")
  counts3 <- example$counts1 * 10000
  rownames(counts3) <- paste0("k", 1:10)
  colnames(counts3) <- paste0("e", 1:20)
  noSubsets <- ExperimentSubset(SingleCellExperiment(list(counts = counts3)))

  combined <- combineCols(example$es1, noSubsets, example$es2)
  expect_equal(dim(combined), c(24L, 60L))
  expectUnionSubsets(combined, example)
  # a plain object, not converted with ExperimentSubset()
  plain <- SingleCellExperiment(list(counts = counts3))
  combined <- combineCols(example$es1, plain, example$es2)
  expect_equal(dim(combined), c(24L, 60L))
  expectUnionSubsets(combined, example)
})

test_that("combineCols, combineRows stop on repeated names (#73)", {
  example <- unionExample("cols")
  expect_error(combineCols(example$es1, example$es1),
               "unique.*'s1' is in more than one.*different names")
  rowsExample <- unionExample("rows")
  expect_error(combineRows(rowsExample$es2, rowsExample$es2),
               "'s2', 's2a', 's2d' are in more than one")
})

test_that("All Subset classes keep subsets in combineCols (#73)", {
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
    expect_s4_class(combined, colsClass[k])
    expect_equal(subsetNames(combined), c("s1", "s2"))
    expect_equal(as.matrix(assay(combined, "s1")), counts[2:5, 1:3])
    expect_equal(as.matrix(assay(combined, "s2")),
                 counts2[1:6, 4:6])
    for (name in subsetNames(combined)) {
      subset <- ExperimentSubset:::.subsets(combined)[[name]]
      internal <- ExperimentSubset:::.internalAssay(subset)
      expect_equal(class(internal)[1],
                   sub("Subset", "", colsClass[k]))
    }

    es3 <- createSubset(ExperimentSubset(third[[k]]), "s3",
                        rows = 6:10, cols = 1:4)
    stacked <- combineRows(es1, es3)
    expect_s4_class(stacked, rowsClass[k])
    expect_equal(subsetNames(stacked), c("s1", "s3"))
    expect_equal(as.matrix(assay(stacked, "s1")), counts[2:5, 1:3])
    expect_equal(as.matrix(assay(stacked, "s3")),
                 counts3[6:10, 1:4])
    for (name in subsetNames(stacked)) {
      subset <- ExperimentSubset:::.subsets(stacked)[[name]]
      internal <- ExperimentSubset:::.internalAssay(subset)
      expect_equal(class(internal)[1],
                   sub("Subset", "", rowsClass[k]))
    }
  }
})

test_that("Combining keeps names set with subsetRownames<- (#73)", {
  example <- combineExample("cols")
  es2 <- example$es2
  subsetRownames(es2, "s2") <- paste0("renamed", 1:5)
  combined <- cbind(example$es1, es2)
  expect_equal(subsetRownames(combined, "s2"),
               paste0("renamed", 1:5))

  example <- unionExample("cols")
  es2 <- example$es2
  ownNames <- c("custom", rownames(example$counts2)[2:6])
  subsetRownames(es2, "s2") <- ownNames
  combined <- combineCols(example$es1, es2)
  expect_equal(subsetRownames(combined, "s2"), ownNames)
  # matched by position, the renamed row keeps its name and the others
  # take the first object's
  combined <- suppressWarnings(combineCols(example$es1, es2,
                                           use.names = FALSE))
  expect_equal(subsetRownames(combined, "s2"),
               c("custom", rownames(example$counts1)[2:6]))

  example <- unionExample("rows")
  es2 <- example$es2
  subsetColnames(es2, "s2") <- c("custom", colnames(example$counts2)[10:20])
  combined <- combineRows(example$es1, es2, use.names = FALSE)
  expect_equal(subsetColnames(combined, "s2"),
               c("custom", colnames(example$counts1)[10:20]))
})

test_that("Combining stops on a subset named as an assay (#73)", {
  example <- unionExample("cols")
  # the second object has a main assay named as the first object's subset
  es2 <- example$es2
  assay(es2, "s1") <- assay(es2, "counts") * 2
  expect_error(combineCols(example$es1, es2),
               "'s1' is a subset in one object and an assay")
  example <- unionExample("rows")
  es1 <- example$es1
  assay(es1, "s2") <- assay(es1, "counts") * 2
  expect_error(combineRows(es1, example$es2),
               "'s2' is a subset in one object and an assay")
})

test_that("use.names = FALSE names subsets of unnamed data (#73)", {
  example <- unionExample("cols")
  # the second object's genes have no names
  counts2 <- unname(example$counts2)
  colnames(counts2) <- paste0("d", 1:20)
  es2 <- ExperimentSubset(SingleCellExperiment(list(counts = counts2)))
  es2 <- createSubset(es2, "s2", rows = 1:6, cols = 9:20)
  combined <- suppressWarnings(combineCols(example$es1, es2,
                                           use.names = FALSE))
  expect_equal(subsetRownames(combined, "s2"),
               rownames(example$counts1)[1:6])

  # the first object's genes have no names, so the result takes the
  # second object's, and the subset keeps them
  counts1 <- unname(example$counts1)
  colnames(counts1) <- paste0("c", 1:20)
  es1 <- ExperimentSubset(SingleCellExperiment(list(counts = counts1)))
  combined <- suppressWarnings(combineCols(es1, example$es2,
                                           use.names = FALSE))
  expect_equal(subsetRownames(combined, "s2"),
               rownames(combined)[1:6])
})

test_that("Combining allows a subset named as its own assay (#73)", {
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
  expect_true("extra" %in% subsetNames(combined))
})

test_that("A subsetName that isn't one string gives a clear error", {
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
    expect_error(metadata(es, subsetName = bad), msg, info = info)
    expect_error(metadata(es, subsetName = bad) <- list(),
                 msg, info = info)
    expect_error(altExp(es, subsetName = bad), msg, info = info)
    expect_error(altExp(es, "a", subsetName = bad) <- es,
                 msg, info = info)
    expect_error(altExps(es, subsetName = bad), msg, info = info)
    expect_error(altExps(es, subsetName = bad) <- list(),
                 msg, info = info)
    expect_error(altExpNames(es, subsetName = bad), msg,
                 info = info)
    expect_error(altExpNames(es, subsetName = bad) <- character(0),
                 msg, info = info)
    expect_error(reducedDim(es, subsetName = bad), msg, info = info)
    expect_error(reducedDim(es, "PCA", subsetName = bad) <- pca,
                 msg, info = info)
    expect_error(reducedDims(es, subsetName = bad), msg,
                 info = info)
    expect_error(reducedDims(es, subsetName = bad) <- list(),
                 msg, info = info)
    expect_error(reducedDimNames(es, subsetName = bad), msg,
                 info = info)
    expect_error(reducedDimNames(es, subsetName = bad) <- "X",
                 msg, info = info)
    expect_error(rowData(es, subsetName = bad), msg, info = info)
    expect_error(rowData(es, subsetName = bad) <- DataFrame(),
                 msg, info = info)
    expect_error(colData(es, subsetName = bad), msg, info = info)
    expect_error(colData(es, subsetName = bad) <- DataFrame(),
                 msg, info = info)
    expect_error(spatialCoords(es, subsetName = bad), msg,
                 info = info)
    expect_error(spatialData(es, subsetName = bad), msg,
                 info = info)
  }
  # NULL is subsetAssayNames()'s default, which lists every subset
  expect_true(all(c("s1", "s2") %in%
                    subsetAssayNames(es, subsetName = NULL)))
  for (bad in badValues[-1]) {
    info <- deparse(bad)
    expect_error(subsetAssayNames(es, subsetName = bad), msg,
                 info = info)
  }
  # (#97) these dispatched on subsetName = "character", so NULL, NA or a
  # number failed with "unable to find an inherited method"
  for (bad in c(badValues, list(1L))) {
    info <- deparse(bad)
    expect_error(subsetDim(es, bad), msg, info = info)
    expect_error(getSubsetAssay(es, bad), msg, info = info)
    expect_error(getSubsetAssay(es, bad, subsetAssayName = "x"),
                 msg, info = info)
    expect_error(subsetColData(es, bad) <- DataFrame(),
                 msg, info = info)
    expect_error(subsetRowData(es, bad) <- DataFrame(),
                 msg, info = info)
    expect_error(
      suppressWarnings(subsetSpatialData(es, bad)), msg, info = info
    )
    expect_error(subsetColnames(es, bad), msg, info = info)
    expect_error(subsetColnames(es, bad) <- letters[1:3],
                 msg, info = info)
    expect_error(subsetRownames(es, bad), msg, info = info)
    expect_error(subsetRownames(es, bad) <- letters[1:5],
                 msg, info = info)
    expect_error(subsetColData(es, bad), msg, info = info)
    expect_error(subsetRowData(es, bad), msg, info = info)
    expect_error(getSubsetExperiment(es, bad), msg, info = info)
    expect_error(setSubsetExperiment(es, bad, es), msg,
                 info = info)
    expect_error(subsetSpatialCoords(es, bad), msg, info = info)
    expect_error(setSubsetAssay(es, bad, pca, "a"), msg,
                 info = info)
  }
  # subsetParent() dispatches on x only, so it gets every bad value
  for (bad in badValues) {
    expect_error(subsetParent(es, bad), msg, info = deparse(bad))
  }
  # the links of a TreeSummarizedExperiment subset (#97)
  tes <- ExperimentSubset(makeReducedDimObjects()$tse)
  tes <- createSubset(tes, "s1", rows = 1:5, cols = 2:4)
  for (bad in c(badValues, list(1L))) {
    info <- deparse(bad)
    expect_error(subsetRowLinks(tes, bad), msg, info = info)
    expect_error(subsetColLinks(tes, bad), msg, info = info)
  }
  # on an object without subsets, that comes first, whatever subsetName is
  sce <- makeReducedDimObjects()$sce
  expect_error(reducedDimNames(sce, subsetName = NULL),
               "which has no subsets")
  expect_error(reducedDimNames(sce, subsetName = NULL) <- "X",
               "which has no subsets")
  expect_error(reducedDimNames(sce, NA_character_),
               "which has no subsets")
})

test_that("An i that isn't one string gives a clear error (#85)", {
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
      expect_error(assay(es, bad), iMsg, info = info)
      expect_error(assay(es, bad) <- m, iMsg, info = info)
      expect_error(getSubsetAssay(es, bad), subsetNameMsg,
                   info = info)
    }
    # an NA subsetAssayName became "s1/NA" and said that isn't an assay
    expect_error(
      assay(es, "s1", subsetAssayName = NA_character_),
      "'subsetAssayName' must be a single character string", info = name
    )
    # a numeric i still goes to the main object's assays
    expect_equal(assay(es, 1), assay(es, "counts"), info = name)
  }
})

test_that("ExperimentSubset() keeps an ExperimentSubset (#89)", {
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
    expect_identical(ExperimentSubset(es), es, info = name)
  }

  # subset = still creates a subset on it
  es <- ExperimentSubset(inputs$sce)
  withSubset <- ExperimentSubset(es, subset = list(subsetName = "s1",
                                                   rows = 1:3,
                                                   cols = 1:4))
  expect_equal(subsetNames(withSubset), "s1")
  expect_equal(subsetDim(withSubset, subsetName = "s1"), c(3, 4))
})

test_that("Subset* subclasses store the base class (#89)", {
  setClass("Issue89SubsetSce", contains = "SubsetSingleCellExperiment")
  on.exit(removeClass("Issue89SubsetSce"), add = TRUE)
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  es <- ExperimentSubset(list(counts = counts))
  mine <- new("Issue89SubsetSce", es)

  # ExperimentSubset() returns it as is
  expect_identical(ExperimentSubset(mine), mine)

  # createSubset() stores a SingleCellExperiment, not the subclass
  mine <- createSubset(mine, "s1", rows = 1:3, cols = 1:4)
  stored <- getSubsetExperiment(mine, "s1")
  expect_equal(as.character(class(stored)), "SingleCellExperiment")

  # cbind() keeps the subsets, still stored as SingleCellExperiment
  other <- es
  colnames(other) <- paste0("other", 1:10)
  other <- createSubset(other, "s2", cols = 1:2)
  other <- new("Issue89SubsetSce", other)
  combined <- cbind(mine, other)
  expect_equal(subsetNames(combined), c("s1", "s2"))
  for (name in c("s1", "s2")) {
    stored <- getSubsetExperiment(combined, name)
    expect_equal(as.character(class(stored)),
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
      expect_equal(as.character(class(stored)),
                   "SingleCellExperiment",
                   info = paste(how, name))
    }
  }
  expect_equal(subsetNames(combinations$rbind), c("s1", "s3"))
  expect_equal(subsetNames(combinations$combineCols),
               c("s1", "s2"))
})

test_that("ExperimentSubset() converts subclasses (#89)", {
  setClass("Issue89Sce", contains = "SingleCellExperiment",
           slots = c(extra = "character"))
  on.exit(removeClass("Issue89Sce"), add = TRUE)
  counts <- matrix(seq_len(60), nrow = 6,
                   dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
  sce <- SingleCellExperiment(assays = list(counts = counts))
  mine <- new("Issue89Sce", sce, extra = "dropped")

  # converted as a SingleCellExperiment, with a warning
  expect_warning(
    es <- ExperimentSubset(mine),
    regexp = "'x' is an object of class Issue89Sce, a subclass of"
  )
  expect_equal(as.character(class(es)),
               "SubsetSingleCellExperiment")
  expect_equal(assay(es, "counts"), counts)

  # subset = still creates a subset
  subset <- list(subsetName = "s1", rows = 1:3)
  expect_warning(es <- ExperimentSubset(mine, subset = subset),
                 regexp = "Issue89Sce")
  expect_equal(subsetNames(es), "s1")
})

test_that("ExperimentSubset() rejects an unsupported x (#89)", {
  counts <- matrix(seq_len(60), nrow = 6)
  expect_error(ExperimentSubset(counts),
               regexp = "'x' must be a SummarizedExperiment")
})

test_that("unlink can name the assays to unlink (#59)", {
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
  messages <- capture_messages(expect_error(
    setSubsetExperiment(es, "s1", changedCounts, unlink = "decontXcounts"),
    regexp = "unlink = c(\"decontXcounts\", \"counts\")", fixed = TRUE
  ))
  expect_length(messages, 0)
  # only the named assay is unlinked, with a message saying so
  messages <- capture_messages(
    es2 <- setSubsetExperiment(es, "s1", sub, unlink = "decontXcounts")
  )
  expect_length(messages, 1)
  expect_match(messages, "'s1' now stores its own 'decontXcounts'")
  expect_equal(subsetAssayNames(es2, "s1"),
               c("counts", "decontXcounts"))
  # s1 stores only decontXcounts; counts stays linked
  expect_equal(subsetAssayNames(es2),
               c("s1", "s1/decontXcounts"))
  expect_equal(as.matrix(assay(es2, "s1/decontXcounts")),
               ownDecontX)
  expect_equal(as.matrix(assay(es2, "s1/counts")),
               simCounts[1:5, 1:8])
  expect_equal(assay(es2, "decontXcounts"), simCounts * 0.5)
  # names that are not linked or not changed are ignored
  messages <- capture_messages(
    es3 <- setSubsetExperiment(es, "s1", sub,
                               unlink = c("decontXcounts", "counts", "other"))
  )
  expect_length(messages, 1)
  expect_equal(subsetAssayNames(es3),
               c("s1", "s1/decontXcounts"))
  expect_equal(as.matrix(assay(es3, "s1/counts")),
               simCounts[1:5, 1:8])
  expect_equal(as.matrix(assay(es3, "s1/decontXcounts")),
               ownDecontX)
  # character(0) unlinks nothing, as FALSE
  expect_error(
    setSubsetExperiment(es, "s1", sub, unlink = character(0)),
    regexp = "'decontXcounts' in subset 's1' is linked"
  )
  # every changed linked assay unlink doesn't allow is named at once
  assay(changedCounts, "decontXcounts") <- ownDecontX + 1
  expect_error(
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
  messages <- capture_messages(expect_error(
    setSubsetExperiment(es, "s1", badName, unlink = TRUE),
    regexp = "can't contain '/'"
  ))
  expect_length(messages, 0)
})

test_that("unlink = TRUE says which assays it unlinks (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  assay(es, "decontXcounts") <- simCounts * 0.5
  sub <- getSubsetExperiment(es, "s1")
  assay(sub, "decontXcounts") <- as.matrix(assay(sub, "counts")) * 0.9
  assay(sub, "counts") <- as.matrix(assay(sub, "counts")) + 1L
  messages <- capture_messages(
    es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE)
  )
  countsMessage <- "'s1' now stores its own 'counts', unlinked from the main"
  expect_true(any(grepl(countsMessage, messages)))
  decontXMessage <- "'s1' now stores its own 'decontXcounts'"
  expect_true(any(grepl(decontXMessage, messages)))
  # storing an own assay again unlinks nothing, so says nothing
  expect_silent(
    es <- setSubsetExperiment(es, "s1", sub, unlink = TRUE)
  )
  # a subset made from s1 names s1's assay as the source
  es <- createSubset(es, "child", rows = 1:2, cols = 1:2, parentAssay = "s1")
  expect_message(
    es <- setSubsetAssay(es, "child", matrix(0, 2, 2), "counts",
                         unlink = TRUE),
    regexp = "'child' now stores its own 'counts', unlinked from 's1/counts'"
  )
})

test_that("setSubsetAssay and assay<- take unlink names (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  own <- simCounts[1:5, 1:8] + 1
  expect_error(
    setSubsetAssay(es, "s1", own, "counts", unlink = "decontXcounts"),
    regexp = "unlink = c(\"decontXcounts\", \"counts\")", fixed = TRUE
  )
  expect_message(
    es2 <- setSubsetAssay(es, "s1", own, "counts", unlink = "counts"),
    regexp = "'s1' now stores its own 'counts'"
  )
  expect_equal(as.matrix(assay(es2, "s1/counts")), own)
  expect_error(
    assay(es, "s1/counts", unlink = "decontXcounts") <- own,
    regexp = "is linked"
  )
  expect_message(assay(es, "s1/counts", unlink = "counts") <- own,
                 regexp = "now stores its own 'counts'")
  expect_equal(as.matrix(assay(es, "s1/counts")), own)
  # on a main assay, unlink names are an error, as TRUE is
  expect_error(
    assay(es, "counts", unlink = "counts") <- simCounts,
    regexp = "only applies"
  )
  # character(0) unlinks nothing, so is allowed there, as FALSE is
  assay(es, "counts", unlink = character(0)) <- simCounts + 1
  expect_equal(assay(es, "counts"), simCounts + 1)
})

test_that("An invalid unlink is an error (#59)", {
  es <- ExperimentSubset(simSce)
  es <- createSubset(es, "s1", rows = 1:5, cols = 1:8)
  sub <- getSubsetExperiment(es, "s1")
  own <- simCounts[1:5, 1:8] + 1
  msg <- "'unlink' must be TRUE, FALSE or a character vector of assay names"
  for (bad in list(NA, c(TRUE, FALSE), 1, NA_character_, c("counts", ""))) {
    expect_error(setSubsetExperiment(es, "s1", sub, unlink = bad),
                 regexp = msg, fixed = TRUE)
    expect_error(setSubsetAssay(es, "s1", own, "x", unlink = bad),
                 regexp = msg, fixed = TRUE)
    expect_error(assay(es, "s1/x", unlink = bad) <- own,
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

test_that("A subset sees the main object's altExps (#58)", {
  es <- makeAdtFixture()
  expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(dim(adt), c(2L, 8L))
  expect_equal(colnames(adt), paste0("c", 1:8))
  expect_equal(assay(adt, "counts"), adtCounts[, 1:8])
  # by position and in altExps(), in the same order
  expect_equal(
    assay(altExp(es, 1, subsetName = "s1"), "counts"), adtCounts[, 1:8])
  expect_equal(
    assay(altExp(es, subsetName = "s1"), "counts"), adtCounts[, 1:8])
  expect_equal(names(altExps(es, subsetName = "s1")), "ADT")
  # and in getSubsetExperiment()
  sub <- getSubsetExperiment(es, "s1")
  expect_equal(altExpNames(sub), "ADT")
  expect_equal(assay(altExp(sub, "ADT"), "counts"),
               adtCounts[, 1:8])
  # a view: a later change to the main object's altExp shows in the subset
  altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts * 2L))
  expect_equal(
    assay(altExp(es, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] * 2L)
  # the subset's own altExps come after the linked ones
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es, "spikes", subsetName = "s1") <- spikes
  expect_equal(altExpNames(es, subsetName = "s1"),
               c("ADT", "spikes"))
  expect_equal(nrow(altExp(es, 2, subsetName = "s1")), 3L)
  # an altExp the subset doesn't see gives SingleCellExperiment's error
  expect_error(altExp(es, "none", subsetName = "s1"),
               regexp = "'none' not in 'altExpNames")
  expect_error(altExp(es, 3, subsetName = "s1"))
  none <- createSubset(ExperimentSubset(simSce), "s", cols = 1:3)
  expect_error(altExp(none, subsetName = "s"),
               regexp = "no available entries")
})

test_that("withColData and withDimnames work on a subset's altExps (#58)", {
  es <- makeAdtFixture()
  colData(es, subsetName = "s1") <- S4Vectors::DataFrame(group = letters[1:8])
  adt <- altExp(es, "ADT", subsetName = "s1", withColData = TRUE)
  expect_equal(adt$group, letters[1:8])
  # the subset's colData columns are dropped again when stored
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1", withColData = TRUE) <- adt
  expect_equal(ncol(colData(ownAltExp(es, "s1", "ADT"))), 0L)
  # a value with other colnames is accepted only without withDimnames
  renamed <- adt
  colnames(renamed) <- paste0("x", 1:8)
  expect_warning(altExp(es, "ADT", subsetName = "s1") <- renamed,
                 regexp = "are not the same as")
  altExp(es, "ADT", subsetName = "s1", withDimnames = FALSE) <- renamed
  expect_equal(colnames(altExp(es, "ADT", subsetName = "s1")),
               paste0("c", 1:8))
})

test_that("A subset stores only the altExp assays it adds (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  messages <- capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- adt)
  expect_length(messages, 0)
  expect_equal(
    assayNames(altExp(es, "ADT", subsetName = "s1")),
    c("counts", "decontPro"))
  expect_equal(assayNames(altExp(es, "ADT")), "counts")
  expect_equal(assayNames(ownAltExp(es, "s1", "ADT")),
               "decontPro")
  # counts stays linked
  altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts + 1L))
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(assay(adt, "counts"), adtCounts[, 1:8] + 1L)
  expect_equal(assay(adt, "decontPro"), adtCounts[, 1:8] / 2)
  # storing what the subset sees stores nothing
  es <- makeAdtFixture()
  altExp(es, "ADT", subsetName = "s1") <- altExp(es, "ADT", subsetName = "s1")
  expect_null(ownAltExp(es, "s1", "ADT"))
})

test_that("A changed linked altExp assay needs unlink (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "counts") <- assay(adt, "counts") + 100L
  messages <- capture_messages(expect_error(
    altExp(es, "ADT", subsetName = "s1") <- adt,
    regexp = paste0("'ADT/counts' in subset 's1' is linked to the main ",
                    "altExp 'ADT'. Set unlink = TRUE or unlink = ",
                    "\"ADT/counts\" to store a separate 'ADT/counts'"),
    fixed = TRUE
  ))
  expect_length(messages, 0)
  expect_null(ownAltExp(es, "s1", "ADT"))
  # a bare name is a main assay, not the altExp's
  expect_error(
    altExp(es, "ADT", subsetName = "s1", unlink = "counts") <- adt,
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  for (unlink in list(TRUE, "ADT/counts")) {
    es2 <- es
    messages <- capture_messages(
      altExp(es2, "ADT", subsetName = "s1", unlink = unlink) <- adt)
    expect_equal(
      messages,
      paste0("Subset 's1' now stores its own 'ADT/counts', unlinked from ",
             "the main altExp 'ADT'.\n"))
    expect_equal(
      assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
      adtCounts[, 1:8] + 100L)
    expect_equal(assay(altExp(es2, "ADT"), "counts"), adtCounts)
  }
  # also through altExps<-
  expect_error(
    altExps(es, subsetName = "s1") <- list(ADT = adt),
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  # unlink only applies in a subset
  expect_error(
    altExp(es, "ADT", unlink = TRUE) <- altExp(es, "ADT"),
    regexp = "'unlink' only applies when storing in a subset.", fixed = TRUE)
  expect_error(
    altExp(es, "ADT", subsetName = "s1", unlink = NA) <- adt,
    regexp = "'unlink' must be TRUE, FALSE", fixed = TRUE)
})

test_that("An altExp of another shape replaces the linked one (#58)", {
  es <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  expect_error(
    altExp(es, "ADT", subsetName = "s1") <- small,
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  suppressMessages(
    altExp(es, "ADT", subsetName = "s1", unlink = "ADT/counts") <- small)
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(rownames(adt), "p9")
  expect_equal(assayNames(adt), "counts")
  expect_equal(assay(altExp(es, "ADT"), "counts"), adtCounts)
})

test_that("A subset's altExp rowData and metadata replace the linked ones (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  rowData(adt)$target <- c("CD3", "CD4")
  metadata(adt)$decontPro <- list(iterations = 10)
  altExp(es, "ADT", subsetName = "s1") <- adt
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(rowData(adt)$target, c("CD3", "CD4"))
  expect_equal(metadata(adt)$decontPro, list(iterations = 10))
  expect_equal(ncol(rowData(altExp(es, "ADT"))), 0L)
  own <- ownAltExp(es, "s1", "ADT")
  expect_length(assayNames(own), 0)
  expect_equal(colnames(rowData(own)), "target")
})

test_that("A chained subset sees altExps through its parent (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  spikes <- SummarizedExperiment(
    list(counts = matrix(seq_len(24), 3, 8,
                         dimnames = list(NULL, paste0("c", 1:8)))))
  altExp(es, "spikes", subsetName = "s1") <- spikes
  es <- createSubset(es, "s2", cols = 3:5, parentAssay = "s1")
  expect_equal(altExpNames(es, subsetName = "s2"),
               c("ADT", "spikes"))
  view <- altExp(es, "ADT", subsetName = "s2")
  expect_equal(colnames(view), paste0("c", 3:5))
  expect_equal(assay(view, "counts"), adtCounts[, 3:5])
  expect_equal(assay(view, "decontPro"), adtCounts[, 3:5] / 2)
  expect_equal(
    as.vector(assay(altExp(es, "spikes", subsetName = "s2"), "counts")),
    as.vector(matrix(seq_len(24), 3, 8)[, 3:5]))
  # a change to an assay of the parent subset names where it is stored
  assay(view, "decontPro") <- assay(view, "decontPro") + 1
  expect_error(
    altExp(es, "ADT", subsetName = "s2") <- view,
    regexp = paste0("'ADT/decontPro' in subset 's2' is linked to altExp ",
                    "'ADT' of subset 's1'"),
    fixed = TRUE)
  messages <- capture_messages(
    altExp(es, "ADT", subsetName = "s2", unlink = TRUE) <- view)
  expect_equal(
    messages,
    paste0("Subset 's2' now stores its own 'ADT/decontPro', unlinked from ",
           "altExp 'ADT' of subset 's1'.\n"))
})

test_that("A subset removes and renames only its own altExps (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  # NULL removes the subset's part; the linked altExp shows again
  altExp(es, "ADT", subsetName = "s1") <- NULL
  expect_null(ownAltExp(es, "s1", "ADT"))
  expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
               "counts")
  # altExps<- can't remove a linked altExp
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  altExps(es, subsetName = "s1") <- list(spikes = spikes)
  expect_equal(altExpNames(es, subsetName = "s1"),
               c("ADT", "spikes"))
  altExps(es, subsetName = "s1") <- list()
  expect_equal(altExpNames(es, subsetName = "s1"), "ADT")
  # renaming: the subset's own altExps only
  altExp(es, "spikes", subsetName = "s1") <- spikes
  altExpNames(es, subsetName = "s1") <- c("ADT", "ercc")
  expect_equal(altExpNames(es, subsetName = "s1"),
               c("ADT", "ercc"))
  expect_equal(altExpNames(es), "ADT")
  expect_error(
    altExpNames(es, subsetName = "s1") <- c("protein", "ercc"),
    regexp = paste0("Subset 's1' can't rename 'ADT', which is linked to ",
                    "the main altExp 'ADT'."),
    fixed = TRUE)
  # names must differ, so an own altExp can't take a linked name, which
  # would hide the linked altExp without unlink
  expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", "ADT"),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", NA),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  expect_error(
    altExpNames(es, subsetName = "s1") <- c("ADT", "ercc", "x"),
    regexp = "one different name for each of the 2 altExps", fixed = TRUE)
  expect_equal(altExpNames(es, subsetName = "s1"),
               c("ADT", "ercc"))
})

test_that("altExpNames<- rejects values that aren't character", {
  # (#105) other values failed at S4 dispatch, before the package's check
  valueMsg <- "'value' must be a character vector of altExp names."
  classes <- c("SingleCellExperiment", "TreeSummarizedExperiment",
               "SpatialExperiment")
  for (class in classes) {
    es <- makeAdtFixture(class)
    # a factor of the right length, NULL and a number
    expect_error(
      altExpNames(es, subsetName = "s1") <- factor("protein"),
      valueMsg, fixed = TRUE, info = class
    )
    expect_error(altExpNames(es, subsetName = "s1") <- NULL,
                 valueMsg, fixed = TRUE, info = class)
    expect_error(altExpNames(es, subsetName = "s1") <- 1L,
                 valueMsg, fixed = TRUE, info = class)
    expect_error(altExpNames(es) <- factor("protein"),
                 valueMsg, fixed = TRUE, info = class)
    expect_error(altExpNames(es) <- NULL,
                 valueMsg, fixed = TRUE, info = class)
    expect_error(altExpNames(es) <- 1L,
                 valueMsg, fixed = TRUE, info = class)
    # nothing was renamed
    expect_equal(altExpNames(es), "ADT", info = class)
    expect_equal(altExpNames(es, subsetName = "s1"), "ADT",
                 info = class)
    # a character value still renames the main altExp
    altExpNames(es) <- "protein"
    expect_equal(altExpNames(es), "protein", info = class)
    expect_equal(altExpNames(es, subsetName = "s1"), "protein",
                 info = class)
  }
})

test_that("Writes after an unlink don't need unlink again (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "counts") <- assay(adt, "counts") + 100L
  suppressMessages(altExp(es, "ADT", subsetName = "s1", unlink = TRUE) <- adt)
  # s1's own counts is no longer linked
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  messages <- capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- adt)
  expect_length(messages, 0)
  expect_equal(assayNames(ownAltExp(es, "s1", "ADT")),
               c("counts", "decontPro"))
  spikes <- SummarizedExperiment(
    list(counts = matrix(1L, 3, 8, dimnames = list(NULL, paste0("c", 1:8)))))
  messages <- capture_messages(
    altExp(es, "spikes", subsetName = "s1") <- spikes)
  expect_length(messages, 0)
  messages <- capture_messages(
    es2 <- setSubsetExperiment(es, "s1", getSubsetExperiment(es, "s1")))
  expect_length(messages, 0)
  expect_equal(
    assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] + 100L)
  # also after storing an altExp of another shape
  es <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  suppressMessages(altExp(es, "ADT", subsetName = "s1", unlink = TRUE) <- small)
  assay(small, "counts") <- assay(small, "counts") + 1L
  messages <- capture_messages(
    altExp(es, "ADT", subsetName = "s1") <- small)
  expect_length(messages, 0)
  expect_equal(
    as.vector(assay(altExp(es, "ADT", subsetName = "s1"), "counts")),
    rep(2L, 8))
})

test_that("A subset's altExp reducedDims replace the linked ones (#58)", {
  es <- makeAdtFixture()
  adtSce <- SingleCellExperiment(list(counts = adtCounts))
  reducedDim(adtSce, "PCA") <- matrix(seq_len(40), 20, 2)
  altExp(es, "ADT") <- adtSce
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(unname(reducedDim(adt, "PCA")),
               matrix(seq_len(40), 20, 2)[1:8, ])
  # storing the same reducedDims stores nothing
  altExp(es, "ADT", subsetName = "s1") <- adt
  expect_null(ownAltExp(es, "s1", "ADT"))
  reducedDim(adt, "PCA") <- matrix(0, 8, 2)
  metadata(adt)$note <- "own"
  altExp(es, "ADT", subsetName = "s1") <- adt
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(unname(reducedDim(adt, "PCA")), matrix(0, 8, 2))
  expect_equal(metadata(adt)$note, "own")
  expect_equal(unname(reducedDim(altExp(es, "ADT"), "PCA")),
               matrix(seq_len(40), 20, 2))
  expect_length(assayNames(ownAltExp(es, "s1", "ADT")), 0)
  # an altExp without reducedDims shows the linked ones and stores nothing
  es <- makeAdtFixture()
  altExp(es, "ADT") <- adtSce
  adt <- altExp(es, "ADT", subsetName = "s1")
  reducedDims(adt) <- list()
  altExp(es, "ADT", subsetName = "s1") <- adt
  expect_null(ownAltExp(es, "s1", "ADT"))
})

test_that("A rename carries what subsets store along (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  es <- createSubset(es, "s2", cols = 1:3, parentAssay = "s1")
  adt <- altExp(es, "ADT", subsetName = "s2")
  assay(adt, "zero") <- assay(adt, "counts") * 0L
  altExp(es, "ADT", subsetName = "s2") <- adt
  altExpNames(es) <- "protein"
  expect_equal(altExpNames(es, subsetName = "s1"), "protein")
  expect_equal(
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
  expect_equal(altExpNames(es, subsetName = "s2"),
               c("protein", "ercc"))
  expect_equal(assayNames(altExp(es, "ercc", subsetName = "s2")),
               c("counts", "logcounts"))
  # a rename onto a name a child already stores is an error
  other <- SummarizedExperiment(
    list(counts = matrix(1L, 1, 3, dimnames = list(NULL, paste0("c", 1:3)))))
  altExp(es, "other", subsetName = "s2") <- other
  expect_error(
    altExpNames(es, subsetName = "s1") <- c("protein", "other"),
    regexp = "Subset 's2' already stores an altExp called 'other'.",
    fixed = TRUE)
  # also when the subset stores only the new name, an altExp of its own
  es <- makeAdtFixture()
  own <- SummarizedExperiment(
    list(counts = matrix(1L, 2, 8, dimnames = list(c("p1", "p2"),
                                                   paste0("c", 1:8)))))
  altExp(es, "X", subsetName = "s1") <- own
  expect_error(
    altExpNames(es) <- "X",
    regexp = "Subset 's1' already stores an altExp called 'X'.",
    fixed = TRUE)
  expect_equal(altExpNames(es), "ADT")
})

test_that("A subset's whole altExp stays whole (#58)", {
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
  expect_warning(
    altExp(es, "ADT") <- twoRows(1L, c("counts", "other")), regexp = NA)
  adt <- altExp(es, "ADT", subsetName = "s1")
  expect_equal(assayNames(adt), "counts")
  expect_equal(as.vector(assay(adt, "counts")), rep(5L, 16))
  # and is kept when the main ADT is reshaped
  expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts)),
    regexp = NA)
  expect_equal(rownames(altExp(es, "ADT", subsetName = "s1")),
               c("r1", "r2"))
  # and stays whole through later writes of other altExps and round trips
  altExp(es, "ADT") <- twoRows(1L, c("counts", "other"))
  altExp(es, "B", subsetName = "s1") <- own
  es <- setSubsetExperiment(es, "s1", getSubsetExperiment(es, "s1"))
  expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
               "counts")
  expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts)),
    regexp = NA)
  expect_equal(rownames(altExp(es, "ADT", subsetName = "s1")),
               c("r1", "r2"))
  # an altExp the main object adds under a name s1 already stores is
  # hidden by s1's
  altExp(es, "X", subsetName = "s1") <- own
  altExp(es, "X") <- twoRows(1L, c("counts", "other"))
  expect_equal(assayNames(altExp(es, "X", subsetName = "s1")),
               "counts")
})

test_that("Objects saved before pieces were recorded read as before (#58)", {
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
  expect_equal(assayNames(adt), c("counts", "other"))
  expect_equal(assay(adt, "counts"), adtCounts[, 1:8] + 5L)
})

test_that("combineCols carries what subsets store along (#58)", {
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
  expect_warning(
    out <- combineCols(es, makeOther(c("p1", "p2"))), regexp = NA)
  expect_equal(assayNames(altExp(out, "ADT", subsetName = "s1")),
               c("counts", "decontPro"))
  # a third protein: ADT gets another row, and s1's part is dropped
  expect_warning(
    out <- combineCols(es, makeOther(c("p1", "p2", "p3"))),
    regexp = "Dropped what subset 's1' stored in altExp 'ADT'", fixed = TRUE)
  expect_equal(dim(altExp(out, "ADT", subsetName = "s1")),
               c(3L, 8L))
  expect_equal(assayNames(altExp(out, "ADT", subsetName = "s1")),
               "counts")
})

test_that("Removing or reshaping an altExp drops what subsets added to it (#58)", {
  es <- makeAdtFixture()
  adt <- altExp(es, "ADT", subsetName = "s1")
  assay(adt, "decontPro") <- assay(adt, "counts") / 2
  altExp(es, "ADT", subsetName = "s1") <- adt
  # new values keep the subset's part
  expect_warning(
    altExp(es, "ADT") <- SummarizedExperiment(list(counts = adtCounts + 1L)),
    regexp = NA)
  expect_equal(assayNames(altExp(es, "ADT", subsetName = "s1")),
               c("counts", "decontPro"))
  dropped <- paste0("Dropped what subset 's1' stored in altExp 'ADT': the ",
                    "altExp it added to was removed, or now has other rows ",
                    "or another class.")
  # other rows
  es2 <- es
  expect_warning(altExp(es2, "ADT") <- altExp(es2, "ADT")[1, ],
                 regexp = dropped, fixed = TRUE)
  expect_null(ownAltExp(es2, "s1", "ADT"))
  expect_equal(dim(altExp(es2, "ADT", subsetName = "s1")),
               c(1L, 8L))
  # removed, one or all
  es2 <- es
  expect_warning(altExp(es2, "ADT") <- NULL,
                 regexp = dropped, fixed = TRUE)
  expect_length(altExpNames(es2, subsetName = "s1"), 0)
  es2 <- es
  expect_warning(altExps(es2) <- list(),
                 regexp = dropped, fixed = TRUE)
  expect_length(altExpNames(es2, subsetName = "s1"), 0)
  # a whole altExp of the subset's own is kept
  es3 <- makeAdtFixture()
  small <- SummarizedExperiment(list(
    counts = matrix(1L, 1, 8, dimnames = list("p9", paste0("c", 1:8)))))
  suppressMessages(altExp(es3, "ADT", subsetName = "s1", unlink = TRUE) <- small)
  expect_warning(altExp(es3, "ADT") <- NULL, regexp = NA)
  expect_equal(rownames(altExp(es3, "ADT", subsetName = "s1")),
               "p9")
  # a parent subset's change reaches its children, also through
  # setSubsetExperiment()
  es4 <- makeAdtFixture()
  es4 <- createSubset(es4, "s2", cols = 1:3, parentAssay = "s1")
  adt <- altExp(es4, "ADT", subsetName = "s2")
  assay(adt, "zero") <- assay(adt, "counts") * 0L
  altExp(es4, "ADT", subsetName = "s2") <- adt
  expect_warning(
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
  expect_warning(
    es5 <- setSubsetExperiment(es5, "s1", sub),
    regexp = "Dropped what subset 's2' stored in altExp 'spikes'",
    fixed = TRUE)
  expect_equal(altExpNames(es5, subsetName = "s2"), "ADT")
})

test_that("setSubsetExperiment stores only changed altExp assays (#58)", {
  es <- makeAdtFixture()
  sub <- getSubsetExperiment(es, "s1")
  # a round trip stores no copy of ADT
  es2 <- setSubsetExperiment(es, "s1", sub)
  expect_null(ownAltExp(es2, "s1", "ADT"))
  expect_equal(altExpNames(es2, subsetName = "s1"), "ADT")
  # a new assay in ADT is stored without unlink
  assay(altExp(sub, "ADT"), "decontPro") <- adtCounts[, 1:8] / 2
  es2 <- setSubsetExperiment(es, "s1", sub)
  expect_equal(assayNames(ownAltExp(es2, "s1", "ADT")),
               "decontPro")
  # changed linked assays, main and altExp, are named in one error
  changed <- sub
  assay(changed, "counts") <- as.matrix(assay(changed, "counts")) + 1
  assay(altExp(changed, "ADT"), "counts") <- adtCounts[, 1:8] + 1L
  messages <- capture_messages(expect_error(
    setSubsetExperiment(es, "s1", changed),
    regexp = paste0("'counts' in subset 's1' is linked to the main assay ",
                    "'counts' and 'ADT/counts' to the main altExp 'ADT'. ",
                    "Set unlink = TRUE or unlink = c(\"counts\", ",
                    "\"ADT/counts\")"),
    fixed = TRUE))
  expect_length(messages, 0)
  expect_error(
    setSubsetExperiment(es, "s1", changed, unlink = "counts"),
    regexp = "'ADT/counts' in subset 's1' is linked", fixed = TRUE)
  messages <- capture_messages(
    es2 <- setSubsetExperiment(es, "s1", changed, unlink = TRUE))
  expect_length(messages, 2)
  expect_match(messages[2], "its own 'ADT/counts'", fixed = TRUE)
  expect_equal(
    assay(altExp(es2, "ADT", subsetName = "s1"), "counts"),
    adtCounts[, 1:8] + 1L)
})

test_that("TreeSE and SpatialExperiment subsets see altExps (#58)", {
  for (class in c("TreeSummarizedExperiment", "SpatialExperiment")) {
    es <- makeAdtFixture(class)
    expect_equal(altExpNames(es, subsetName = "s1"), "ADT",
                 info = class)
    expect_equal(
      assay(altExp(es, "ADT", subsetName = "s1"), "counts"),
      adtCounts[, 1:8], info = class)
    sub <- getSubsetExperiment(es, "s1")
    expect_equal(altExpNames(sub), "ADT", info = class)
  }
})
