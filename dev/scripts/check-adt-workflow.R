# RNA and ADT workflows on one ExperimentSubset object (issue #58)
#
# Real-data check of subsets that see the main object's altExps, on the
# 10x 5k PBMC CITE-seq data set (RNA + 32 antibodies) from
# DropletTestFiles. Runs the RNA workflow of the PBMC article (cell calls,
# QC, decontX, clustering, SingleR, T-cell re-analysis) and an ADT
# workflow (decontPro, CLR) on the same object, and prints PASS/FAIL for
# each check.
#
# Not a test: it downloads about 120 MB (cached by BiocFileCache) and takes
# about 8 minutes, most of it in decontX and decontPro. Run it from the
# package folder, which it loads with devtools::load_all():
#
#   Rscript dev/scripts/check-adt-workflow.R
#
# Needs DropletTestFiles, DropletUtils, scuttle, scDblFinder, celda,
# decontX, scrapper, SingleR and celldex. On 2026-11-05 (1.23.37) all 33
# checks passed. The workflow article (vignettes/articles/pbmc-workflow.Rmd)
# now shows these steps (#93), except the last section's carry-along
# checks, which wait for subsets made from an altExp (#92).

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) > 0) args[1] else "."
options(warn = 1)
t0 <- Sys.time()

suppressPackageStartupMessages({
  devtools::load_all(repo, quiet = TRUE)
  library(DropletTestFiles)
  library(DropletUtils)
  library(scuttle)
  library(scDblFinder)
  library(celda)
  library(decontX)
  library(scrapper)
  library(SingleR)
  library(celldex)
})

results <- data.frame(check = character(0), status = character(0))
check <- function(what, ok) {
  status <- if (isTRUE(ok)) "PASS" else "FAIL"
  results[nrow(results) + 1, ] <<- c(what, status)
  message(status, ": ", what)
}
step <- function(what) {
  message("\n== ", what, " (", format(round(Sys.time() - t0, 1)), ")")
}
mb <- function(x) round(as.numeric(object.size(x)) / 2^20, 1)
sameValues <- function(a, b) {
  isTRUE(all.equal(as.matrix(a), as.matrix(b), check.attributes = FALSE))
}
# messages other than the methods package's notes on its class cache,
# which SeuratObject (loaded by decontX) sets off on S4 dispatch
ownMessages <- function(messages) {
  messages[!grepl("^(Found more than one class|Also defined by)", messages)]
}
# what a subset stores itself of an altExp (NULL if nothing)
ownAltExp <- function(es, subsetName, name) {
  internal <- .internalAssay(.subsets(es)[[subsetName]])
  if (!name %in% altExpNames(internal))
    return(NULL)
  altExp(internal, name, withDimnames = FALSE)
}

# ---------------------------------------------------------------------
step("Data")
rawFile <- getTestFile("tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/raw.tar.gz")
filteredFile <- getTestFile(
  "tenx-3.1.0-5k_pbmc_protein_v3/1.0.0/filtered.tar.gz")
dataFolder <- file.path(tempdir(), "pbmc5k")
untar(rawFile, exdir = dataFolder)
untar(filteredFile, exdir = dataFolder)
raw <- read10xCounts(file.path(dataFolder, "raw_feature_bc_matrix"),
                     col.names = TRUE)
filtered <- read10xCounts(file.path(dataFolder,
                                    "filtered_feature_bc_matrix"),
                          col.names = TRUE)
# droplets with any counts; the rest carry no information
raw <- raw[, Matrix::colSums(counts(raw)) > 0]
rownames(raw) <- uniquifyFeatureNames(rowData(raw)$ID, rowData(raw)$Symbol)
# RNA in the main object, antibodies in altExp "ADT"
raw <- splitAltExps(raw, rowData(raw)$Type, ref = "Gene Expression")
altExpNames(raw) <- "ADT"
rownames(altExp(raw, "ADT")) <- sub("_TotalSeqB$", "",
                                    rownames(altExp(raw, "ADT")))
rawAdt <- counts(altExp(raw, "ADT"))
cellBarcodes <- colnames(filtered)
message("droplets: ", ncol(raw), "; cells: ", length(cellBarcodes),
        "; genes: ", nrow(raw), "; antibodies: ", nrow(rawAdt))
check("Cell Ranger's cells are among the kept droplets",
      all(cellBarcodes %in% colnames(raw)))

# ---------------------------------------------------------------------
step("One object: ExperimentSubset with RNA and ADT")
es <- ExperimentSubset(raw)
sizes <- c(start = mb(es))
es <- createSubset(es, subsetName = "filtered", cols = cellBarcodes)
check("filtered sees the main object's ADT",
      identical(altExpNames(es, subsetName = "filtered"), "ADT"))
filteredAdt <- altExp(es, "ADT", subsetName = "filtered")
check("filtered's ADT is 32 antibodies x its cells",
      identical(dim(filteredAdt), c(32L, length(cellBarcodes))))
check("filtered's ADT counts are the raw ADT counts of its cells",
      sameValues(counts(filteredAdt), rawAdt[, cellBarcodes]))
check("filtered's ADT has the subset's cell names",
      identical(colnames(filteredAdt), cellBarcodes))
check("creating a subset copies no ADT", is.null(ownAltExp(es, "filtered", "ADT")))
sizes["filtered subset"] <- mb(es)

# ---------------------------------------------------------------------
step("RNA: cell QC with getSubsetExperiment()/setSubsetExperiment()")
filteredCells <- getSubsetExperiment(es, subsetName = "filtered")
check("getSubsetExperiment() includes ADT",
      identical(altExpNames(filteredCells), "ADT"))
isMito <- grepl("^MT-", rownames(filteredCells))
filteredCells <- addPerCellQCMetrics(filteredCells,
                                     subsets = list(mito = isMito))
set.seed(100)
filteredCells <- scDblFinder(filteredCells)
check("QC used the ADT too (altexps_ADT_sum)",
      "altexps_ADT_sum" %in% colnames(colData(filteredCells)))
es <- setSubsetExperiment(es, subsetName = "filtered", value = filteredCells)
check("the round trip stores no copy of ADT",
      is.null(ownAltExp(es, "filtered", "ADT")))
check("the QC columns are stored",
      all(c("subsets_mito_percent", "scDblFinder.class") %in%
            colnames(colData(es, subsetName = "filtered"))))

# ---------------------------------------------------------------------
step("RNA: decontX")
decontamination <- celda::decontX(assay(es, "filtered"), verbose = FALSE)
es <- setSubsetAssay(es, subsetName = "filtered",
                     inputMatrix = decontamination$decontXcounts,
                     subsetAssayName = "decontaminated")
colData(es, subsetName = "filtered")$contamination <-
  decontamination$contamination

# ---------------------------------------------------------------------
step("RNA: high-quality cells")
qc <- colData(es, subsetName = "filtered")
keep <- qc$detected >= 300 & qc$subsets_mito_percent <= 15 &
  qc$contamination <= 0.3 & qc$scDblFinder.class == "singlet"
message("kept ", sum(keep), " of ", length(keep), " cells")
es <- createSubset(es, subsetName = "clean", cols = which(keep),
                   parentAssay = "filtered/decontaminated")
cleanBarcodes <- subsetColnames(es, "clean")
cleanAdt <- altExp(es, "ADT", subsetName = "clean")
check("clean, a subset of filtered, sees the main ADT for its cells",
      sameValues(counts(cleanAdt), rawAdt[, cleanBarcodes]))

# ---------------------------------------------------------------------
step("RNA: normalise, variable genes, PCA, clusters, UMAP")
cleanCounts <- assay(es, "clean")
sizeFactors <- centerSizeFactors(colSums(cleanCounts))
lognormalized <- normalizeCounts(cleanCounts, sizeFactors, delayed = FALSE)
es <- setSubsetAssay(es, subsetName = "clean", inputMatrix = lognormalized,
                     subsetAssayName = "lognormalized")
variances <- modelGeneVariances(lognormalized)
variableGenes <- chooseHighlyVariableGenes(variances$statistics$residuals,
                                           top = 2000)
es <- createSubset(es, subsetName = "clean.variable", rows = variableGenes,
                   parentAssay = "clean/lognormalized")
check("clean.variable (a subset of rows) still sees all 32 antibodies",
      identical(dim(altExp(es, "ADT", subsetName = "clean.variable")),
                c(32L, length(cleanBarcodes))))
scaled <- t(scale(t(as.matrix(assay(es, "clean.variable")))))
pca <- runPca(scaled)
reducedDim(es, "reduction.pca", subsetName = "clean") <- t(pca$components)
clusters <- clusterGraph(buildSnnGraph(pca$components))
reducedDim(es, "embedding.umap", subsetName = "clean") <-
  runUmap(pca$components)

step("RNA: SingleR labels")
monaco <- MonacoImmuneData()
predictions <- SingleR(test = assay(es, "clean/lognormalized"), ref = monaco,
                       labels = monaco$label.main)
colData(es, subsetName = "clean")$cluster <- clusters$membership
colData(es, subsetName = "clean")$label <- predictions$labels
print(table(predictions$labels))

# ---------------------------------------------------------------------
step("ADT: decontPro on the clean cells")
# ambient ADT from droplets that aren't cells
notCell <- !colnames(es) %in% cellBarcodes
adtTotals <- Matrix::colSums(rawAdt)
ambientPool <- which(notCell & adtTotals > 0)
set.seed(100)
ambientDroplets <- sample(ambientPool, min(20000, length(ambientPool)))
ambientAdt <- as.matrix(rawAdt[, ambientDroplets])
cellType <- as.integer(factor(clusters$membership))
sizeBefore <- mb(es)
mainAdtBefore <- altExp(es, "ADT")
decontProTime <- system.time(
  pro <- decontPro(as.matrix(counts(cleanAdt)), cellType,
                   ambient_counts = ambientAdt)
)
message("decontPro took ", round(decontProTime[["elapsed"]]), " s")
decontaminatedAdt <- pro$decontaminated_counts
dimnames(decontaminatedAdt) <- dimnames(cleanAdt)
adt <- altExp(es, "ADT", subsetName = "clean")
assay(adt, "decontPro") <- decontaminatedAdt
messages <- character(0)
withCallingHandlers(
  altExp(es, "ADT", subsetName = "clean") <- adt,
  message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
messages <- ownMessages(messages)
message("messages: ", paste(messages, collapse = " | "))
check("adding decontPro to clean's ADT needs no unlink and says nothing",
      length(messages) == 0)
check("clean stores only the decontPro assay of ADT",
      identical(assayNames(ownAltExp(es, "clean", "ADT")), "decontPro"))
check("clean's ADT has counts (linked) and decontPro",
      identical(assayNames(altExp(es, "ADT", subsetName = "clean")),
                c("counts", "decontPro")))
check("the main object's ADT is unchanged",
      identical(altExp(es, "ADT"), mainAdtBefore))
sizes["after decontPro"] <- mb(es)
growth <- sizes[["after decontPro"]] - sizeBefore
message("object grew by ", growth, " MB; the decontPro matrix is ",
        mb(decontaminatedAdt), " MB")
check("the object grows by about the decontPro matrix, not a copy of ADT",
      growth < 1.5 * mb(decontaminatedAdt) + 1)

step("ADT: CLR normalisation")
clr <- function(m) {
  logged <- log1p(m)
  sweep(logged, 2, colMeans(logged))
}
adt <- altExp(es, "ADT", subsetName = "clean")
assay(adt, "clr") <- clr(assay(adt, "decontPro"))
altExp(es, "ADT", subsetName = "clean") <- adt
check("clean stores decontPro and clr; counts stays linked",
      identical(assayNames(ownAltExp(es, "clean", "ADT")),
                c("decontPro", "clr")))

# ---------------------------------------------------------------------
step("Both: T cells from RNA labels, checked with protein")
tcellLabels <- c("CD4+ T cells", "CD8+ T cells", "T cells")
isTcell <- predictions$labels %in% tcellLabels
es <- createSubset(es, subsetName = "tcells", cols = which(isTcell),
                   parentAssay = "clean")
tcellBarcodes <- subsetColnames(es, "tcells")
tcellAdt <- altExp(es, "ADT", subsetName = "tcells")
check("tcells sees counts, decontPro and clr through clean",
      identical(assayNames(tcellAdt), c("counts", "decontPro", "clr")))
check("tcells' ADT counts are the raw ADT counts of its cells",
      sameValues(counts(tcellAdt), rawAdt[, tcellBarcodes]))
check("tcells' decontPro is clean's, for its cells",
      sameValues(assay(tcellAdt, "decontPro"),
                 decontaminatedAdt[, tcellBarcodes]))
cleanClr <- assay(altExp(es, "ADT", subsetName = "clean"), "clr")
cd3InT <- median(cleanClr["CD3", isTcell])
cd3NotT <- median(cleanClr["CD3", !isTcell])
message("median CD3 CLR: T cells ", round(cd3InT, 2), ", others ",
        round(cd3NotT, 2))
check("RNA-labelled T cells have more CD3 protein", cd3InT > cd3NotT)
isCd4 <- predictions$labels == "CD4+ T cells"
isCd8 <- predictions$labels == "CD8+ T cells"
check("CD4+ T cells have more CD4 than CD8a protein; CD8+ the reverse",
      median(cleanClr["CD4", isCd4]) > median(cleanClr["CD8a", isCd4]) &&
        median(cleanClr["CD8a", isCd8]) > median(cleanClr["CD4", isCd8]))

step("RNA: T-cell re-analysis with getSubsetExperiment()")
tcells <- getSubsetExperiment(es, subsetName = "tcells")
check("the T cells' SingleCellExperiment has ADT with all three assays",
      identical(assayNames(altExp(tcells, "ADT")),
                c("counts", "decontPro", "clr")))
tcellLognormalized <- assay(tcells, "lognormalized")
tcellVariances <- modelGeneVariances(tcellLognormalized)
tcellVariable <- chooseHighlyVariableGenes(
  tcellVariances$statistics$residuals, top = 2000)
tcellScaled <- t(scale(t(as.matrix(tcellLognormalized[tcellVariable, ]))))
tcellPca <- runPca(tcellScaled)
reducedDim(tcells, "reduction.pca") <- t(tcellPca$components)
tcells$cluster <- clusterGraph(buildSnnGraph(tcellPca$components))$membership
sizeBefore <- mb(es)
es <- setSubsetExperiment(es, subsetName = "tcells", value = tcells)
check("storing the T-cell results copies no ADT",
      is.null(ownAltExp(es, "tcells", "ADT")))
check("storing the T-cell results copies no RNA assay",
      length(assayNames(.internalAssay(.subsets(es)[["tcells"]]))) == 0)
message("object grew by ", round(mb(es) - sizeBefore, 1), " MB")

step("ADT: decontPro again on the T cells alone")
tcellType <- as.integer(factor(tcells$cluster))
tcellPro <- decontPro(as.matrix(counts(tcellAdt)), tcellType,
                      ambient_counts = ambientAdt)
tcellDecontaminated <- tcellPro$decontaminated_counts
dimnames(tcellDecontaminated) <- dimnames(tcellAdt)
adt <- altExp(es, "ADT", subsetName = "tcells")
assay(adt, "decontPro") <- tcellDecontaminated
blocked <- tryCatch({
  altExp(es, "ADT", subsetName = "tcells") <- adt
  "stored"
}, error = function(e) conditionMessage(e))
check("replacing the linked decontPro without unlink is an error naming it",
      grepl("'ADT/decontPro' in subset 'tcells' is linked to altExp 'ADT' of subset 'clean'",
            blocked, fixed = TRUE))
messages <- character(0)
withCallingHandlers(
  altExp(es, "ADT", subsetName = "tcells", unlink = "ADT/decontPro") <- adt,
  message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
messages <- ownMessages(messages)
message("messages: ", paste(messages, collapse = " | "))
check("unlink = \"ADT/decontPro\" stores it, with a message",
      length(messages) == 1 &&
        grepl("now stores its own 'ADT/decontPro'", messages))
check("tcells stores only its decontPro; counts and clr stay linked",
      identical(assayNames(ownAltExp(es, "tcells", "ADT")), "decontPro"))
check("clean's decontPro is unchanged",
      sameValues(assay(altExp(es, "ADT", subsetName = "clean"), "decontPro"),
                 decontaminatedAdt))
sizes["end of workflow"] <- mb(es)

# ---------------------------------------------------------------------
step("Changes to the main ADT carry along")
renamed <- es
altExpNames(renamed) <- "protein"
check("renaming the main ADT renames what clean and tcells store of it",
      identical(altExpNames(renamed, subsetName = "tcells"), "protein") &&
        identical(assayNames(altExp(renamed, "protein", subsetName = "tcells")),
                  c("counts", "decontPro", "clr")))
isotype <- grepl("control", rownames(altExp(es, "ADT")))
noIsotype <- es
warnings <- character(0)
withCallingHandlers(
  altExp(noIsotype, "ADT") <- altExp(noIsotype, "ADT")[!isotype, ],
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
message(paste(warnings, collapse = "\n"))
check("dropping the isotype controls drops clean's and tcells' parts, with warnings",
      length(warnings) == 2 && is.null(ownAltExp(noIsotype, "clean", "ADT")) &&
        is.null(ownAltExp(noIsotype, "tcells", "ADT")))
check("clean then sees the 29 remaining antibodies' counts only",
      identical(dim(altExp(noIsotype, "ADT", subsetName = "clean")),
                c(29L, length(cleanBarcodes))) &&
        identical(assayNames(altExp(noIsotype, "ADT", subsetName = "clean")),
                  "counts"))

# ---------------------------------------------------------------------
step("Summary")
print(subsetSummary(es))
message("\nobject size (MB):")
print(sizes)
print(results, row.names = FALSE)
message("\n", sum(results$status == "PASS"), " of ", nrow(results),
        " checks passed in ", format(round(Sys.time() - t0, 1)))
