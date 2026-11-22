# Makes es-1.23.54.rds: an object saved by ExperimentSubset 1.23.54,
# before AssaySubset gained the altExpName slot (ADR 0007). The tests
# check that such objects still work, with and without updateObject().
# It was run once, from the package root on devel at 1.23.54, with
#   devtools::load_all()
#   source("tests/testthat/fixtures/make-es-1.23.54.R")
# Don't run it again on a later version: the point is the old class.
stopifnot(packageVersion("ExperimentSubset") == "1.23.54")
set.seed(1)
counts <- matrix(rpois(60, lambda = 5), nrow = 6,
                 dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
es <- ExperimentSubset(list(counts = counts))
es <- createSubset(es, subsetName = "first", rows = 1:4, cols = 1:8)
es <- setSubsetAssay(es, "first", assay(es, "first") + 1, "plusOne")
es <- createSubset(es, subsetName = "second", rows = 2:3, cols = 3:6,
                   parentAssay = "first/plusOne")
saveRDS(es, "tests/testthat/fixtures/es-1.23.54.rds", version = 3)
