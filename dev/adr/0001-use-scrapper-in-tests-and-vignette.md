# 0001. Use scrapper in place of deprecated scran and scuttle functions

- **Status:** Accepted
- **Date:** 2026-10-21
- **Deciders:** Joshua D. Campbell

## Context

The tests and the vignette's PBMC example use scran and scuttle (through
scater) to normalize counts, choose highly variable genes and compute QC
metrics: `scater::normalizeCounts()`, `scran::modelGeneVar()`,
`scran::getTopHVGs()` and `scater::perCellQCMetrics()`. scuttle 1.22
deprecates every function that has a counterpart in scrapper, and scran
deprecates `getTopHVGs()` and the helpers `modelGeneVar()` calls, so the
test suite raises five deprecation warnings (issue #9). Once these are
defunct, the tests fail and readers copying the vignette get errors.

ExperimentSubset only calls these packages in tests and the vignette,
never in package code, so the dependency is a Suggests entry. The
vignette's PBMC chunks are `eval = FALSE`; the tests are the only place
the workflow runs.

## Decision

We will use scrapper's matrix-level functions (`centerSizeFactors()`,
`normalizeCounts()`, `modelGeneVariances()`,
`chooseHighlyVariableGenes()`, `computeRnaQcMetrics()`, `runPca()`) in
the tests and the vignette, add scrapper to Suggests, and remove scran
from Suggests, since nothing uses it any more. scater stays in Suggests:
the tests use `scater::calculatePCA()`, which is not deprecated.

The tests' highly variable gene fixture uses the same scrapper calls as
the vignette, so the suite exercises the workflow the vignette shows.

## Consequences

- The deprecation warnings go away, and the tests and vignette follow
  the maintained single-cell toolkit.
- scrapper's matrix-level functions take and return matrices, which fits
  ExperimentSubset's pattern of reading `assay(es, "<subset>")` and
  storing results back into a subset. Its `.se` wrappers write to the
  top-level object, not a subset, so they don't fit the subset steps.
- The vignette no longer loads scater, which avoids readers having to
  know that scater's and scrapper's `normalizeCounts()` mask each other.
- `runPca()` in the vignette is not run by any test.

## Alternatives considered

- **Choose subset rows without any variance modelling** (for example the
  top genes by row sum). This avoids a new Suggests entry but leaves the
  vignette showing deprecated functions.
- **scrapper's `.se` wrappers** (`quickRnaQc.se()`,
  `normalizeRnaCounts.se()`, `chooseRnaHvgs.se()`, `runPca.se()`). They
  read and write the top-level object rather than a subset, so the
  vignette would need to build a temporary SingleCellExperiment for each
  subset step and copy results back.
