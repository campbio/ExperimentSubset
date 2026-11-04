# A common ExperimentSubset class for the Subset* classes (#5)

- **Status:** Draft for review
- **Date:** 2026-11-04
- **Branch:** `feature/experimentsubset-class`, from `devel`
- **Lands:** one PR to `devel` (1.23.34); devel only, no release port
- **ADR:** [0005](../adr/0005-experimentsubset-class-union.md)

## Intent

Issue #5: other packages, decontX first, need to tell whether an object
is any ExperimentSubset object, whatever class it extends. Today the
five `Subset*` classes share no superclass, so that takes a list of
five classes or a check for the `subsets` slot:

```r
library(ExperimentSubset)
counts <- matrix(rpois(200, 5), 10, 20,
  dimnames = list(paste0("g", 1:10), paste0("c", 1:20)))
spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = counts))
ess <- ExperimentSubset(spe)
is(ess, "SubsetSingleCellExperiment")
#> [1] FALSE
```

After this change, `is(ess, "ExperimentSubset")` is `TRUE`.

This is a new feature, not a fix for wrong results, so it is devel
only.

## Changes

### The class (`R/allClasses.R`)

- After the five class definitions:
  `setClassUnion("ExperimentSubset", c("SubsetSummarizedExperiment",
  "SubsetRangedSummarizedExperiment", "SubsetSingleCellExperiment",
  "SubsetTreeSummarizedExperiment", "SubsetSpatialExperiment"))`.
- Exported with `@exportClass ExperimentSubset`, as the other classes
  are, giving `exportClasses(ExperimentSubset)` in NAMESPACE next to the
  existing `export(ExperimentSubset)` for the constructor.
- The five class definitions are unchanged.

### Internal use (`R/ExperimentSubset.R`)

- `.hasSubsets()` tests `methods::is(x, "ExperimentSubset")` for each
  object, in place of `methods::.hasSlot(x, "subsets")`.
- `.combineSubsets()` tests `methods::is(out, "ExperimentSubset")` in
  place of `methods::.hasSlot(out, "subsets")`.
- No change in behavior: only these classes and their subclasses have
  the `subsets` slot.

### Docs

- A roxygen block on the union gives a help page,
  `?"ExperimentSubset-class"`: what the class is for
  (`is(x, "ExperimentSubset")`), its five member classes, and that it is
  virtual, so there are no objects of the class itself. `@return` isn't
  needed for a class. The example builds a small SpatialExperiment, like
  the issue's, converts it with `ExperimentSubset()`, and shows that it
  is an `ExperimentSubset` but not a `SubsetSingleCellExperiment`, and
  that the plain SpatialExperiment isn't an `ExperimentSubset`. It
  follows the standards for user-facing code: one step per line, a
  comment per step, `set.seed()` before `rpois()`, and no `pkg::`
  prefixes (SpatialExperiment is in Depends).
- `_pkgdown.yml` already lists `ends_with("-class")` under "Constructor
  and classes", so the page is picked up without an edit;
  `make site-check` confirms it.
- The constructor's help page (`?ExperimentSubset`) gets a one-line
  `@seealso` pointing to `ExperimentSubset-class`, so readers of either
  page find the other.
- `make docs` regenerates `man/` and NAMESPACE.
- No vignette or article change: they don't test for the class, and
  novices don't need to.

### NEWS and version

- DESCRIPTION: 1.23.33 to 1.23.34.
- NEWS.md, under `# ExperimentSubset 1.23.34 (<date>)`: a new virtual
  class, `ExperimentSubset`, that all the `Subset*` classes extend, so
  `is(x, "ExperimentSubset")` tests for any ExperimentSubset object,
  whatever class it extends (#5).

### ADR

- `dev/adr/0005-experimentsubset-class-union.md`, status Proposed, and
  its row in the index in `dev/adr/README.md`. The maintainer changes it
  to Accepted on approval.

## Testing

In `tests/testthat/test-ExperimentSubset.R`, one new `test_that()`,
written first:

- `is(x, "ExperimentSubset")` is `TRUE` for the result of
  `ExperimentSubset()` on a SummarizedExperiment, a
  RangedSummarizedExperiment, a SingleCellExperiment, a
  TreeSummarizedExperiment, a SpatialExperiment, and a plain list of
  assays, each from small simulated counts.
- It is `FALSE` for a plain SingleCellExperiment and a plain
  SummarizedExperiment.
- SummarizedExperiment, SingleCellExperiment, TreeSummarizedExperiment
  and SpatialExperiment are all in Depends, so the test calls their
  constructors directly and needs no skip, like the existing tests
  for those classes.

Before the class exists, `is()` with an unknown class returns `FALSE`
without an error, so the `TRUE` expectations fail on their assertions.

The existing combine and `subsetName` tests cover the switch in
`.hasSubsets()` and `.combineSubsets()`.

Then `make test`, `make coverage` (must not drop), `make lint` on the
changed files, `make site-check`, and before a PR `make check-full` and
`make bioccheck`.

## Out of scope

- An `isExperimentSubset()` predicate (ADR 0005, alternatives).
- Methods defined on the union.
- Any change to decontX; that is its own work in celda.
