# Require assay names everywhere, not only in the constructor (#124)

- **Status:** Approved (2026-11-25)
- **Date:** 2026-11-25
- **Branch:** `fix/assay-names`, from `devel`
- **Lands:** one PR to `devel` (1.23.60); devel only, no release port.
- **ADR:** none; the object model doesn't change.

## Intent

The package looks assays up by name, and `ExperimentSubset()` refuses
main assays without names, but other entry points don't check. On devel
(a21cb67) these accept assays without names:

- `assays(es) <- list(counts = m, m)` and `assay(es, 2) <- m`;
- `ExperimentSubset(sce)` where an altExp of `sce` has unnamed assays;
- `altExp<-` and `altExps<-` on the main object;
- `altExp<-` and `altExps<-` with `subsetName`, and
  `setSubsetExperiment()`, for a new altExp or one whose shape changes.

`assays(es) <- list(m)` fails with "attempt to set an attribute on
NULL": `.assays<-` takes it for a rename (#126) and `.renamed()` breaks
on `NULL` names.

Assay names can also be `NA`, which the constructor accepts.

## Rule

An assay without a name (`NULL`, `""` or `NA`) can't be added, and can't
stay through a change to the assays of the object that holds it. Objects
saved with such assays load and work until their assays change: the
#121 rule, everywhere.

## Changes

1. `.hasUnnamedAssays()` counts `NA` as unnamed; `.checkAssayNames()`
   uses it, so the constructor also refuses `NA` names, with its message.
2. The constructor also checks every altExp of a SingleCellExperiment
   input: "altExp 'HTO' has assays without names. Name them with
   assayNames<- first."
3. `.assays<-` (which `assay(x, i) <-` with a numeric `i` reaches) stops
   when the result has unnamed assays and its assays changed: "Assays
   without names are not supported. Name every assay, with assayNames<-
   or assay(x, \"name\") <- value." With subsets, a rename keeps #126's
   message, which comes first; `.followAssayRenames()` also checks,
   which fixes the crash.
4. `.checkAltExpAssayNames(current, value)`: for each altExp in `value`
   with unnamed assays, stops if `current` has none of that name (item
   2's message), or else with `.checkSameAssays()` (#121's message) if
   its assays changed. Called by the main `altExp<-` and `altExps<-`
   (current: `x`'s altExps), by both with `subsetName` and by
   `setSubsetExperiment()` (current: the altExps the subset sees), and
   by the constructor (current: none). Not in `.planAltExps()`, so
   `dependents = "keep"` copies of a saved object's altExps still work.

## Tests

Test first, each failing on its assertion: every path above, `NA`
names, the `list(m)` crash, and a saved object (built with
`as(sce, "SubsetSingleCellExperiment")`, without the constructor's
checks) that still works while unchanged. The #36 and #121 tests that
put unnamed assays in through `assays<-` or `ExperimentSubset()` build
their object that way instead; their assertions stay.

## Docs, NEWS, version

One sentence on the `ExperimentSubset()` help page and on the assay and
altExp setter pages; `make docs`. NEWS entry under 1.23.60, marked as a
change in behavior, with the crash fix.

## Out of scope

`setSubsetExperiment()` with an unnamed main assay stops with a
confusing message; logged as its own issue.
