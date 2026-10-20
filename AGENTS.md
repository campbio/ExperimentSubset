# ExperimentSubset: notes for coding agents

The shared development standards (r-bioc-dev-standards) load automatically
at session start in Claude Code. This file adds only what is specific to
this package. Other agents (for example through GEMINI.md) don't get them
automatically: read
https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/v1/standards.md
(or the cached copy in `~/.cache/r-bioc-dev-standards/v1/`) before
starting, and follow it; those agents also aren't bound by
`.claude/settings.json`.

In Claude Code, if the standards aren't in your context at session start,
this clone has no `.claude/settings.json`, so neither the standards nor
the permission guardrails are active: stop, and ask the developer to run
`make claude-setup` in a terminal outside Claude Code (not with `!`, which
runs inside it and is refused), then restart `claude`.

## About

ExperimentSubset adds subset support to Bioconductor Experiment objects:
many subsets of the same data (rows, columns, and their own assays,
rowData, colData, reducedDims, and so on) live inside one object instead
of separate copies. Supported classes are SummarizedExperiment,
RangedSummarizedExperiment, SingleCellExperiment, TreeSummarizedExperiment
and SpatialExperiment. Distributed through Bioconductor since 3.12. There
is no Shiny app and no pkgdown site, so those parts of the standards don't
apply, including `make site-check`.

## Layout

- `R/allClasses.R`: the S4 classes (see Object model).
- `R/allGenerics.R`: every package generic, with its roxygen docs and
  examples.
- `R/ExperimentSubset.R`: the `ExperimentSubset()` constructor, which
  wraps an input object in the matching `Subset*Experiment` class, plus
  the shared internal implementations (`.createSubset()`,
  `.setSubsetAssay()`, `.subsetParent()`, ...) and the internal slot
  accessors.
- `R/Subset<Class>.R`, one per supported class: thin `setMethod()`
  wrappers that forward to the shared internals. A change to a method
  usually needs the same change in all five files.
- `Collate:` in DESCRIPTION fixes the load order (classes, then
  generics, then the rest). Add new R files there too.
- `vignettes/ExperimentSubset.Rmd`, with the figure `es_structure.png`.
- Style: the code uses **2-space indentation**, set in `.lintr`
  (`indentation_linter(indent = 2L)`). Converting to 4 spaces would be a
  separate maintainer decision and its own PR.

## Object model

- `AssaySubset` (internal): one subset. Slots `subsetName`, `rowIndices`
  and `colIndices` (indices into the parent), `parentAssay` (the assay or
  subset it was made from, so subsets can chain), and `internalAssay` (a
  SummarizedExperiment-like object holding the subset's own assays and
  metadata).
- `SubsetSummarizedExperiment`, `SubsetRangedSummarizedExperiment`,
  `SubsetSingleCellExperiment`, `SubsetTreeSummarizedExperiment` and
  `SubsetSpatialExperiment`: each extends its Bioconductor class with one
  slot, `subsets`, a SimpleList of `AssaySubset`.
- Never use `@`. Inside the package, use the internal accessors in
  `R/ExperimentSubset.R` (`.subsets()`, `.subsets<-`, `.subsetName()`,
  `.rowIndices()`, `.colIndices()`, `.parentAssay()`, `.internalAssay()`,
  `.internalAssay<-`). Users see the exported API (`createSubset()`,
  `setSubsetAssay()`, `getSubsetAssay()`, `subsetNames()`,
  `subsetRowData()`, `subsetColData()`, ...) and the base methods
  (`assay()`, `reducedDims()`, `metadata()`, `altExps()`) called with
  `subsetName =`.

## Tests

- One file, `tests/testthat/test-ExperimentSubset.R` (testthat edition
  2, with `context()`), so `make test-one FILTER=ExperimentSubset` runs
  the whole suite.
- It's slow on first run: it downloads `TENxPBMCData("pbmc4k")` from
  ExperimentHub, then reuses the local cache.
- It loads scds, TENxPBMCData, scater, scran and airway with `library()`
  at the top, so it fails, rather than skips, when any is missing.

## Extra make targets

None.

## Setup in a new worktree

- Install every package in Suggests (`devtools::install_dev_deps()`);
  checks don't relax missing Suggests.
- DESCRIPTION has `RoxygenNote: 7.1.1`. `make docs` with a current
  roxygen2 rewrites every file in `man/` and updates that field. Keep
  such a refresh in its own PR rather than mixing it into a functional
  change.

## Related packages

None. No Bioconductor package depends on, imports or suggests
ExperimentSubset.

## Overrides

None.
