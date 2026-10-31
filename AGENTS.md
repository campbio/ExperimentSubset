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
is no Shiny app, so that part of the standards doesn't apply.

The pkgdown site is at https://campbio.github.io/ExperimentSubset/: the
release version at the root, the devel version under `/dev/`. CI builds
and deploys both (`.github/workflows/pkgdown.yaml`, ADR 0002); never build
or deploy the site locally. Add new exports to `_pkgdown.yml` and run
`make site-check`.

## Layout

- `R/allClasses.R`: the S4 classes (see Object model).
- `R/allGenerics.R`: every package generic, with its roxygen docs and
  examples.
- `R/ExperimentSubset.R`: the `ExperimentSubset()` constructor, which
  wraps an input object in the matching `Subset*Experiment` class (a
  plain list of assays becomes a `SubsetSingleCellExperiment`) and can
  create a first subset through `subset =`, plus the shared internal
  implementations (`.createSubset()`, `.setSubsetAssay()`,
  `.subsetParent()`, ...) and the internal slot accessors.
- `R/Subset<Class>.R`, one per supported class: thin `setMethod()`
  wrappers that forward to the shared internals. A change to a method
  usually needs the same change in all five files.
- `Collate:` in DESCRIPTION fixes the load order (classes, then
  generics, then the rest). Add new R files there too.
- `vignettes/ExperimentSubset.Rmd`, with the figure `es_structure.png`:
  the short Bioconductor vignette, which is the site's "Get started" page.
- `vignettes/articles/`: longer pkgdown-only articles, not in the
  Bioconductor build (`.Rbuildignore`). Their packages go in
  `Config/Needs/website` in DESCRIPTION, not Suggests. Render one with
  `make article FILTER=<name>`, which needs the package installed.
- Style: the code uses **2-space indentation**, set in `.lintr`
  (`indentation_linter(indent = 2L)`). Converting to 4 spaces would be a
  separate maintainer decision and its own PR.

## Object model

- `AssaySubset` (internal): one subset. Slots `subsetName`, `rowIndices`
  and `colIndices` (indices into the parent), `parentAssay` (the name of
  the assay or subset it was made from, so subsets can chain), and
  `internalAssay` (a SummarizedExperiment-like object holding the
  subset's own assays and metadata).
- `SubsetSummarizedExperiment`, `SubsetRangedSummarizedExperiment`,
  `SubsetSingleCellExperiment`, `SubsetTreeSummarizedExperiment` and
  `SubsetSpatialExperiment`: each extends its Bioconductor class with one
  slot, `subsets`, a SimpleList of `AssaySubset`.
- `R/allClasses.R` exports each of these classes with
  `@exportClass <Name>`, not `@export`. Since roxygen2 7.2.1, `@export`
  on `.Foo <- setClass("Foo")` also exports the generator, which would
  add dot-prefixed exports to NAMESPACE. Follow the same pattern for any
  new class.
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
- It runs offline: a small simulated data set stands in for the PBMC
  data the vignette downloads.
- The whole file skips when scds (for `sce_chcl`) is missing; tests that
  need airway, scrapper or scater skip on their own.

## Extra make targets

None.

## Setup in a new worktree

- Install every package in Suggests (`devtools::install_dev_deps()`);
  checks don't relax missing Suggests.
- DESCRIPTION records the roxygen2 version in
  `Config/roxygen2/version: 8.1.0`. With roxygen2 8.1.0, `make docs`
  changes only the Rd files whose roxygen comments changed. If a newer
  roxygen2 rewrites many files in `man/` and updates that field, keep
  that refresh in its own PR rather than mixing it into a functional
  change.

## Related packages

None. No Bioconductor package depends on, imports or suggests
ExperimentSubset.

## Overrides

None.
