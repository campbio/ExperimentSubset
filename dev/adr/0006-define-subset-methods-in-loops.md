# 0006. Define the Subset* methods once, in loops over the classes

- **Status:** Accepted
- **Date:** 2026-11-19
- **Deciders:** Joshua D. Campbell

## Context

Each of the five `Subset*` classes has its own file,
`R/Subset<Class>.R`, holding thin `setMethod()` wrappers that forward to
shared helpers in `R/ExperimentSubset.R`. The five files total 2,706
lines and 255 method definitions, but only about 70 distinct methods:
37 are repeated in all five files, and 19 more (altExp and reducedDim)
in the three files for classes that extend SingleCellExperiment. A
change to a method has to be repeated in every copy, and the copies
have drifted: the `altExp<-` bug in #21 was in three files (#36).

The methods must stay on the concrete classes. Most are for generics
from SummarizedExperiment or SingleCellExperiment (`assay`, `colData`,
`altExp`, ...), and the helpers fall back to the parent method with
`callNextMethod()` when there's no `subsetName`. ADR 0005 keeps methods
off the `ExperimentSubset` union because, for those generics, a union
method and the parent class's method can be the same distance away.

roxygen documents a method from the `setMethod()` call after its block.
Each method needs an `\alias{}` in an Rd file, or R CMD check reports it
as undocumented.

## Decision

We will define each shared method once, inside a loop over an internal
vector of class names: `.esClasses` (the five classes) for the base
methods and `.esSCEClasses` (the three SingleCellExperiment-based
classes) for altExp and reducedDim. Both vectors live in
`R/allClasses.R` after the union. The loops live in a new file,
`R/SubsetMethods.R`. The methods of a single class stay written out in
`R/SubsetSpatialExperiment.R` and `R/SubsetTreeSummarizedExperiment.R`.
`R/SubsetSummarizedExperiment.R`, `R/SubsetRangedSummarizedExperiment.R`
and `R/SubsetSingleCellExperiment.R` are removed.

The aliases and `\S4method` usage entries for looped methods are
generated at documentation time by an internal helper, `.rdMethods()`,
called from roxygen `@eval` tags. It reads the methods registered after
`load_all()`.

## Consequences

- Each shared method is written once, so a fix reaches every class.
- Signatures, formals and bodies don't change, and every method stays
  on its concrete class, so dispatch and `callNextMethod()` are
  unchanged. A test snapshots the full list of methods to confirm it.
- Help pages keep their aliases and usage entries; only their order
  changes.
- A new shared method goes inside the right loop, and its Rd page needs
  an `@eval .rdMethods(...)` block. A plain `@rdname` block above a loop
  documents nothing. AGENTS.md says so.
- A new `Subset*` class must be added to `.esClasses` (and to
  `.esSCEClasses` if it extends SingleCellExperiment) as well as to the
  union. A test checks that `.esClasses` matches the union's members.
- The docs depend on a package helper run by roxygen. If it breaks,
  `make docs` fails or R CMD check reports undocumented methods; it
  can't fail silently.
- Follow-up in #36: the duplicated `arglist` block and near-copy helpers
  in `R/ExperimentSubset.R`, and a single parent-chain resolver.

## Alternatives considered

- **Methods on the `ExperimentSubset` union for the package's own
  generics, loops for imported generics.** The package's own generics
  have no competing methods, so a union method would be safe. But it
  reverses part of ADR 0005, changes the help aliases (for example to
  `subsetDim,ExperimentSubset-method`), and leaves two patterns side by
  side for little less code.
- **Generate the five files from a template with a script.** roxygen
  and dispatch wouldn't change, but the 2,700 duplicated lines would
  stay in the repo, and every edit would have to go through the
  generator.
- **Loops with aliases only, no generated usage.** A simpler helper, but
  the help pages would lose their per-class `\S4method` usage lines, a
  visible change in what is otherwise a pure refactor.
