# 0005. Add an ExperimentSubset class union over the Subset* classes

- **Status:** Accepted
- **Date:** 2026-11-04
- **Deciders:** Joshua D. Campbell

## Context

The package has five classes, `SubsetSummarizedExperiment`,
`SubsetRangedSummarizedExperiment`, `SubsetSingleCellExperiment`,
`SubsetTreeSummarizedExperiment` and `SubsetSpatialExperiment`. Each
extends its Bioconductor class directly and adds one slot, `subsets`.
They have no common superclass.

Other packages therefore can't test for "any ExperimentSubset object"
without listing all five classes or checking for the `subsets` slot
(issue #5). decontX's planned ExperimentSubset support needs that test:
`SubsetSingleCellExperiment`, `SubsetTreeSummarizedExperiment` and
`SubsetSpatialExperiment` all dispatch to its SingleCellExperiment
methods, so the method has to branch on whether its input is an
ExperimentSubset object. Testing one class is not enough:
`is(x, "SubsetSingleCellExperiment")` is `FALSE` for a
`SubsetSpatialExperiment`.

The package has the same need internally. `.hasSubsets()` and
`.combineSubsets()` in `R/ExperimentSubset.R` test it with
`methods::.hasSlot(x, "subsets")`.

Objects of the five classes may be saved to disk by users, so a change
to their definitions has to keep those objects valid.

## Decision

We will define a class union, `ExperimentSubset`, whose members are the
five `Subset*` classes, in `R/allClasses.R` after their definitions,
and export it with `@exportClass ExperimentSubset`. The five class
definitions do not change.

`is(x, "ExperimentSubset")` becomes the supported test, for other
packages and inside this one: `.hasSubsets()` and `.combineSubsets()`
use it in place of the `subsets` slot check.

The class shares its name with the `ExperimentSubset()` constructor.
R keeps classes and functions apart, as for `SummarizedExperiment`.

## Consequences

- `is(x, "ExperimentSubset")` is `TRUE` for an object of any of the
  five classes, or of a subclass of one, and `FALSE` otherwise.
- Objects saved with earlier versions pass the test too, because class
  union membership is read from the class definitions loaded at run
  time, not stored in the object.
- No method in this package is defined on the union, so dispatch in
  this package is unchanged.
- Other packages can define methods on `ExperimentSubset`, but a
  package that defines methods for both `ExperimentSubset` and
  `SingleCellExperiment` gets different methods for different member
  classes, with no message. For a `SubsetSingleCellExperiment`, both
  are one step away, and R picks the `SingleCellExperiment` method
  because that parent comes first. For a `SubsetSpatialExperiment` or
  `SubsetTreeSummarizedExperiment`, `SingleCellExperiment` is two steps
  away, so R picks the `ExperimentSubset` method. Branching with `is()`
  inside the SingleCellExperiment method, as decontX plans to, gives
  the same result for every member class.
- A new `Subset*` class must be added to the union's members, in
  `R/allClasses.R`, as well as following the other steps for a new
  class. A test fails if a `Subset*` class in the namespace doesn't
  extend the union.

## Alternatives considered

- **An exported predicate only, `isExperimentSubset(x)`.** Answers the
  question, but other packages can't dispatch on a function, and it
  would need its own definition of "ExperimentSubset object" (the slot
  check).
- **The class union plus an `isExperimentSubset()` predicate.** The
  predicate would only wrap `is()`, adding a second exported name for
  the same test.
- **A virtual superclass holding the `subsets` slot**, with each class
  changed to, for example,
  `contains = c("SingleCellExperiment", "ExperimentSubset")`. It rewrites
  all five class definitions and moves the shared slot under multiple
  inheritance, which risks objects saved by earlier versions and
  changes dispatch, for nothing issue #5 needs. Other packages'
  methods would dispatch as they do with the union.
