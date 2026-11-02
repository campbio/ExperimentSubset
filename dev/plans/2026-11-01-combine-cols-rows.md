# Keep subsets in combineCols() and combineRows() (#73)

- **Status:** Draft for review
- **Date:** 2026-11-01
- **Branch:** `fix/combine-cols-rows`, from `devel`
- **Lands:** one PR to `devel` (1.23.29); devel only, no release port

## Intent

`cbind()` and `rbind()` keep the subsets of every object (#67), but
SummarizedExperiment's other combining functions, `combineCols()` and
`combineRows()`, return a plain object and drop every subset without a
message. These functions are for objects whose rows (for `combineCols`)
or columns (for `combineRows`) differ: they take the union by name and
fill the gaps. After this change they keep the subsets of every object
too, with their own assays, rowData, colData, reducedDims and altExps,
matching each subset's rows (or columns) by name into the union.

Success: for all five `Subset*` classes, `combineCols()` and
`combineRows()` return an `ExperimentSubset` object holding every
object's subsets, each with the same rows, columns and values as before
combining. Nothing that works today changes, except that subsets are kept
instead of dropped.

Not a release port: no results are wrong, and the loss shows in the
plain class of the result.

## What upstream does

Measured with the current SummarizedExperiment and SingleCellExperiment
methods (TreeSummarizedExperiment and SpatialExperiment have no methods
of their own):

| Class            | `combineCols()` returns | `combineRows()` returns |
|------------------|-------------------------|-------------------------|
| SE               | SE                      | SE                      |
| RSE              | RSE                     | RSE                     |
| SCE              | SCE                     | RSE                     |
| TreeSE           | SCE                     | RSE                     |
| SpatialExperiment| SCE                     | RSE                     |

- The union keeps the first object's rows (for `combineCols`) in their
  order, then appends each new name in order of first appearance.
- With `use.names = TRUE` (the default) upstream refuses `NULL` or
  duplicated row names (`combineCols`) or column names (`combineRows`),
  so matching by name is always unambiguous.
- With `use.names = FALSE` upstream requires the same number of rows
  (`combineCols`) or columns (`combineRows`) and combines them by
  position, keeping the first object's names.
- The combined dimension is concatenated by position, as in `cbind()`
  and `rbind()`; repeated names there are allowed.
- `combineCols()` on SpatialExperiment objects does not rename repeated
  `sample_id`s (it returns a SingleCellExperiment), unlike `cbind()`.
- Methods dispatch on the first object only.

## Decision: the class of the result

The result is the `Subset*` version of whatever the parent method
returns ("same as upstream, plus the subsets"). For example
`combineCols()` on a `SubsetTreeSummarizedExperiment` returns a
`SubsetSingleCellExperiment`, and `combineRows()` on a
`SubsetSingleCellExperiment` a `SubsetRangedSummarizedExperiment`.

Each subset's `internalAssay` has the class the main object extends
(`.createSubset()` makes it with `as(internal, <base class>)`), so when
the result's base class differs from a subset's `internalAssay` class,
the `internalAssay` is converted with `as()` to the result's base class.
That drops what the new class cannot hold, as upstream does for the main
object: for example a subset's reducedDims and altExps after
`combineRows()` on SCE, or the tree of a TreeSE subset.

Rejected: an error when upstream changes the class (stops calls that run
today), and rebuilding the input class from upstream's result (patches an
upstream gap from here; fragile).

## Design

### Methods

In each `R/Subset<Class>.R`, next to `cbind`/`rbind`, a `combineCols`
and a `combineRows` method for that class, with `@rdname
ExperimentSubset-combine` and, like `cbind`/`rbind`, no `@export`: the
generics are imported, so the methods dispatch once the package is
loaded. Each:

1. collects `objects <- list(x, ...)`;
2. checks subset names are unique across the objects (shared helper);
3. calls `callNextMethod()`, passing `use.names`, `delayed` and `fill`
   through unchanged;
4. hands over the subsets with `.combineSubsets(out, objects, along,
   use.names)`.

Method formals: `function(x, ..., delayed = TRUE, fill = NA, use.names =
TRUE)`. `use.names` is needed to build the maps; `delayed` and `fill` are
formals so that `list(x, ...)` doesn't collect them as objects.

`combineCols` and `combineRows` are generics in S4Vectors (already in
Imports); import them with `@importMethodsFrom S4Vectors combineCols
combineRows`, as `cbind`/`rbind` are imported. No DESCRIPTION change.

### Shared internals (`R/ExperimentSubset.R`)

- **`.combineSubsets(out, objects, along, useNames = FALSE)`**: `cbind`
  and `rbind` call it as now; `combineCols`/`combineRows` pass
  `use.names`. It wraps `out` with `ExperimentSubset()` when it isn't
  already an `ExperimentSubset` object (upstream builds a new object), so
  the methods stay short. Subsets made from a main-object assay send
  their `rowIndices` and `colIndices` through the maps from
  `.combineMaps()`; subsets of subsets keep theirs, which are relative
  to their parent subset. Parent names are made explicit first with
  `.qualifyParents()`, as now. If the class of a subset's
  `internalAssay` is not the base class of `out`, it is converted with
  `as()`. The SpatialExperiment `sample_id` step stays as it is; it runs
  only when `out` is a SpatialExperiment, so it does not affect the new
  functions.
- **`.combineMaps(objects, out, along, useNames)`**: one helper for all
  four functions. Along `along`, a shift by the size of the earlier
  objects; on the other dimension, `match(<names in that object>, <names
  in out>)` when `useNames` is `TRUE`, unchanged positions otherwise
  (which is what `cbind`/`rbind` need).
- **`.followMatchedNames(out, objects, sources, along)`**: with
  `use.names = FALSE`, the result takes the first object's names on the
  matched dimension, so a later object's subsets would keep names the
  main object no longer has. For each subset whose object's names differ
  from `out`'s on that dimension, this sets the stored names to `out`'s
  at the subset's positions, but only where they still equal the
  object's names there: names set with `subsetRownames<-` or
  `subsetColnames<-` are kept, as `x[i, j]` keeps them. Missing (`NULL`)
  names are left as they are. For `cbind`/`rbind` and for
  `use.names = TRUE` the object's names already match, so nothing
  changes.
- **`.checkSubsetNames(objects)`**: the unique-names check, moved out of
  `.checkCombine()` unchanged (same message). `.checkCombine()` calls it
  and keeps `cbind()`'s same-rows check, so `cbind`/`rbind` behave
  exactly as now.

### Docs

`?"ExperimentSubset-combine"` (in `R/ExperimentSubset.R`) gains:

- a paragraph on `combineCols()` and `combineRows()`: they keep every
  object's subsets, matching each subset's rows (or columns) by name into
  the union, or by position with `use.names = FALSE`;
- the class rule and what converting a subset can drop;
- that methods dispatch on the first object, so a plain first object
  (such as a `SingleCellExperiment`) gives a plain result without
  subsets; convert it with `ExperimentSubset()` first;
- `@param x` (the first object, for `combineCols`/`combineRows`),
  `@param use.names` (passed on; default `TRUE`) and `@param
  delayed,fill` (passed on; defaults `TRUE` and `NA`);
- an example: two objects with overlapping but different genes, a subset
  in each, `combineCols()`, then `subsetNames()` and `subsetRownames()`.

Then `make docs` and `make site-check` (the page is already in
`_pkgdown.yml`). The vignette and articles don't mention combining.

### Tests (`tests/testthat/test-ExperimentSubset.R`)

Written first; each must fail on its assertion before the fix.

- A variant of `combineExample()` whose objects share only some rows
  (for `combineCols`) or some columns (for `combineRows`), with names
  in a different order, so mapping by name differs from mapping by
  position.
- `combineCols()` and `combineRows()` keep every subset: names, the
  subset's row and column names, its own assay values (`assay(x, "s2")`
  and `"s2/doubled"`), its own colData column, the subset of a subset
  and the subset of a subset assay.
- `use.names = FALSE` leaves positions unchanged on the matched
  dimension.
- Repeated subset names across the objects are an error, with the
  existing message.
- An object without subsets in the middle shifts later subsets
  correctly.
- All five classes: the class of the result (SubsetSCE for TreeSE and
  SpatialExperiment `combineCols`, SubsetRSE for SCE, TreeSE and
  SpatialExperiment `combineRows`, otherwise the input class), that
  subsets are kept, and that each subset's `internalAssay` has the
  result's base class.
- Names set with `subsetRownames<-` / `subsetColnames<-` survive
  `cbind`, `combineCols` and `combineRows`, with `use.names = TRUE` and
  `FALSE`; with `FALSE`, unchanged names follow the first object's.
- A plain (unconverted) object in the middle of `combineCols()`.
- Existing `cbind`/`rbind` tests stay unchanged and must pass after the
  refactor.

### NEWS and version

Version 1.23.29; a NEWS entry describing the new behavior and the class
rule, ending `(#73)`.

## Checks

While developing: `make test-one FILTER=ExperimentSubset`. Before
hand-off: `make test`, `make coverage` (must not drop), `make lint` on
the touched files. Before the PR: `make check-full` and `make bioccheck`.
