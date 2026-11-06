# Clear errors for a bad subsetName in every generic (#97)

- **Status:** Approved in chat
- **Date:** 2026-11-06
- **Branch:** `fix/set-subset-assay-null`, from `devel`
- **Lands:** one PR to `devel` (1.23.39); devel only, no release port

## Intent

`setSubsetAssay(es, "s1", m, subsetAssayName = NULL)` fails with "unable
to find an inherited method for function 'setSubsetAssay'" instead of the
package's own message (#97). The cause is S4 dispatch: the methods require
`subsetName = "character"` and `subsetAssayName = "character"`, so any
other value never reaches the checks in `.setSubsetAssay()`.

The same happens in every package generic whose methods require
`subsetName = "character"`. Measured on devel (1.23.38), `subsetName =
NULL`, `NA` and `1L` all give the dispatch error in each of them. After
this change, these reach the package's checks and get the messages from
#83: "'subsetName' must be one subset name, a single character string."

Success: every call that works today gives the same result; every call
listed below that now fails at dispatch gives an error naming the
argument and the problem instead.

Not a release port: no result is wrong; only the error is unclear.

## Changes

### Dispatch on `x` only

Add `signature = "x"` to these generics in `R/allGenerics.R`, as
`createSubset()` already has, and drop every argument but `x` from their
methods' signatures in the five `R/Subset<Class>.R` files:

`setSubsetAssay`, `getSubsetAssay`, `getSubsetExperiment`,
`setSubsetExperiment`, `subsetDim`, `subsetColData`, `subsetRowData`,
`subsetColData<-`, `subsetRowData<-`, `subsetColnames`, `subsetRownames`,
`subsetColnames<-`, `subsetRownames<-`, `subsetColLinks`,
`subsetRowLinks`, `subsetSpatialCoords`, `subsetSpatialData`.

The other package generics (`subsetAssayNames`, `subsetParent`,
`reducedDimNames<-`, `subsetSpatialData<-`, ...) already have methods
that require only `x`, so they can't fail this way and are unchanged.
Generics from other packages (`assay`,
`reducedDim`, `altExp`, `colData<-`, ...) are unchanged: their dispatch
on `i`, `type` or `e` is theirs, and a numeric `i` deliberately reaches
the parent class's method (#85).

Each internal implementation must check `subsetName` (and
`subsetAssayName` where it applies) before using it. Most already do via
`.isSubset()` or `.checkSubsetName()`; any that don't get the check.

### `value` in `subsetColData<-` and `subsetRowData<-`

These methods also required `value = "DataFrame"`, so a `data.frame`
failed at dispatch. They now accept what `colData(x, subsetName =) <-`
and `rowData(x, subsetName =) <-` accept: a `DataFrame`; a `data.frame`,
converted with `DataFrame(value, check.names = FALSE)`; or `NULL`, which
removes the subset's own columns. Anything else gives "'value' must be a
DataFrame, a data.frame or NULL." The conversion and check move into the
shared internal setters (`.subsetColData<-`, `.subsetRowData<-`), so the
`colData<-`/`rowData<-` path and the `subsetColData<-` path use one
implementation; `rowData(x, subsetName =) <- 1:3` then gets the same
message instead of SummarizedExperiment's.

## Tests (first, in `tests/testthat/test-ExperimentSubset.R`)

- One loop over the generics above, calling each with `subsetName =
  NULL`, `NA` and `1L` and expecting the `subsetName` message (the
  spatial ones in the SpatialExperiment block, the links ones in the
  TreeSummarizedExperiment block).
- `setSubsetAssay()`: add `NULL`, `NA` and `1L` to the existing
  `subsetAssayName` loop.
- `subsetColData<-` and `subsetRowData<-`: a `data.frame` is stored and
  reads back equal to the `DataFrame`; `NULL` removes the subset's own
  columns; `1:10` gives the `value` message.

Each must fail on devel on its assertion.

## Docs, NEWS, version

- `make docs`: only the `\S4method{}` usage lines change.
- `subsetColData<-`/`subsetRowData<-` `@param value`: say a
  `data.frame` or `NULL` is accepted.
- NEWS for 1.23.39, one entry citing #97.
- Version 1.23.38 -> 1.23.39.

## Checks

`make test-one FILTER=ExperimentSubset` while developing; `make test`,
`make coverage`, `make lint` before hand-off; `make check-full` and
`make bioccheck` before the PR.

## Trade-offs

- A later method can't dispatch on `subsetName`'s class (for example a
  numeric subset index) without changing the generics again. Nothing
  does today.
- `showMethods()` output and the Rd usage lines change. No other package
  defines methods on these generics.
