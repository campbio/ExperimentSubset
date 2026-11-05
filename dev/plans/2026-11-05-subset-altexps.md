# Subsets see the main object's altExps (#58)

- **Status:** Approved and implemented
- **Date:** 2026-11-05
- **Branch:** `feature/subset-altexps`, from `devel`
- **Lands:** one PR to `devel` (1.23.37); devel only, no release port
- **ADR:** none. This extends the linked-assay model (#2, #53) to
  altExps; the classes don't change.

## Intent

Issue #58: with `subsetName =`, `altExp()` and `altExpNames()` see only
the altExps stored on the subset, and `getSubsetExperiment()` leaves the
main object's altExps out. decontX needs a subset's ADT counts:
decontPro reads `assay(altExp(x, "ADT", subsetName = s), "counts")` and
writes its result as a new assay, `decontPro`, into that altExp, with
`altExp(x, "ADT", subsetName = s) <-` (decontX ADR 0005).

After this change a subset sees the main object's altExps, cut to the
subset's cells, as it sees linked assays. It stores only what it
changes, assay by assay, so decontPro's write-back stores the 30 x 300
`decontPro` matrix and leaves ADT `counts` linked.

Decided in chat (2026-11-05):

- Linking is per assay inside an altExp, not per whole altExp.
- `unlink` names an altExp's assay as `"altExpName/assayName"`
  (`"ADT/counts"`) everywhere. A bare name is always a main assay.
- An altExp is cut by columns only. A subset's `rows` pick genes of the
  main object, never features of an altExp, and a subset sees all of an
  altExp's features and assays. Subsets made from an altExp, with their
  own features, are #92 and out of scope here.
- reducedDims and metadata of a subset itself stay the subset's own:
  an embedding computed on all cells is not the subset's embedding.

## Example

```r
# main: 20,000 genes x 1,000 cells; altExp "ADT": 30 proteins x 1,000 cells
es <- createSubset(es, "tcells", rows = hvgs, cols = tcellIdx)  # 2,000 x 300
altExpNames(es, subsetName = "tcells")                  # "ADT"
dim(altExp(es, "ADT", subsetName = "tcells"))           # 30 x 300

adt <- altExp(es, "ADT", subsetName = "tcells")
assay(adt, "decontPro") <- decontaminated
altExp(es, "ADT", subsetName = "tcells") <- adt         # no unlink needed
assayNames(altExp(es, "ADT", subsetName = "tcells"))    # "counts", "decontPro"
assayNames(altExp(es, "ADT"))                           # "counts"

assay(adt, "counts") <- fixed
altExp(es, "ADT", subsetName = "tcells") <- adt
#> Error: 'ADT/counts' in subset 'tcells' is linked to the main altExp
#> 'ADT'. Set unlink = TRUE or unlink = "ADT/counts" ...
```

## Changes (`R/ExperimentSubset.R`)

### Reading

- **Which altExps a subset sees.** A new `.altExpOwners(x, s)`, like
  `.assayOwners()`: the main object's altExp names in their order, then
  the names stored only on subsets down the chain (`.subsetChain()`), each
  subset's in its stored order.
- **One altExp as the subset sees it.** A new
  `.subsetAltExpView(x, s, name, own = TRUE)` builds it from the top of
  the chain down:
  - Start from the main object's altExp of that name, cut to the subset's
    columns with `.indicesIn(x, s, NULL, .colIndices)`. Rows are never
    cut. If the main object has none, start from the first subset down the
    chain that stores one, cut with `.indicesIn()` into that subset.
  - At each later subset down the chain to `s` that stores its own part
    of the altExp (to `s`'s parent only, when `own = FALSE`, which gives
    the linked altExp for the write checks), that part, cut to `s`'s
    columns, is merged in:
    - **Same structure** (same class, number of rows, and rownames): its
      assays are added, each hiding a linked assay of the same name.
      Its rowData, colData, metadata and, for a `SingleCellExperiment`,
      reducedDims each replace the linked one whole when not empty
      (zero columns or zero length means "use the linked one").
      Everything else (nested altExps, `rowPairs`, `colPairs`,
      `mainExpName`) comes from the linked altExp.
    - **Different structure:** the stored part replaces the linked altExp
      entirely. This reads objects saved before this change, whose subsets
      hold full copies of the same shape the same way as before, apart
      from also seeing linked assays the copy lacks.
  - The view's colnames are the subset's colnames.
- **Getters** (`.altExp()`, `.altExps()`, `.altExpNames()` with
  `subsetName`): a new `.withSubsetAltExps(x, s)` returns a copy of the
  subset's internal object holding every visible altExp
  (`altExps(copy, withDimnames = FALSE) <- views`), and the
  `SingleCellExperiment` method runs on it. So `e` by name or position
  (positions follow the visible order), `withDimnames` and
  `withColData` work as now. `.altExpNames()` uses
  `names(.altExpOwners())` and builds no views.

### Writing

- **One storing path.** A new `.storeAltExps(x, s, value, unlink)`
  takes the full list of altExps the subset should have:
  - For each altExp in `value` with no linked altExp of that name: it is
    stored whole as the subset's own.
  - With a linked altExp (`.subsetAltExpView(x, s, name, own = FALSE)`):
    - Same structure: each assay that is not linked, or linked but
      changed (`.sameAssay()`), is stored; unchanged linked assays are
      not. rowData, colData, metadata and reducedDims are stored only when
      they differ from the linked altExp's; otherwise the stored part has
      them empty.
    - Different structure: stored whole; every linked assay of that
      altExp counts as changed.
    - Changed linked assays need `unlink`: `TRUE`, or names including
      `"ADT/counts"`. Otherwise it is an error naming every blocked
      `"altExp/assay"`, with its source ("the main altExp 'ADT'" or
      "altExp 'ADT' of subset 'tcells'"), and suggesting the full
      `unlink = c(...)`, as `.checkLinked()` does for assays. All checks
      run before anything is stored. An assay the subset already stores
      itself (after an unlink, or in an altExp of another shape) is no
      longer linked and needs no `unlink` again; this was added after the
      spec review found that later writes otherwise failed.
  - The subset's own altExps (or own parts) that `value` doesn't have are
    removed. Linked altExps that `value` doesn't have stay visible: a
    subset can't remove what it doesn't store, as for linked assays.
  - Each unlinked assay gives a message once everything is stored:
    "Subset 'tcells' now stores its own 'ADT/counts', unlinked from the
    main altExp 'ADT'."
- **`.altExp<-()` and `.altExps<-()`** apply the `SingleCellExperiment`
  setter to `.withSubsetAltExps(x, s)`, which keeps its checks and
  `e`/`withDimnames`/`withColData` handling, then pass
  `altExps(result, withDimnames = FALSE)` to `.storeAltExps()`.
  `altExp(x, "ADT", subsetName = s) <- NULL` removes the subset's own
  part, so the linked ADT shows again.
- **`unlink`** reaches these setters through `...`. It is checked with
  `.checkUnlink()` and removed from the arguments passed on
  (`.otherArgs()` drops `unlink` as well as `subsetName`). Without
  `subsetName`, a non-`FALSE` `unlink` is an error, as for `assay<-`:
  "'unlink' only applies when storing in a subset." The class files'
  replacement methods need no change, since `unlink` is in `...`.
- **`.altExpNames<-()`** with `subsetName`: `value` gives one name per
  visible altExp, as the getter returns them, with no `NA` or duplicate
  (a duplicate would let an own altExp take a linked name and hide it
  without `unlink`). Renaming a linked altExp is
  an error ("'ADT' is linked to the main altExp 'ADT'; rename it in the
  main object"); the subset's own altExps are renamed.

### Changes to a linked altExp (added after code review)

A subset's part that adds to a linked altExp only makes sense while that
altExp exists with the same shape. Decided in chat (2026-11-05), "carry
along":

- After any altExp change at a holder (`altExp<-`, `altExps<-`,
  `altExpNames<-` on the main object, the same with `subsetName`, and
  `setSubsetExperiment()`), `.carryAltExpParts()` visits the subsets below
  that holder, parents first.
- **Rename:** their parts under the old name are renamed. A rename onto a
  name a subset below already stores is an error, and nothing changes.
- **Removed or reshaped** (other rows, rownames or class): a part that
  added to it is dropped, with a warning naming the subset and altExp. A
  part that already replaced it whole is kept as the subset's own altExp.
- New values of the same shape change nothing; subsets see them.

A subset can't clear a linked altExp's rowData, colData, metadata or
reducedDims: an empty one shows the linked one, as a subset can't remove
a linked assay. This is documented.

### `getSubsetExperiment()` and `setSubsetExperiment()`

- `.getSubsetExperiment()`: `altExps(out)` becomes every visible altExp.
- `.setSubsetExperiment()`: `altExps(value)` goes through
  `.storeAltExps()`, so a round trip doesn't copy linked altExps. Its
  `unlink` checks for assays and altExp assays run together before
  anything is stored, with one error naming every blocked assay, main and
  altExp; the messages print after everything is stored.

### Not changed

- The main object's own altExp methods (no `subsetName`).
- `[`, `cbind` and `rbind`: a subset's stored altExp parts are cut and
  combined with its internal object, as now; the main object's altExps
  are handled by `SingleCellExperiment`.
- SE and RSE subsets have no altExps.

## Behavior changes (NEWS)

- With `subsetName =`, `altExp()`, `altExps()` and `altExpNames()` include
  the main object's altExps and those of parent subsets, and
  `subsetSummary()` lists them.
- `altExpNames(x, subsetName =) <-` takes one name per visible altExp, and
  can't rename linked ones.
- `altExps(x, subsetName =) <-` can't remove a linked altExp.
- `getSubsetExperiment()` includes the visible altExps;
  `setSubsetExperiment()` no longer stores copies of unchanged ones.

## Docs

- `?"ExperimentSubset-altExp"`: linking, cutting by columns only, the
  `unlink` rule with `"altExp/assay"` names (`unlink` in the `...` list),
  and the removal and renaming rules. The example adds an altExp to the
  main object and reads it from the subset, then adds an assay to it
  through the subset.
- `?getSubsetExperiment` `@return` and `?setSubsetExperiment` `unlink`
  and `@return`.
- Vignette: the altExp section says a subset sees the main object's
  altExps.
- NEWS, version 1.23.37.

## Tests

Test first; each new test must fail on its assertion on devel.

- The issue's reproducer: `altExpNames()` and `altExp()` for `s1` give
  ADT for the 8 cells, with the subset's colnames; also by position and
  in `getSubsetExperiment()`.
- Live view: changing the main ADT after creating the subset shows in the
  subset.
- Adding an assay through the subset stores only that assay (checked by
  the subset's internal part), leaves the main ADT unchanged, and needs
  no `unlink`.
- A changed linked assay: an error without `unlink`, naming
  `'ADT/counts'`, with nothing stored; stored with `unlink = TRUE` and
  with `unlink = "ADT/counts"`, with the message; a bare `"counts"` does
  not unlink it.
- Chained subset: sees main assays and the parent subset's own altExp
  assays, cut to its cells.
- A stored altExp of a different shape hides the linked one.
- Own rowData or metadata on the altExp replace the linked ones; unchanged
  ones aren't stored.
- `altExp(..., subsetName =) <- NULL` shows the linked one again;
  `altExps<-` without a linked altExp leaves it visible.
- `altExpNames<-`: renaming an own altExp works; renaming a linked one is
  an error.
- `setSubsetExperiment()`: a round trip stores nothing for ADT; a changed
  ADT `counts` and a changed main `counts` give one error naming both;
  `withColData` and `withDimnames` still work.
- `unlink` on the main object's `altExp<-` is an error.
- TreeSE and SpatialExperiment subsets see the main object's altExps too.
- Existing tests: in "Testing altExperiments", the main object's altExps
  (`a1`, `a2`) are now also visible in `hvg1000`. The
  `altExpNames(es, subsetName = "hvg1000") <- c("a1", "a2")` call and
  the `altExpNames()` expectation after `altExps(..., subsetName =) <-`
  change to match the new rules; nothing else in them is weakened.

## Checks

`make test`, `make coverage` (no drop), `make check-full`,
`make bioccheck`, `make lint` on changed lines, `make docs`.
