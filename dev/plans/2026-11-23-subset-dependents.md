# What happens to subsets when what they were made from goes (#122)

- **Status:** Draft, for review
- **Date:** 2026-11-23
- **Branch:** `feature/subset-dependents`, from `devel`
- **Lands:** one PR to `devel` (1.23.57); devel only, no release port.
- **ADR:** none; the object model doesn't change.
- **Follow-ups:** #126 (renaming a main assay leaves subsets pointing at
  the old name), #127 (`dependents = "keep"` when a whole altExp is
  removed).

## Intent

A subset's parent assay can go away while the subset stays. Today this
is handled three ways (#122):

- removing a main assay (`assay(x, "counts") <- NULL`,
  `assays(x) <- ...`) leaves the subset pointing at nothing:
  `assay(x, "s")` errors, while its own assays still work;
- removing a subset's own assay that another subset was made from
  (`setSubsetExperiment()`) stops (`.checkNotParents()`), and so does
  storing a subset's altExp without the assay an altExp-rooted subset was
  made from (`.checkAltExpParents()`);
- giving the main altExp other rows, which drops what a container subset
  stored in it (the #58 rule), leaves altExp-rooted subsets without their
  parent assay; removing a whole altExp drops its subsets with a warning.

No exported function removes a subset, so a user who wants to drop
`counts` to save memory has no way to clean up first.

Goal: every path behaves the same way, chosen with one argument; nothing
is lost unless the user asks; and subsets can be removed.

Decided in chat (2026-11-22 to 2026-11-23):

- **Default:** stop with an error, on every path, before anything
  changes.
- **Argument:** `dependents = c("error", "keep", "remove")`.
  - `"keep"`: the subsets stay intact, along with every subset below
    them, by storing their own slice of the data that is going away.
  - `"remove"`: those subsets, and every subset below them, are removed,
    with a warning.
- **Same logic everywhere:** main assays, subsets' own assays,
  altExp-rooted subsets, removing a whole altExp, and the new
  `removeSubset()`.
- **Exception:** `"keep"` when a whole altExp is removed stops with "not
  supported yet" (#127), because those subsets' rows would index nothing
  without a change to the object model.
- **Renames:** out of scope (#126).

## User-facing API

### `dependents`

`dependents` is a new argument with the values `"error"` (the default),
`"keep"` and `"remove"`. It is taken by:

| Path | How it reaches the argument |
|---|---|
| `assay(x, i) <- NULL`, `assays(x) <- value` on the main object | new `assays<-` method; SummarizedExperiment's `assay<-` calls `assays(x, withDimnames, ...)[[i]] <- value`, so `dependents` passes through `...` |
| `setSubsetExperiment(x, subsetName, value, unlink, dependents)` | new formal argument |
| `altExp(x, e) <-`, `altExps(x) <-`, with or without `subsetName =` | through `...`, as `unlink` is today |
| `altExpNames(x) <-` | not affected: a rename is followed (#92) |
| `removeSubset(x, subsetName, dependents)` | new function |

Any other value of `dependents` gives the error "'dependents' must be
"error", "keep" or "remove"." It is checked on every path, even when no
subset is affected, so a typo is always caught.

`combineRows()`, `combineCols()`, `cbind()`, `rbind()` and `x[i, j]`
don't take it and are unchanged: combining isn't a removal the user
asks for, and `combineRows()` keeps removing the affected subsets with
the current warning.

### Which subsets are affected

A subset is affected if what it was made from is visible before the
change and not after:

- its parent assay (`.primaryAssayName()`): the main object's assay,
  another subset's own assay, or an assay of the altExp it was made from;
- for a subset made with `createSubset(altExp =)`, the altExp itself;
- for `removeSubset()`, the subset it was made from.

Only the top affected subsets are acted on: those that lost what they
were made from directly, not only because a subset above them did.
Subsets below them follow: they are removed with `"remove"`, and with
`"keep"` they read through the top subset as before.

A subset that only sees the removed assay as an extra linked assay, such
as one made from `logcounts` that also sees `counts`, is not affected. It
stops seeing `counts`, as today.

### `"error"`

Stops before anything changes, naming what goes and the subsets made
from it:

```
Error: Assay 'counts' is the parent assay of subset 's'. Use
dependents = "keep" to keep it, with its own copy of that data, or
dependents = "remove" to remove it and the subsets made from it.
```

The wording varies with what goes away:

- "Assay 'counts' is the parent assay of subsets 's' and 'u'."
- "Assay 'scaled' of subset 's' is the parent assay of subset 't'."
- "Assay 'ADT/clr' of subset 'clean' is the parent assay of subset
  'clean.antibodies'."
- "Subset 's' is the parent of subset 't'." (from `removeSubset()`)
- "altExp 'ADT' is the parent of subset 'adt.sub'." For a removed altExp,
  the error offers only `"remove"`.

These messages replace the current `.checkNotParents()` and
`.checkAltExpParents()` error, "'value' has no assay ..., which subsets
... were made from. Keep it in 'value'."

### `"keep"`

Each top affected subset stores, as its own assay under the same name,
its slice of the data it is losing, taken from the object as it was
before the change. A message is shown for each, like the existing "now
stores its own rowTree" messages:

```
Subset 's' now stores its own assay 'counts'.
```

- The subset keeps its `parentAssay`, so `subsetSummary()` and
  `subsetParent()` still show `counts`, and its row and column indices
  don't change. Its own assays take priority over linked ones of the same
  name (`.assayOwners()`), so `assay(x, "s")` returns the kept copy.
- The copy is not linked: if an assay named `counts` is later added to
  the main object, `s` keeps its own.
- **When the main altExp gets other rows**, the kept slice is taken from
  the old altExp and then follows the rows, as stored assays already do:
  the features whose names remain are kept.
- **When a whole altExp is removed**, `"keep"` stops with:
  "dependents = "keep" isn't supported yet when an altExp is removed
  (#127). Use dependents = "remove"."
- **Removing a subset** has its own rules; see `removeSubset()` below.

### `"remove"`

Removes the top affected subsets and every subset below them
(`.subsetsBelow()`), with one warning:

```
Warning: Removed subsets 's' and 't': the assay 'counts' they were made
from was removed.
```

This generalizes the current `.dropSubsets()` warning ("Dropped subsets
... : the altExp 'ADT' they were made from was removed."), which becomes
"Removed ..." so that every path uses the same wording.

### `removeSubset(x, subsetName, dependents = "error")`

- `subsetName` is a character vector of one or more subset names; an
  unknown name is an error, as in other functions.
- The subsets named are removed. With subsets made from them that aren't
  also named:
  - `"error"` stops, naming them;
  - `"remove"` removes them too, with the warning;
  - `"keep"` re-roots each child of a removed subset onto that subset's
    parent (below).
- Returns `x` without the subsets, keeping its class.

**Re-rooting with `"keep"`.** Removing `s` from a chain p → s → c, where
p is the main object or a subset:

1. `c`'s row and column indices, which index `s`, are composed through
   `s`'s indices, so they index p.
2. `c` stores, as its own, its slice of each assay and each altExp that
   it saw from `s`'s own storage, and that it doesn't already store. An
   assay `s` linked from p stays linked, since `c` now reads p directly.
3. `c`'s `parentAssay` names the same assay one level up. `"s"` or
   `"s/scaled"` becomes `"scaled"` if p is the main object, or
   `"p/scaled"` if p is a subset; `c` now stores `scaled` itself. If `c`'s
   primary assay was one `s` linked from p, such as `counts`, it becomes
   `s`'s own `parentAssay`.
4. A message for each: "Subset 'c' now stores its own assays 'scaled' and
   'X'."
5. An altExp-rooted child of `s` (made with `createSubset(altExp =)` from
   `s`'s view of an altExp) is re-rooted on p's view of that altExp, with
   `s`'s stored parts of it copied into the child as in step 2. If p's
   view has other rows than `s`'s, because `s` stored a whole replacement
   altExp, it stops: "Subset 'c' can't be kept: altExp 'ADT' of subset
   's' has other rows than that of its parent. Use dependents =
   "remove"."
6. Chains: if several subsets in one chain are removed, the remaining
   children re-root onto the nearest parent that is kept, by applying
   steps 1 to 5 one removed subset at a time, from the top down.
7. `c`'s rowData and colData are its own and don't change. A subset
   doesn't see its parent's rowData or colData (#50), so there is
   nothing to copy.

## Internals

### One shared check

The new `.resolveDependents(before, after, dependents, cause)`, in
`R/ExperimentSubset.R`, runs on every path after the new object is built
and before it is returned:

1. Find the affected subsets: those in `before` and in `after` whose
   parent assay, found with `.primaryAssayName()` and visible in
   `.assayOwners(before, name)`, isn't in `.assayOwners(after, name)`, or
   whose root altExp (`.rowRoot()`, `.rootShape()`) is gone. Keep only
   the top ones: drop any whose parent subset is also affected.
2. `"error"`: stop with the message above. The object is unchanged,
   because R doesn't assign the result of a replacement function that
   errors.
3. `"remove"`: `after` without the top subsets and those below them,
   with the warning (a generalized `.dropSubsets()`).
4. `"keep"`: for each top subset, store `.subsetAssayView(before, name,
   primary)` as its own assay `primary` in `after`'s internal object,
   with the message. For an altExp-rooted subset whose container lost the
   assay, store it in the subset's own internal object (whose rows are
   the altExp's), before rows are followed.

`cause` gives the words for the messages ("the assay 'counts'", "altExp
'ADT'", and so on).

It replaces `.checkNotParents()`, `.checkAltExpParents()` and
`.stopUsedParents()`, and the drop step of `.carryRowRoots()` for paths
that take `dependents`. `.carryCombined()` (combineRows/combineCols)
keeps calling the drop step directly, with the "Removed" wording.

### Paths

- **New `assays<-` method**, in the `.esClasses` loop in
  `R/SubsetMethods.R`, forwarding to `.assays<-` in
  `R/ExperimentSubset.R`. It takes `dependents` from `...`, calls the
  next method without it, then `.resolveDependents(x, out, ...)`.
  `.assay<-` passes `dependents` through `callNextMethod(x, i, ...)`.
  Check during implementation that `dependents` is dropped before
  SummarizedExperiment's `assays()` getter, which `assay<-` also calls
  with `...`.
- **`setSubsetExperiment()`** gains `dependents`. `.checkNotParents()` is
  removed; `.resolveDependents()` runs at the end on `before` and the
  result. Its own assays that `value` lacks are still removed as today.
- **`altExp<-` and `altExps<-`** take `dependents` from `...`, with or
  without `subsetName`, through `.otherArgs()` and `.unlinkArg()` as
  `unlink` is today. `.applyAltExps()`, `.carryAltExps()` and
  `.carryRowRoots()` get `dependents`, so that the #58 drop of container
  parts and a removed altExp reach `.resolveDependents()`.
- **`removeSubset()`**: generic in `R/allGenerics.R`; one method in the
  `.esClasses` loop; the internal `.removeSubset()` in
  `R/ExperimentSubset.R`, with a `.reRootChild()` helper for `"keep"`.
  Rebuild subsets with `.copySubset()`.

### Edge cases

- Removing several assays at once, such as `assays(x) <- SimpleList()`,
  covers all affected subsets in one error, warning or set of messages.
- `"keep"` on a subset whose parent is a qualified name (`"s/scaled"`)
  stores `scaled` in the top affected subset, the child of `s`.
- A subset made from the removed assay that already stores an assay of
  that name: it isn't affected, since the assay stays visible to it.
- Objects saved before this change need no update: no slot changes.

## Docs

- New help page for `removeSubset()` with an example, its generic named
  in the `@eval .rdMethods(...)` block.
- `dependents` documented on the help pages for `assay<-`/`assays<-`,
  `setSubsetExperiment()`, `altExp<-` and `altExps<-`. `assays` is added
  to the `.rdMethods()` call of the `assay` page.
- `_pkgdown.yml`: add `removeSubset`; run `make site-check`.
- The vignette is unchanged; the CITE-seq article is unchanged, since it
  doesn't remove parent assays (to check with `make article`).
- NEWS, 1.23.57, marking the behavior changes:
  - Removing a main assay that subsets were made from now stops by
    default; before, it succeeded and left them broken.
  - Removing a whole altExp that subsets were made from now stops by
    default; before, it dropped them with a warning. `dependents =
    "remove"` gives the old behavior.
  - Giving the main altExp other rows, when that drops an assay an
    altExp-rooted subset was made from, now stops by default.
  - The `setSubsetExperiment()` error message changes.
  - The "Dropped subsets" warning now reads "Removed subsets".
- Version: 1.23.56 → 1.23.57 in DESCRIPTION.

## Testing

Test first, in `tests/testthat/test-ExperimentSubset.R`; each new test
must fail on its assertion, not on a missing function.

- **Main assays:** `assay(x, "counts") <- NULL` and `assays(x) <- ...`,
  each with `"error"` (object unchanged, message), `"keep"` (values
  identical to before, `subsetParent()` and `subsetSummary()` unchanged,
  subsets below unchanged, message), and `"remove"` (subsets and those
  below gone, warning). An invalid `dependents` is an error.
- **Unaffected subsets:** a subset made from `logcounts` stays and loses
  `counts`, with no error.
- **Not relinked:** after `"keep"`, adding a new `counts` to the main
  object doesn't change `assay(x, "s")`.
- **Subsets' own assays:** `setSubsetExperiment()` without `scaled`, with
  all three values. The existing tests of the old message are updated to
  the new one, and nothing else in them changes.
- **altExp:**
  - storing `clean`'s ADT without `clr`, with all three values;
  - giving the main ADT other rows, with all three values; with `"keep"`,
    the kept `clr` has the remaining features;
  - removing ADT: `"error"`, `"remove"` (the old drop behavior, with the
    new wording), and `"keep"` (the #127 error).
- **`removeSubset()`:**
  - a leaf;
  - several names at once;
  - an unknown name;
  - a parent with `"error"` and with `"remove"`;
  - `"keep"` on the main object → s → c, and p → s → c with p a subset:
    `c`'s values, row and column names, `subsetParent()`, linked assays
    still linked, and its own altExp parts;
  - removing two subsets of one chain;
  - an altExp-rooted child, both kept and refused.
- **Classes:** at least one test on a SummarizedExperiment-based object,
  since `assays<-` is defined for all five classes.
- **combineRows():** the existing tests still pass, with the "Removed"
  wording.

Then: `make test`, `make check`, `make bioccheck`, `make lint` and
`make site-check`, all clean. `make article` for the CITE-seq article, to
confirm it renders.
