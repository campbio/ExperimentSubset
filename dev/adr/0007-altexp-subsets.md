# 0007. Subsets made from an altExp: an altExpName slot on AssaySubset

- **Status:** Accepted
- **Date:** 2026-11-21
- **Deciders:** Joshua D. Campbell

## Context

Every subset's rows index the main object's features, and every chain
of subsets ends at the main object. A subset can't hold only some of an
altExp's features, such as the ADT proteins without the isotype
controls (#92). The CITE-seq workflow article (#93) needs that.

To support it, a subset has to record that its rows index an altExp's
features rather than the main object's. Its columns still index the
main object's cells, through the chain of subsets it was made from. One
subset in a chain records this, and the subsets below it inherit it
through the chain.

Constraints:

- `AssaySubset` is an S4 class stored in users' saved objects. Adding a
  slot means objects saved before the change don't have it. In R,
  reading a missing slot is an error, and `validObject()` fails.
- `parentAssay` is a free string. altExp names can be anything,
  including a subset's name or a string containing `/`, the separator in
  `"subset/assay"`.
- The package already keeps one hidden marker in a subset's internal
  object: the names of altExp "pieces", in `int_metadata` (#58).
- Bioconductor's convention for class changes is an `updateObject()`
  method.

## Decision

We will add a slot `altExpName` (`NullOrCharacter`, prototype `NULL`) to
`AssaySubset`. It's set only on the subset made with
`createSubset(altExp =)` and holds the altExp's name.

The internal accessor `.rootAltExp()` returns `NULL` when the slot is
missing, so objects saved before the change keep working without
`updateObject()`. `updateObject()` methods for `AssaySubset` and the
five `Subset*` classes add the slot.

Every place that rebuilds an `AssaySubset` goes through one internal
helper that copies all slots, so a rebuild can't drop the new one.

## Consequences

- The record is explicit: the class shows what a subset holds, and the
  accessor sits with the others.
- Objects saved with 1.23.54 or earlier work as they are, and
  `updateObject()` brings them up to date. A test reads a small object
  saved before the change.
- The package gains `updateObject()` methods, which the S4 method
  snapshot test lists.
- A future slot follows the same pattern: a tolerant accessor, a line in
  `updateObject()`, and the copy helper.
- Internals can no longer assume that a subset's rows index the main
  object. Code that pairs row indices with `rownames(x)` or `rowData(x)`
  must use the subset's root object. AGENTS.md says so.

## Alternatives considered

- **A key in the internal object's `int_metadata`**, like altExp pieces.
  No class change and no `updateObject()`. Not chosen:
  - the record is hidden;
  - it survives only while nothing rebuilds the internal object, and
    `combineRows()` converts SingleCellExperiment objects to
    RangedSummarizedExperiment, dropping `int_metadata`;
  - it still changes the object model.
- **Encoding the altExp in `parentAssay`** (such as `"ADT::clean/clr"`).
  Fragile, since altExp names can contain any separator, as #92 notes.
- **Keeping an altExp's subsets inside the altExp**, by making the
  altExp itself an ExperimentSubset object. No class change. Not chosen:
  - a subset couldn't start from another subset's view of the altExp,
    which the article needs (its `decontPro` and `clr` are stored on the
    `clean` subset's ADT);
  - every function taking `subsetName` would have to search the altExps.
