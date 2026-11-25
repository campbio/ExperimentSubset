# Subsets follow a rename of the main assays they were made from (#126)

- **Status:** Approved (2026-11-25)
- **Date:** 2026-11-25
- **Branch:** `fix/assay-renames`, from `devel`
- **Lands:** one PR to `devel` (1.23.59); devel only, no release port.
- **ADR:** none; the object model doesn't change.

## Intent

`.assays<-` returns early when the main assays were only renamed
(`.renamedOnly()`), so subsets keep the old names as their parents. On
devel (f192b6f):

- the issue's example: after `assayNames(es)[1] <- "raw"`, `u`'s bare
  parent `"counts"` falls back to the one subset that stores a
  `"counts"`, `q`, and `assay(es, "u")` silently reads `q`'s copy;
- swapping two names makes a subset silently read the other assay;
- a subset made from `"p/logcounts"` errors once main `logcounts` is
  renamed;
- a main assay can be renamed to a subset's name, after which a parent
  naming that subset means the assay.

## Changes

1. `.assays<-` first runs `.qualifyParents(x)`, as other paths that
   store an assay do, so a bare stored-assay parent in an object saved
   before 1.23.15 can't be captured by a main assay renamed to its name.
2. When only names changed, old names map to new by position, all at
   once (so swaps work). Each subset whose rows come from the main
   object, read in the object before the change:
   - a bare parent `a` that was a main assay becomes its new name,
     unless the subset stores its own `a` (as `dependents = "keep"`
     leaves it, #122), which it keeps reading;
   - `"P/a"`, where `P` saw `a` linked from the main object, becomes
     `"P/<new>"`; where `P` stores its own `a`, it stays.
3. Errors, before anything changes:
   - a main assay can't take a subset's name, by rename or by
     `assays<-` adding one: "Assay names must differ from subset names;
     'p' is a subset.";
   - a subset that sees a renamed main assay linked and stores its own
     assay of the new name would hide the linked one: "Subset 'q' already
     stores an assay called 'raw'. Rename it first."

The fallback for bare stored-assay parents (objects saved before
1.23.15) stays as it is: objects record no version, and with the above
a current object no longer reaches it. Decided in chat: no
version-dependent behavior.

## Tests

In `tests/testthat/test-ExperimentSubset.R`, written first:

- the issue's example: `u` reads `m[1:3, 1:6]`, its parent is `"raw"`;
  `q` keeps its own copy and its parent `"counts"`;
- a swap: subsets follow their data;
- `"p/logcounts"` becomes `"p/lc"`; `"P/a"` where `P` stores `a` stays;
- `assays(x) <- setNames(...)` as `assayNames<-`;
- both errors;
- an object with a bare stored-assay parent, renamed into that name,
  still reads the subset's assay;
- subsets made from an altExp are unaffected.

## Docs, NEWS, version

One sentence in the `assay<-` details; `make docs`. NEWS entry; version
1.23.59.
