# One parent-chain walk, and createSubset() without copying the parent (#36, part 2)

- **Status:** Approved
- **Date:** 2026-11-20
- **Branch:** `fix/dedup-chain-walks`, from `devel`
- **Lands:** one PR to `devel` (1.23.54); devel only, no release port
- **ADR:** none; only internal functions in `R/ExperimentSubset.R` change

## Intent

Part 1 of #36 (#118) defined the shared `Subset*` methods once. This is
the rest of #36, as far as it is still worth doing.

Since #36 was filed, #2, #20 and #28 have rewritten much of what it
describes:
- `.assay()` no longer walks the chain by recursion;
- `.getParentAssayName()` is gone. Name lookups no longer depend on
  first versus last: `.parentSubsetName()` stops when two subsets store
  a name, and `.qualifyParents()` rewrites bare names;
- every walk has a cycle guard.

What is left:

1. **Three walks up the parent chain.** `.subsetParent()`,
   `.subsetChain()` and `.toRootIndices()` each follow
   `.parentSubsetName()` upward in their own loop, each with its own
   cycle guard. The guards use two wordings. `.rootIndices(x, s, f)`
   computes the same thing as `.indicesIn(x, s, NULL, f)`.
2. **`createSubset()` copies the parent assay twice.** `.rows()` and
   `.cols()` each call `assay(x, parentAssay)`, which builds the whole
   assay with dimnames, only to read its size and row or column names.
   For a subset parent, it also slices the assay it links to.
3. **`.rows()` and `.cols()` are near copies**, 33 lines each.

Not changed, by the developer's choice:
- the `arglist`/`callNextMethod()` block repeated in 19 internals;
- `.subsetColnames()`/`.subsetRownames()` and their setters.

`callNextMethod()` has to be called in the method's own frame, so a
helper could replace only the `subsetName` lookup, saving about one
line per block. The colnames/rownames pairs are four three-line
functions, and merging them would only add indirection. #36 is closed
after this lands, with a comment saying why these two stay.

**Success:**
- the S4 method snapshot from #118 is unchanged;
- every existing test passes unchanged;
- the vignette output is unchanged;
- coverage doesn't drop;
- `createSubset()` no longer builds the parent assay.

## Design

### One walk: `.containerChain()`

```r
# Subsets that parentAssay is inside, from the outermost down to the one
# it names or points into; character(0) for a main assay. Each step moves
# up one subset, so a longer walk means a cycle. label names the start
# in the error.
.containerChain <- function(x, parentAssay, label = parentAssay) {
  chain <- character(0)
  for (step in seq_len(subsetCount(x) + 1)) {
    container <- .parentSubsetName(x, parentAssay)
    if (is.null(container))
      return(chain)
    chain <- c(container, chain)
    parentAssay <- .parentAssay(.subsets(x)[[container]])
  }
  stop("The parent chain of '", label, "' does not lead to a main assay.")
}
```

The other functions are built on it:

- **`.subsetChain(x, s)`** becomes
  `c(.containerChain(x, .parentAssay(.subsets(x)[[s]]), s), s)`. It
  gives the same result and the same error.
- **`.liftIndices(x, parentAssay, indices, getIndices, container = NULL,
  label = parentAssay)`** replaces `.toRootIndices()`. It maps indices
  relative to `parentAssay` up the chain, stopping at `container` (a
  subset name), or at the main object when `container` is `NULL`.
  `.createSubset()` calls it for the new subset's root rows and
  columns.
- **`.indicesIn(x, subsetName, getIndices, container = NULL)`** is one
  call into `.liftIndices()`, starting from the subset's own indices
  and its `parentAssay`. The argument order changes so that `container`
  can default to `NULL`.
- **`.rootIndices()` is removed.** Its 12 call sites become
  `.indicesIn(x, name, .rowIndices)` or `.indicesIn(x, name,
  .colIndices)`. The two calls in `.subsetAssayView()` pass
  `container =`.
- **`.subsetParent()`** keeps its checks and its output exactly: `NULL`
  for a subset made with no parent, otherwise a list of names. It builds
  that list by walking `.subsetChain(x, ref$subset)` from the subset
  upward, instead of following `.parentSubsetName()` itself. At each
  step, it adds:
  - the `parentAssay` when it differs from the container (a
    `"P/assay"`, or a legacy bare name);
  - then the container;
  - for the outermost subset, its `parentAssay` (a main assay, or a
    name that resolves to nothing, as today).

  For a `"subset/assay"` input, the list starts with the subset, as
  today.

**Error wording:** one guard means one wording: "The parent chain of
'<label>' does not lead to a main assay." `.toRootIndices()` said "of
the subset"; `.liftIndices()` in `createSubset()` now names the
`parentAssay`. For a `"S/assay"` input, `subsetParent()` names `S`, as
on devel: `.resolveAssay()` walks the chain of `S` first, so the error
comes from there (corrected during review). Only an object whose subsets form a cycle reaches this error, and no
existing test checks its text.

`.primaryAssayName()`, `.assayOwners()`, `.nearestOwner()`,
`.altExpHolders()`, `.subsetsBelow()` and `.remapSubsets()` keep calling
`.subsetChain()` and don't change.

### One index function, without the copy: `.selectIndices()`

`.rows()` and `.cols()` become
`.selectIndices(x, select, parentAssay, dim)`, where `dim` is `"row"`
or `"col"`. The words in each message come from `dim`: `'rows'`/`'cols'`,
`rows`/`columns` and `row`/`column`. Every message stays word for word:

| Check | Message (rows; columns the same) |
|---|---|
| parent has none | Parent '<p>' has no rows, so a subset can't be made from it. |
| names, but parent has none | 'rows' are given as names, but the parent assay has no row names. |
| out of range, NA, or too many | Selected rows not available in the specified assay. |
| empty | 'rows' must select at least one row. |
| not whole or below 1 | 'rows' must be positive whole numbers or row names. |
| duplicates | 'rows' must not contain duplicates. |

The checks keep their order, so the same input gives the same error.

**Size and names without building the parent:**
- **A main assay:** the main object's `nrow(x)`/`ncol(x)` and
  `rownames(x)`/`colnames(x)`. `assay(x, i)` returns exactly these,
  because it gives the assay the object's dimnames.
- **A subset, or `"subset/assay"`:** the subset's internal object, which
  has the subset's dimensions and names.
- **Fallback:** names are read only when `select` is character. If the
  subset's internal object then has no names of that dimension, names
  come from `assay(x, parentAssay)`, as today.

  That case is real. After `rownames(x) <- NULL`, the assays stored in
  `x` keep their row names. A subset made with no names then shows them
  through `assay()`, and `rows = "g1"` matches them. The fallback keeps
  that working, and copies only in this case.

`parentAssay` comes from `.resolveAssay()`, which `.createSubset()`
already calls when a `parentAssay` is given.

**No assays:** the default parent of an object with no assays is
`assayNames(x)[1]`, which is `NULL`. Today `createSubset()` then fails
with an S4 dispatch error: "unable to find an inherited method for
function 'assay' for signature 'x = …, i = "NULL"'". The new code
stops first, with "'x' has no assays to make a subset from.", and a new
test checks it. This is the only change to an error that a user can
reach.

## Tests

Tests are written first where the change can be observed, and run
against the current code before it changes.

1. **The parent assay isn't built** (fails today). Results are computed
   first. Then `local_mocked_bindings(.assay = function(...) stop(...))`
   makes any `assay(x, ...)` on the object stop. `createSubset()` from a
   main assay, from a subset and from a `"subset/assay"` must still
   succeed, with the same row and column names and indices as the
   unmocked run.
2. **Names stored only on the assay** (passes today; guards the
   fallback). After `rownames(x) <- NULL`, a subset with no row names,
   and a subset of it selected by the stored row names, gives the same
   indices as selecting by position.
3. **No assays** (fails today on the message). `createSubset()` on an
   object with no assays stops with "has no assays".
4. **Cycle guard** (passes today; the path had no test). An object whose
   subsets name each other as parents is built with the internal setters
   (`.subsets<-`, `AssaySubset()`). Then `subsetParent()`, `assay()` of
   the subset, and `createSubset(parentAssay =)` from it each stop with
   "does not lead to a main assay" instead of hanging.
5. **Unchanged:** the S4 method snapshot from #118, the 10
   `subsetParent()` expectations, the `.rows()`/`.cols()` message tests,
   and the rest of the suite.

## Housekeeping

- `NEWS.md`, 1.23.54: "`createSubset()` no longer builds the parent
  assay to check `rows` and `cols`, and stops with a clear message when
  there is no named first assay. The parent-chain walks share one internal
  function. No change in results (#36)."
- `DESCRIPTION`: version 1.23.54. Keep `Config/roxygen2/version` at
  8.1.0.
- `AGENTS.md`: only if it names a removed helper. It lists
  `.subsetParent()`, which stays.
- **Checks:**
  - `make test`;
  - `make coverage` (must not drop from devel);
  - `make check-full`, `make bioccheck`;
  - `make lint`, `make site-check`;
  - the built vignette compared with devel's.
- **Review:** a spec review by a fresh subagent, then `/code-review`.
- **After merge:** close #36 with a comment on the two pieces left as
  they are.

## Found during implementation

- **`.selectIndices()` takes `ref` too.** It is `.resolveAssay(x,
  parentAssay)`, which `.createSubset()` already has, so the parent isn't
  resolved again. The size and names come from a small helper,
  `.parentExtent()`.
- **No assays, or an unnamed first assay.** The error is "'parentAssay'
  is needed: 'x' has no named first assay to make the subset from." It
  covers three cases:
  - no assays;
  - assays that all lost their names (`assays(x) <- list(m)`, where
    `assayNames(x)[1]` is `NULL`);
  - an unnamed first assay (`""`).

  For `""`, devel failed in `assay()` with "invalid subscript". Without
  this check, the branch would have made a subset with
  `parentAssay = ""`. The spec review found that case. One test covers
  all three.
- **A parent subset that lost its assay.** Because `devel` built
  `assay(x, parentAssay)`, it stopped with "'counts' is not an assay of
  subset 's1'." when the assay a parent subset was made from had been
  removed (`assay(es, "counts") <- NULL`). Without that call, the branch
  made the subset and failed only later, at `assay()`.

  `.createSubset()` now makes the same check with the same message,
  through `.primaryAssayName()` and `.assayOwners()`, without copying.
  A new test covers it, and the spec review found this case too.

  An old subset with no `parentAssay` reads its own "counts". When its
  internal object has none, the check now stops with "'counts' is not
  an assay of subset 's2'." `devel` stopped with SummarizedExperiment's
  subscript error instead.

  Without the check, the branch made the child, and `assay()` of the
  child then quietly read the main object's counts. `/code-review`
  found this.

  `.checkSeesAssay()` holds the rule and its message for both
  `createSubset()` and `assay()`.
- **Review findings left as they are:**
  - Repeated chain walks inside `createSubset()`: a few lookups on short
    chains, small next to the copies removed.
  - A no-names case in the mock test: that fallback builds the assay on
    purpose, and test 2 checks its result.
  - A clear error for a `container` not in the chain: no caller passes
    one, and devel behaves the same.
- **`force(label)` in `.containerChain()`.** The default
  `label = parentAssay` is evaluated lazily, so without `force()` the
  error would name the subset the loop had reached, not where the walk
  started. No caller relies on the default today. The spec review found
  this.
- **Two more `.indicesIn()` calls.** `.subsetAltExpView()` calls it too,
  besides `.subsetAssayView()`. Both calls now use the new argument
  order.
- **Old and new code compared.** A scratch script compared the devel
  versions of the changed functions with the new ones on 60 random
  objects (chains of 3-7 subsets, parents that are subsets or
  `"subset/assay"`, rows and columns by position or name, every fifth
  object with its row names removed). It found the same results in 1,910
  comparisons:
  - the subset chain;
  - indices up to the main object and to every container;
  - `subsetParent()` for every subset and `"subset/assay"`;
  - the stored subsets after `createSubset()`;
  - the error for 7 kinds of bad `rows` or `cols`.
