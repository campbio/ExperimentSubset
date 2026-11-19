# Define the Subset* methods once, in loops over the classes (#36, part 1)

- **Status:** Approved
- **Date:** 2026-11-19
- **Branch:** `fix/dedup-subset-methods`, from `devel`
- **Lands:** one PR to `devel` (1.23.53); devel only, no release port
- **ADR:** [0006](../adr/0006-define-subset-methods-in-loops.md)

## Intent

Issue #36: the five `R/Subset<Class>.R` files hold 2,706 lines and 255
`setMethod()`/`setReplaceMethod()` calls, but only about 70 distinct
methods. Each is a thin wrapper around a `.helper()` in
`R/ExperimentSubset.R`, written out once per class:

| Set | Methods | Classes |
|---|---|---|
| Base: metadata, assay, createSubset, show, subset*, rowData/colData, `[`, cbind/rbind, combineCols/combineRows, ... | 37 | all five |
| altExp/altExps/altExpNames and reducedDim/reducedDims, getters and setters | 19 | SCE, TreeSE, SpatialExperiment |
| subsetSpatialCoords, subsetImgData, spatialCoords, subsetSpatialData, spatialData | 7 | SpatialExperiment |
| subsetRowLinks/ColLinks, subsetRowTree/ColTree | 6 | TreeSE |

A change to a method has to be made up to five times, and the copies
drift: the `altExp<-` bug in #21 was in three files.
`SubsetSummarizedExperiment.R` and `SubsetRangedSummarizedExperiment.R`
differ only in the class name and the order of their methods.

The goal is that each wrapper is written once, with no user-visible
change: the same exported names, the same S4 methods on the same
classes, the same dispatch, and the same help pages and aliases.

This spec covers only the class files. The rest of #36 gets its own
specs and PRs:
- the `arglist`/`callNextMethod` block repeated about 20 times in
  `R/ExperimentSubset.R`, and the near-copies `.subsetColnames()`/
  `.subsetRownames()` and `.rows()`/`.cols()`;
- the issue comment: one parent-chain resolver, slicing a chained
  subset once instead of copying at each level, and the first-versus-
  last name lookup. That changes internal behavior and performance, so
  it isn't a pure refactor.

## Changes

### Class vectors (`R/allClasses.R`)

After the `ExperimentSubset` union, two internal character vectors:

- `.esClasses`: the five `Subset*` class names, in the union's order.
- `.esSCEClasses`: `"SubsetSingleCellExperiment"`,
  `"SubsetTreeSummarizedExperiment"`, `"SubsetSpatialExperiment"`, the
  classes that extend SingleCellExperiment.

### Shared methods (`R/SubsetMethods.R`, new)

- `for (cls in .esClasses) { ... }` holds the 37 base methods.
- `for (cls in .esSCEClasses) { ... }` holds the 19 altExp and
  reducedDim methods.
- Each `setMethod()`/`setReplaceMethod()` call keeps today's generic,
  signature (with `cls` in place of the class name), formals and body,
  copied unchanged. This includes the four `colData<-` value
  signatures (DataFrame, NULL, data.frame, ANY) and their #99 comment,
  reworded so it doesn't name SummarizedExperiment or
  SpatialExperiment only.
- The `@importMethodsFrom SummarizedExperiment cbind rbind` and
  `@importMethodsFrom S4Vectors combineCols combineRows` tags move with
  their methods.
- Before the old files are deleted, a scratchpad script compares every
  base method's body across the five files, and every altExp/reducedDim
  body across the three, after replacing the class name. Today's diffs
  show only whitespace and line-wrapping differences. If any body
  differs in substance, implementation stops and the difference goes
  back to the maintainer, because choosing one copy would change
  behavior for some class.

### Class-specific methods

- `R/SubsetSpatialExperiment.R` keeps only its 7 SpatialExperiment
  methods, and `R/SubsetTreeSummarizedExperiment.R` its 6 TreeSE
  methods, written out by hand as today; each applies to one class.
- `R/SubsetSummarizedExperiment.R`, `R/SubsetRangedSummarizedExperiment.R`
  and `R/SubsetSingleCellExperiment.R` are deleted (approved with this
  spec).
- DESCRIPTION `Collate:`: `allClasses.R`, `allGenerics.R`,
  `ExperimentSubset.R`, `SubsetMethods.R`, `SubsetSpatialExperiment.R`,
  `SubsetTreeSummarizedExperiment.R`.

### Documentation (`.rdMethods()`)

roxygen documents a method from the `setMethod()` call that follows its
block, and it can't see calls inside a loop. Without help, R CMD check
reports every looped method as an undocumented S4 method.

- A new internal helper, `.rdMethods(generics)`, in
  `R/SubsetMethods.R`, returns roxygen tag lines for every method of
  `generics` whose first signature class is in `.esClasses`:
  - an `@aliases` line with `<generic>,<sig1>[,<sig2>...]-method` for
    each, as roxygen writes them today (for example
    `colData<-,SubsetSummarizedExperiment,DataFrame-method`);
  - an `@usage` line with an `\S4method{<generic>}{<signature>}(<formals>)`
    entry for each, from the registered method's formals, with
    defaults deparsed (`withDimnames = TRUE`) and replacement methods
    written as `\S4method{colData}{...}(x, ...) <- value`, as roxygen
    writes them today.
  - Methods are listed by class in `.esClasses` order, then in the
    order of `generics`. The signature is the one given to
    `setMethod()` (the method's defined signature), not the one padded
    to the generic's full signature, as roxygen writes it today:
    `rowData<-,SubsetSummarizedExperiment-method`, but
    `colData<-,SubsetSummarizedExperiment,ANY-method` where `"ANY"` was
    written. The helper reads it through the methods package's
    functions, not `@`.
- roxygen loads the package with `load_all()` before it reads the
  blocks, so the methods exist when `@eval` runs. The helper reads them
  through `methods::findMethods()`/`methods::getMethod()`, so the docs
  follow the loops: a class added to `.esClasses` is documented without
  further edits.
- Each Rd page that documents looped methods gets one code-free block
  in `R/SubsetMethods.R`, next to its loop:
  ```r
  #' @rdname ExperimentSubset-rowData-colData
  #' @eval .rdMethods(c("rowData", "colData", "rowData<-", "colData<-"))
  NULL
  ```
  The Rd names, titles, descriptions, `@param`s and examples stay where
  they are now (`R/allGenerics.R`, `R/ExperimentSubset.R`).
- The Spatial and Tree files keep their ordinary `@rdname` blocks.
- `make docs` regenerates `man/`. Expected changes are the order of
  aliases and usage entries (by `.esClasses` instead of by file) and the
  "Please edit documentation in" header lines. A scratchpad script
  compares the sorted `\alias{}` lines and the sorted `\S4method` usage
  entries of every Rd file before and after; they must match exactly.
  If roxygen wraps usage lines, the comparison ignores whitespace.
- NAMESPACE must not change, apart from the order of import lines.

### AGENTS.md

- Layout: replace the `R/Subset<Class>.R` bullet with `R/SubsetMethods.R`
  (the shared wrappers, looped over `.esClasses`/`.esSCEClasses`, plus
  `.rdMethods()`), and the two class-specific files. Replace "A change
  to a method usually needs the same change in all five files" with
  where a new method goes and that its Rd page needs an `@eval
  .rdMethods()` block.
- Object model: the class list stays; add that a new `Subset*` class
  goes into `.esClasses` (and `.esSCEClasses` if it extends
  SingleCellExperiment) as well as the union.

### NEWS and version

- DESCRIPTION: 1.23.52 to 1.23.53.
- NEWS.md, under `# ExperimentSubset 1.23.53 (<date>)`: the methods
  shared by the `Subset*` classes are now defined once, in loops over
  the classes, in `R/SubsetMethods.R`, rather than copied into one file
  per class; `SubsetSummarizedExperiment.R`,
  `SubsetRangedSummarizedExperiment.R` and `SubsetSingleCellExperiment.R`
  are gone. No change in behavior (#36).

### ADR

- `dev/adr/0006-define-subset-methods-in-loops.md` (status Proposed)
  and its row in `dev/adr/README.md`, proposed with this spec. The
  maintainer accepts it before implementation.

## Testing

In `tests/testthat/test-ExperimentSubset.R`:

1. **Method inventory snapshot, written and committed first, on today's
   code.** `expect_snapshot()` of the sorted
   `<generic>,<signature>` strings of every method in the package's
   namespace whose first signature class is a `Subset*` class. This is
   a characterization test: it passes on today's code by design rather
   than failing first, because a refactor has no new behavior to fail
   on. Its job is to fail if the refactor drops, adds or retypes a
   method. The snapshot file (`tests/testthat/_snaps/`) is committed
   with it and doesn't change during the refactor.
2. `.esClasses` has the same members as the `ExperimentSubset` union,
   and every `.esSCEClasses` class is in `.esClasses` and extends
   SingleCellExperiment. This fails first because `.esClasses` doesn't
   exist yet (an object-not-found error rather than an assertion; the
   vectors are new).
3. `.rdMethods()` returns, for `"subsetDim"`, the alias
   `subsetDim,SubsetSummarizedExperiment-method` and the usage entry
   `\S4method{subsetDim}{SubsetSummarizedExperiment}(x, subsetName)`;
   for `"colData<-"`, the alias and the replacement usage of the
   DataFrame method; for `"reducedDim<-"`, a usage entry with
   `withDimnames = TRUE` and only the three SCE-based classes.

The existing suite, unchanged, covers the behavior of every method.

Then `make test`, `make coverage` (must not drop: the deleted lines
were covered copies, and `.rdMethods()` is covered by test 3),
`make lint`, `make docs`, `make site-check`, and before a PR
`make check-full` and `make bioccheck`. R CMD check's codoc and
undocumented-methods checks confirm the generated docs.

## Out of scope

- Methods defined on the `ExperimentSubset` union (ADR 0005 stands).
- The `arglist`/`callNextMethod` duplication and the near-copy helpers
  in `R/ExperimentSubset.R`.
- The parent-chain resolver, slice-once copying and name lookup from
  the issue comment.
- Any change to a method's signature, formals or body.

The PR says "Part of #36"; the issue stays open for the other two
pieces.

## Found during implementation

- `R/SubsetTreeSummarizedExperiment.R` had `@usage NULL` on its
  `subsetColData<-` and `subsetRowData<-` methods (#8), added when
  their three-class signatures gave usage lines too wide for the PDF
  manual. Their signatures are now `x` only, so `.rdMethods()` writes
  their usage lines like the other classes': 257 `\S4method` entries
  instead of 255, with all 312 aliases unchanged. `make check-full`
  reports no usage-width NOTE.
- S4 stores a method with arguments its generic lacks (`combineCols()`
  and `combineRows()`: `delayed`, `fill`) wrapped in `.local()`, so
  `.rdUsage()` reads the formals through `methods::unRematchDefinition()`,
  as roxygen does.
- `.rdMethods()` stops when a named generic has no methods for the
  classes, so a misspelled generic can't silently drop docs (from the
  spec review).
