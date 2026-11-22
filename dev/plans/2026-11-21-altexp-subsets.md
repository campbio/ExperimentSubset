# Subsets made from an altExp: their own features and cells (#92)

- **Status:** Approved
- **Date:** 2026-11-21
- **Branch:** `feature/altexp-subsets`, from `devel`
- **Lands:** one PR to `devel` (1.23.55); devel only, no release port.
  The CITE-seq article section is a follow-up PR.
- **ADR:** [0007](../adr/0007-altexp-subsets.md), Accepted.

## Intent

Today every subset's rows are the main object's features. A subset can't
hold only some of an altExp's features, such as some ADT proteins (#92).

The CITE-seq article (#93) drives this. Its follow-up section will:
- make a subset of the high-quality cells without the three isotype
  controls, which sees ADT `counts`, `decontPro` and `clr`;
- show that renaming the main ADT, or dropping the isotype controls from
  it, carries along to that subset.

In the article, `decontPro` and `clr` are stored on the `clean` subset's
view of ADT, not on the main ADT. So the subset has to be able to start
from a subset's view of an altExp, not only from the main altExp.

Decided in chat (2026-11-21):
- **Starting point:** a subset made with `altExp =` can start from the
  main object's altExp or from a subset's view of it.
- **Class:** its internal object has the class of the altExp it starts
  from, so a subset of a SingleCellExperiment altExp is a
  SingleCellExperiment. Subsets made from it keep that class, as subsets
  of subsets already do.
- **Recording:** a new slot on `AssaySubset`, not a hidden
  `int_metadata` key or an encoded `parentAssay`.
- **When that altExp gets other rows:** the subset follows its features
  by name, as `x[i, ]` does today.
- **When that altExp disappears:** the subset and those made from it are
  dropped, with a warning.
- **The article's new section:** a follow-up PR.
- **Trees and spatial data:** row-tree accessors stop; column-side
  readers (coordinates, images, column trees) work through the column
  chain, since its cells are the main object's; column-side setters
  work only when the subset's class can hold the data.

**Success:**
- the article's case works on a small simulated object:
  - an ADT subset of `clean`'s cells without the isotype controls sees
    `counts`, `decontPro` and `clr`;
  - it stores its own assays and a PCA;
  - dropping the isotype controls from the main ADT, or renaming the
    main ADT, carries along to it;
- objects saved before this change still work, with and without
  `updateObject()`;
- the existing tests pass unchanged, apart from the S4 method snapshot,
  which gains the new `updateObject` methods;
- the vignette output is unchanged, and coverage doesn't drop.

## Interface

`createSubset(x, subsetName, rows = NULL, cols = NULL, parentAssay =
NULL, altExp = NULL)`. `altExp` goes last in the generic's formals, so
existing positional calls don't change. The `signature` stays `"x"`.

`altExp = "ADT"` is for the three classes that extend
SingleCellExperiment. On the other two it stops with "'altExp' needs an
object with altExps, such as a SubsetSingleCellExperiment; 'x' is a
SubsetSummarizedExperiment." `altExp` must be `NULL` or one non-empty
string.

With `altExp = "ADT"`, `parentAssay` names ADT's assays, in the grammar
the package already uses:

| `parentAssay` | starts from | `cols` relative to |
|---|---|---|
| `NULL` | the main ADT, its first assay | the main object's cells |
| `"counts"` | the main ADT's `counts` | the main object's cells |
| `"clean"` | `clean`'s view of ADT, its first assay | `clean`'s cells |
| `"clean/clr"` | `clean`'s view of ADT, its `clr` | `clean`'s cells |

- **`rows`** picks ADT's features as the starting point sees them, by
  position or name. A subset sees all of an altExp's features, so these
  are the main ADT's features, unless `clean` stores an ADT with other
  rows.
- **What is stored:**
  - the new slot, `altExpName = "ADT"`;
  - `parentAssay`, resolved to an explicit assay: `"counts"` or
    `"clean/clr"`, never a bare `"clean"`, so that reordering ADT's
    assays later can't change what the subset points to.
- **Errors**, checked before anything is stored:
  - with no `parentAssay`, the main object has no altExp called `altExp`:
    "'ADT' is not an altExp of the main object.";
  - with a subset as the parent, that subset doesn't see the altExp:
    "'ADT' is not an altExp of subset 'clean'.";
  - the assay isn't one of the altExp's assays there: "'clr' is not an
    assay of altExp 'ADT'." or "... of altExp 'ADT' of subset 'clean'.";
  - the parent is itself a subset made from an altExp: "Subset 'x' is
    made from altExp 'ADT', so a subset of it can't use 'altExp'; leave
    'altExp' out to subset its rows.";
  - the `rows`/`cols` checks and messages of `.selectIndices()`, worded
    as they are now.
- **Subsets made from it:** `createSubset(es, "x", parentAssay =
  "clean.antibodies", rows = 1:5)`, with no `altExp`, picks rows of
  `clean.antibodies`. `parentAssay = "clean.antibodies/clr"` works the
  same way.

After that it is an ordinary subset, whose root is the altExp, not the
main object:

- **Rows** index the altExp's features. The row chain stops at the
  subset made with `altExp`.
- **Columns** are relative to the parent, and the column chain goes all
  the way up: `clean.antibodies` → `clean` → `filtered` → the main cells.
- **Assays it sees:**
  - every assay of the altExp as its starting point sees it (for
    `clean`: `counts` linked from the main ADT, and `decontPro`/`clr` from
    `clean`), linked, cut to its rows and columns;
  - then its own, stored with `setSubsetAssay()` and the usual `unlink`
    rules. Messages name the source as "altExp 'ADT'" or "altExp 'ADT' of
    subset 'clean'".
- **Its own data:** reducedDims, metadata and altExps, as for any subset.
  It doesn't see the main object's altExps (there is no ADT inside ADT),
  only those stored by the subset that starts the root and the subsets
  below it.
- **Its internal object** is built by cutting the starting altExp to the
  subset's rows and columns and emptying its assays, rowData, colData,
  reducedDims, altExps and metadata. That keeps the altExp's class
  exactly, subclasses included. A subset of a plain SummarizedExperiment
  altExp can't hold reducedDims, as for a subset of a
  `SubsetSummarizedExperiment` today.
- **Tree and spatial accessors**, for these subsets and those below them,
  on a SpatialExperiment or TreeSummarizedExperiment main object:
  - **Row side: stop.** `subsetRowTree()`, `subsetRowTree<-` and
    `subsetRowLinks()` stop with "Subset 'x' is made from altExp 'ADT';
    its rows aren't the main object's, so the row tree doesn't apply."
    The main object's row tree links genes, and reading it at ADT
    positions would give the wrong features without an error.
  - **Column side: read as for any subset.** `subsetColTree()`,
    `subsetColLinks()`, `subsetSpatialCoords()`, `spatialCoords()` with
    `subsetName`, `subsetImgData()`, and the deprecated
    `subsetSpatialData()`/`spatialData()` with `subsetName` all work from
    `.indicesIn(x, s, .colIndices)`. The column chain ends at the main
    object's cells, so they return the main object's data at the subset's
    cells, or what the nearest subset up the column chain stores (such as
    coordinates or a column tree on `clean`).
  - **Column-side setters: only when the class can hold it.**
    `subsetSpatialCoords<-`, `subsetColTree<-`, `subsetSpatialData<-`, and
    `setSubsetExperiment()` storing changed coordinates or trees, work
    when the subset's internal object can hold them (ADT is itself a
    SpatialExperiment or TreeSummarizedExperiment). Otherwise they stop
    with "Subset 'x' is a SingleCellExperiment, which can't store
    spatialCoords." (or "a colTree"), not with an S4 dispatch error.

**Other functions with `subsetName =`:**
- **`rowData(parentRowData = TRUE)`:** the altExp's rowData, as the root
  sees it, at the subset's rows.
- **`colData(parentColData = TRUE)`:** the main object's colData at its
  cells, as today.
- **`getSubsetExperiment()`:** starts from the root's altExp object, not
  `x`. It's cut to the subset's rows and its columns relative to the
  root's container. Then it gets the assays the subset sees (its primary
  assay first), rowData, colData, and its own reducedDims, altExps and
  metadata. The class is the altExp's. The spatial and column-tree steps
  run only when that class can hold them; the row-tree step never runs.
  `setSubsetExperiment()` works through the same pieces.
- **`subsetParent()`:** names the altExp root as `"altExp/assay"`, the
  form `unlink` already uses, prefixed with the container when there is
  one. Then it continues up the column chain as today:
  - `list("clean/ADT/clr", "clean", "filtered/decontaminated",
    "filtered", "counts")` for `clean.antibodies`;
  - `list("ADT/counts")` for a subset of the main ADT;
  - for a subset made from `clean.antibodies` with
    `parentAssay = "clean.antibodies"`: `"clean.antibodies"`, then the
    list above.

  `subsetSummary()`'s Parent column joins these with " -> ", as now.
- **Left out (YAGNI):**
  - an altExp stored only on a subset as a root (it must be the main
    object's altExp, as the main object or a subset sees it);
  - `altExp` in `ExperimentSubset(subset = list(...))`;
  - logical `rows`/`cols`.

## Design

### Recording: the `altExpName` slot (ADR 0007)

- **The slot:** `AssaySubset` gains `altExpName = "NullOrCharacter"`,
  with prototype `NULL`, and `AssaySubset()` gains
  `altExpName = NULL`. Only the subset made with `altExp =` has it set;
  subsets below it find the root through the chain.
- **The accessor:** `.altExpName(subset)` in `R/ExperimentSubset.R`
  returns `NULL` when the slot is missing (objects saved before the
  change: `methods::.hasSlot()`), so old objects work without
  `updateObject()`.
- **`updateObject()`:**
  - a method for `AssaySubset` rebuilds an object without the slot, with
    `altExpName = NULL`;
  - a method for each class in `.esClasses`, defined in the loop in
    `R/SubsetMethods.R`, calls `callNextMethod()`, then `updateObject()`
    on each subset;
  - both are documented on one help page, with `@eval .rdMethods(...)`
    for the looped ones.
- **One copy helper:** `.copySubset(subset, ...)` returns `subset` with
  the named slots replaced, through `AssaySubset()`, copying the rest,
  `altExpName` included. It replaces the field-by-field
  `AssaySubset(...)` rebuilds in:
  - `.qualifyParents()`;
  - `.remapSubsets()`;
  - `.combineSubsets()`.

### Where a subset's rows come from

- **`.rowRoot(x, subsetName)`** walks `.subsetChain(x, subsetName)` from
  the subset upward and stops at the first subset with an `altExpName`.
  It returns:
  - `NULL` when there is none, which means the main object, as today;
  - otherwise `list(start =, altExp =, container =)`: the subset that
    starts the root, the altExp's name, and the subset its `parentAssay`
    points into (`NULL` for the main object).
- **`.rootObject(x, root)`** returns the object whose rows the subset
  indexes:
  - `x`, for `NULL`;
  - otherwise `altExp(x, root$altExp)` (main) or
    `.subsetAltExpView(x, root$container, root$altExp)`, each with
    `withDimnames = FALSE` as far as the callers need it.

### The chain walk

- **`.parentSubsetName(x, parentAssay, altExp = NULL)`:** with an
  `altExp`, a bare name is an assay of the altExp, so it never matches a
  main assay or the legacy rule for a stored-assay name. Only `"C"` or
  `"C/assay"` for an existing subset gives a container; anything else
  ends the chain.
- **`.containerChain()`** passes each step's `.altExpName()` to it.
  `.liftIndices()` gets an `altExp` argument for the first step, since
  `createSubset()` calls it before the subset exists.
- **`.qualifyParents()`** skips subsets with an `altExpName`.
- **Columns:** the column chain doesn't change, so
  `.indicesIn(x, s, .colIndices)` doesn't change.
- **Rows:** `.indicesIn(x, s, .rowIndices)` and `.liftIndices()` for
  rows stop at the subset that starts the root, so `container = NULL`
  means "in the root's rows". With a `container` below the root, they
  stop there, as now.
- **Callers that pair those rows with the main object's** `rownames(x)`,
  `rowData(x)` or `nrow(x)` use `.rootObject()` instead:
  - `.createSubset()`, for row names;
  - `.parentExtent()`;
  - `.getSubsetExperiment()`;
  - `.subsetRowData()`.

### Assays seen

- **`.assayOwners()`** starts from the root object's assay names; `""`
  still means "the root".
- **`.subsetAssayView()`** reads an assay owned by `""` from the root
  object. It uses rows in the root's rows and columns relative to the
  root's container (`.indicesIn(x, s, .colIndices, root$container)`).
- **`.primaryAssayName()` and `.resolveAssay()`** follow from those.
  `"clean.antibodies/clr"` resolves as an assay that
  `clean.antibodies` sees.
- **`.linkedSources()`** names the root, for an altExp root, as
  "altExp 'ADT'" or "altExp 'ADT' of subset 'clean'", through
  `.altExpSource()`.
- **`.altExpHolders()`** walks the chain only from the subset that starts
  the root.

### Trees and spatial data

- **Owner lookups:** `.nearestOwner()`'s checks (in `.treeOwner()` and
  `.spatialCoordsOwner()`) skip internal objects whose class can't hold
  a tree or coordinates, so a SingleCellExperiment internal object is
  passed over and the walk goes on up the column chain. The spatial
  check already wraps its test in `isTRUE()`. The tree check
  (`length(tree(internal, NULL)) > 0`) gets a class test first.
- **Row trees:** `.subsetTree()`/`.subsetLinks()` with `dim = "row"`
  and `subsetRowTree<-` stop for a subset with an altExp root, before any
  lookup.
- **Positions:** the column-side readers already use only column
  indices, so they don't change.
- **Setters:** `.subsetSpatialCoords<-`, `.subsetColTree<-` (through
  `.replaceTree()`), `.storeChangedSpatialCoords()` and
  `.storeChangedTrees()` check that the internal object is a
  SpatialExperiment or TreeSummarizedExperiment, and stop with the
  message in Interface otherwise.

### Building the subset

For an altExp root, `.createSubset()`:
1. checks `altExp` and resolves `parentAssay` to the container and the
   explicit assay;
2. reads the size and names from the root object (`.parentExtent()`)
   or, for `cols` relative to a container, from the container;
3. builds the internal object: `.rootObject()` cut to the new subset's
   rows and its root columns, with assays, rowData, colData,
   reducedDims, altExps and metadata emptied;
4. stores the subset with `altExpName`.

A subset made from an altExp-rooted parent, without `altExp`, builds its
internal object the same way from its parent's: the parent's internal
object, cut and emptied. So it keeps the class. A main-rooted subset is
built as today.

### Changes that carry along

**`.carryRowRoots(before, after, holder, renamed)`** runs after
`.carryAltExpParts()`, in each place that already calls it:
- `altExp<-`, `altExps<-` and `altExpNames<-`, on the main object and
  with `subsetName`;
- `setSubsetExperiment()`;
- `.carryCombined()`.

It runs after `.carryAltExpParts()`, because dropping `clean`'s part of
ADT can change the rows `clean` sees. For each subset that starts an
altExp root, parents first:

1. **Renamed.** The slot takes the new name when `holder` is the main
   object, the root's container, or a subset the container was made
   from.
2. **Gone.** If the container no longer sees that altExp, the subset and
   every subset made from it are removed. The warning: "Dropped subsets
   'clean.antibodies', 'x': the altExp 'ADT' they were made from was
   removed." If `clean` stores its own whole ADT, that one is kept, so
   the root still sees an ADT and goes on to step 3.
3. **Other rows.** If the root object's rownames differ before and after,
   the subset keeps the features whose names remain, in its own order.
   Its internal object is cut to them, and the subsets below follow, as
   for `x[i, ]`:
   - no rownames on either side and the same number of rows: nothing
     changes;
   - otherwise, rows that can't be matched are an error, raised before
     anything is stored, so the setter changes nothing: "Subset 's' can't
     follow the new rows of altExp 'ADT': they have no names to match.";
   - a subset left with no rows is kept with none.
4. **Other class.** Each internal object down that chain is converted to
   the new class with `as()`.

**The remapping helper.** `.remapSubsets()`'s loop moves into
`.remapFrom(x, out, keptRows, keptCols)`, which starts from given kept
positions for some subsets. `x[i, j]` seeds it from the main object's
`i` and `j`, and step 3 seeds it from the root's kept features.

**`x[i, j]`:** `i` never cuts the rows of an altExp root. Its kept rows
are `NULL`, since `[` doesn't change altExp rows. `j` cuts its columns
through the column chain, as today.

**Combining:**
- **`cbind`:** only the column map applies.
- **`rbind`:** SingleCellExperiment requires identical altExps (checked
  in R), so the row map is skipped for altExp roots, and their rows and
  columns don't change.
- **`combineCols()`:** the column map applies, and the main object's row
  map is skipped for altExp roots. The altExp's rows become the union
  (checked in R), and step 3 matches them by name, or by position with
  `use.names = FALSE`.
- **`combineRows()`:** the result has no altExps, so step 2 drops altExp
  roots with the warning.
- **Classes:** `.combineSubsets()` converts each internal object to the
  result's class, except for altExp roots and the subsets below them.

## Tests

Tests are written first and fail on their assertions against devel.
The fixture, built with a helper in the test file:
- 10 genes × 20 cells;
- an ADT `SingleCellExperiment` of 4 proteins, the last an isotype
  control;
- subsets `clean` (cells 1–12), whose ADT stores its own `clr`, and
  `tcells` (cells 1–6 of `clean`);
- a variant whose ADT is a `SummarizedExperiment`.

1. **Creating:**
   - from the main ADT and from `"clean/clr"`, with `rows` and `cols` as
     positions and names;
   - the default `parentAssay` (the first assay, stored explicitly);
   - dims, row and column names, and `subsetParent()`.
2. **Errors:** every new message in Interface, and `altExp` on a
   `SubsetSummarizedExperiment`.
3. **Assays seen:**
   - `assay()` of the subset and of `"clean.antibodies/counts"`, `/clr`
     equals the altExp view cut to its rows and columns;
   - `subsetAssayNames()`.
4. **Storing:**
   - an own assay;
   - a linked one needs `unlink`, and the message names "altExp 'ADT' of
     subset 'clean'";
   - a PCA in an SCE root;
   - `reducedDim<-` failing in the SE variant.
5. **Below a root:** a subset of `clean.antibodies` by rows and by
   columns, and its assays, `getSubsetExperiment()` and class.
6. **Data:** `rowData`/`colData` with parent data; `getSubsetExperiment()`
   class, dims and assays (primary first); a `setSubsetExperiment()`
   round trip.
7. **altExps:** `altExpNames()` of the subset doesn't list the main
   object's altExps, and an altExp the subset stores is seen below it.
8. **Summary:** `subsetSummary()` prints the Parent chain.
9. **Tree and spatial,** on a SpatialExperiment and a
   TreeSummarizedExperiment main object, each with an ADT
   SingleCellExperiment altExp and a subset made from it:
   - `subsetRowTree()`, `subsetRowTree<-` and `subsetRowLinks()` stop;
   - `subsetColTree()`, `subsetColLinks()`, `subsetSpatialCoords()`,
     `spatialCoords(subsetName =)` and `subsetImgData()` equal the main
     object's at the subset's cells;
   - with coordinates or a column tree stored on the container (`clean`),
     the ADT subset sees those instead;
   - `subsetSpatialCoords<-` and `subsetColTree<-` stop with the class
     message;
   - with an ADT that is a SpatialExperiment, `subsetSpatialCoords<-`
     stores coordinates on the ADT subset.
10. **`[`:** `es[1:5, ]` leaves the root's rows alone; `es[, 3:20]` cuts
    its cells and those below.
11. **Other rows:**
    - dropping the isotype control from the main ADT, and from `clean`'s
      ADT stored whole: the subset and those below keep the remaining
      features, with internal data cut to match;
    - reordered rows;
    - unnamed rows of the same number;
    - unmatched unnamed rows: an error, with the object unchanged.
12. **Rename:** `altExpNames(es) <- "protein"` updates the slot, and the
    subset still reads its assays.
13. **Removed:** `altExp(es, "ADT") <- NULL` and `altExps<-` without ADT
    drop the root and those below, with the warning; the case where
    `clean` stores its own whole ADT keeps them.
14. **Class change:** replacing ADT with a `SummarizedExperiment`
    converts the internal objects.
15. **Combining:** `cbind`, `rbind`, `combineCols()` (union rows,
    `use.names = FALSE`) and `combineRows()` (dropped).
16. **Old objects:** `tests/testthat/fixtures/es-1.23.54.rds`, saved from
    devel before the class change by
    `tests/testthat/fixtures/make-es-1.23.54.R`, an object with two
    chained subsets:
    - `assay()` of each subset and `createSubset()` work as loaded;
    - after `updateObject()`, `.hasSlot(subset, "altExpName")` is `TRUE`
      and the results are the same.
17. **Updated on purpose:** the S4 method snapshot from #118 gains the
    `updateObject` methods. The rest of the suite stays unchanged.

## Housekeeping

- **`DESCRIPTION`:** version 1.23.55. `Config/roxygen2/version` stays
  8.1.0.
- **`NEWS.md`, 1.23.55:** `createSubset(altExp =)`, what such a subset
  sees and stores, how altExp changes carry along, the new slot and
  `updateObject()` (#92).
- **Docs:**
  - `createSubset()`: `@param altExp` and a section on subsets made from
    an altExp, with an example on simulated data;
  - `ExperimentSubset-altExp` details: a sentence and a link;
  - the `updateObject` page.
- **`_pkgdown.yml`:** the `updateObject` page, then `make site-check`.
- **`AGENTS.md`:** the Object model section gets the `altExpName` slot,
  `.altExpName()` and `.copySubset()`.
- **ADR 0007** and the README index.
- **Checks:**
  - `make test`;
  - `make coverage` (must not drop from devel);
  - `make lint`, `make site-check`;
  - `make check-full`, `make bioccheck`;
  - the built vignette compared with devel's.
- **Review:** a spec review by a fresh subagent, then `/code-review`.
- **Follow-up PR:** the article section, and the last two steps of
  `dev/scripts/check-adt-workflow.R`.

## Found during implementation

- **The accessor is `.rootAltExp()`, not `.altExpName()`.** An internal
  `.altExpName(e, visible)` already resolves an altExp name or position,
  so the slot's accessor has another name. The slot is still
  `altExpName`.
- **Small helpers the design implied:**
  - `.altExpView(x, container, name)`: an altExp as the main object or a
    subset sees it, or `NULL`; `.rootObject()` and the carry step use it;
  - `.rootChain()`: the part of a chain whose data a subset sees, from
    the subset that starts the root down;
  - `.altExpStart()` and `.createAltExpSubset()`: the checks and the
    building for `createSubset(altExp =)`;
  - `.emptyCut()`: an altExp, or a parent's internal object, cut and
    emptied;
  - `.carryAltExps()`: runs `.carryAltExpParts()` then
    `.carryRowRoots()`, in every place that called the first;
  - `.checkTreeDim()`, `.checkCanStore()` and `.treeDims()` for trees and
    spatial data.
- **A subset of an altExp-rooted subset** gets its internal object, and
  so its row and column names, from its parent's internal object, cut
  and emptied. That keeps the class, as the spec says. Main-rooted
  subsets still take their names from the main object.
- **Link messages name each assay's own source.** `.linkedSources()` now
  takes `x`. For a subset made from `clean`'s view of ADT, `counts` is
  "linked to the main altExp 'ADT'" and `clr` "to altExp 'ADT' of subset
  'clean'", from the view's owners, not one wording for the whole root.
- **`updateObject()`** comes from S4Vectors, which re-exports it, so
  there's no new dependency. The `AssaySubset` method also runs
  `updateObject()` on the subset's stored object.
- **Combining:** `cbind()` and `rbind()` match altExp rows by position
  (their altExps must match), as `combineCols(use.names = FALSE)` does.
  The SpatialExperiment step that copies `sample_id` into each subset
  skips subsets made from an altExp, which store no copy.
- **`subsetSpatialData<-` is already defunct**, so it needs no class
  check.
- **`as(<SingleCellExperiment>, "SummarizedExperiment")` drops the
  rownames** (in SingleCellExperiment, not this package). Replacing ADT
  with such an object is then "other rows without names", an error. The
  class-change test sets the rownames again.
- **roxygen2 8.1.1 is installed.** `make docs` rewrote only the changed
  Rd files and bumped `Config/roxygen2/version`, which is set back to
  8.1.0, as AGENTS.md asks.
- **The full suite takes over 10 minutes** here, so it ran in the
  background.
- **Only the main object's altExps can be roots, now checked.** The spec
  left out an altExp stored only in a subset as a root, but nothing
  stopped it. `createSubset()` now stops with "'own' is stored only in
  subsets; a subset can be made only from an altExp of the main
  object." So a subset can't rename a root's altExp (renaming a linked
  altExp is already an error), and `.carryRowRoots()` follows only
  renames on the main object.
- **A dead branch removed.** `getSubsetExperiment()` converted its result
  when its class differed from the subset's internal object's, but the
  carry step keeps them the same.
- **Coverage:** the first full run gave 98.52% against devel's 98.69%
  (Codecov, 3c7b519). Tests were added for the paths it missed.
- **An existing bug, not from this branch:** storing an altExp in a
  subset (`altExp(x, e, subsetName = s) <-`) fails with "cannot subset by
  character when names are NULL" when another altExp the subset sees has
  an unnamed assay. To be logged as an issue.
