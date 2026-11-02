# Stop masking SingleCellExperiment's reducedDims<- (#26)

- **Status:** Draft for review
- **Date:** 2026-11-02
- **Branch:** `fix/reduceddims-generics`, from `devel`
- **Lands:** one PR to `devel` (1.23.30); devel only, no release port

## Intent

`R/allGenerics.R` defines a package generic `reducedDims<-` with the
signature `(x, subsetName, value)`. Because it is exported, it masks
SingleCellExperiment's `reducedDims<-` `(x, withDimnames = TRUE, ...,
value)` whenever ExperimentSubset is attached, so ordinary SCE code
breaks:

```r
library(ExperimentSubset)
reducedDims(sce, withDimnames = FALSE) <- l
#> Error: unused argument (withDimnames = FALSE)
```

This was reproduced on `devel` at 35a5daf for a plain SCE and for a
`SubsetSingleCellExperiment`. `SingleCellExperiment::reducedDims(...)`
works, so code in other packages that uses the namespace is not
affected.

SCE's generic already has `...`, so the package generic isn't needed:
`reducedDim<-`, `reducedDims` and `altExps<-` already take `subsetName`
through `...` on SCE's own generics. After this change, `reducedDims<-`
does too.

`reducedDimNames()` and `reducedDimNames<-` are different. In
SingleCellExperiment 1.34.0 their generics are `(x)` and `(x, value)`,
with no `...`, so `subsetName` has no way through. The package keeps its
own generics for them (decided in the #26 evaluation). They mask SCE's,
but they forward every call without `subsetName` to SCE unchanged, so no
SCE code breaks. This change only cleans up their methods.

Not a release port: no results are wrong, and the failing calls give an
error rather than a wrong answer.

## Changes

### `reducedDims<-` (code)

- Remove the `reducedDims<-` generic from `R/allGenerics.R` and the
  `"ANY"` method from `R/ExperimentSubset.R`.
- Add a shared internal `.reducedDims<-`(x, ..., value) in
  `R/ExperimentSubset.R`, the same as `.altExps<-`: without `subsetName`
  it calls `callNextMethod()`; with it, it checks the name with
  `.isSubset()` and sets the reducedDims of the subset's internal
  object, passing the other arguments (`withDimnames`, ...) through.
- Add `setReplaceMethod("reducedDims", ...)` with
  `function(x, withDimnames = TRUE, ..., value)` for
  `SubsetSingleCellExperiment`, `SubsetTreeSummarizedExperiment` and
  `SubsetSpatialExperiment`, each forwarding to `.reducedDims<-`. No
  method for the SE and RSE classes: their subsets hold a
  SummarizedExperiment, which has no reducedDims, today and after.
- NAMESPACE (through roxygen): `importMethodsFrom(SingleCellExperiment,
  "reducedDims<-")`; `export("reducedDims<-")` goes, and
  `exportMethods("reducedDims<-")` stays, which still exports SCE's
  generic with the new methods.
- Behavior change: `reducedDims(x, subsetName = ) <- value` on an object
  that is not an `ExperimentSubset` (for example a plain SCE) now goes
  to SCE's method, which ignores `subsetName`, as `altExps<-` already
  does. Before, it gave "subset does not exist". Calls on
  `ExperimentSubset` objects keep their current results.

### `reducedDimNames()` and `reducedDimNames<-` (code)

- Keep both generics in `R/allGenerics.R` unchanged, including the
  `(x, subsetName, value)` signature of `reducedDimNames<-`, so
  positional calls keep working.
- Replace the three identical copies of each `"ANY"` method (in
  `R/SubsetSingleCellExperiment.R`, `R/SubsetTreeSummarizedExperiment.R`
  and `R/SubsetSpatialExperiment.R`; the last file in Collate wins) with
  one of each in `R/ExperimentSubset.R`.
- The setter checks the name with `.isSubset()` instead of
  `x@subsets[[subsetName]]`. Its error for an unknown subset changes
  from "`<name>` does not exist in the subsets slot of the object." to
  the message every other method gives: "`<name>` subset does not
  exist. 'subsetName' must be the name of a subset."

### Docs

- `?"ExperimentSubset-reducedDim"` documents the new `reducedDims<-`
  methods, with an example of `reducedDims(es, subsetName = ) <- list`.
  Its line saying `reducedDims<-` has its own page goes.
- The generated `man/reducedDims-set.Rd` and
  `man/reducedDims-set-ANY-method.Rd` go (`make docs`). `?"reducedDims<-"`
  then opens SCE's page, as `?"altExps<-"` does.
- `_pkgdown.yml`: the reference section lists
  `ExperimentSubset-reducedDim` and `starts_with("reducedDim")`, which
  still match; run `make site-check`.
- Vignette: check its list of supported methods (line 323 and the
  numbered list) still reads correctly; no code there changes.

### Tests (`tests/testthat/test-ExperimentSubset.R`), written first

1. New: `reducedDims<-` takes `withDimnames` (#26). For a plain SCE and
   for `SubsetSingleCellExperiment`, `SubsetSpatialExperiment` and
   `SubsetTreeSummarizedExperiment` (the objects of the #77 test),
   `reducedDims(x, withDimnames = FALSE) <- list(PCA = pca)` runs and
   stores `pca`, on the main object and, for the `Subset*` classes, on
   a subset with `subsetName =`. Fails today with "unused argument".
2. New: the package's `reducedDims<-` is SCE's generic
   (`attr(ExperimentSubset::\`reducedDims<-\`, "package")` is
   `"SingleCellExperiment"`). Fails today.
3. New: for the three `Subset*` SCE classes, `reducedDims<-`,
   `reducedDimNames()` and `reducedDimNames<-` with an unknown
   `subsetName` give the `.isSubset()` error, and `reducedDimNames<-`
   with a known one renames the subset's reducedDims. Guards the moved
   methods for each class, not only the last in Collate.
4. Changed: the existing expectation at line 312 for
   `reducedDimNames(es, subsetName = "hvg123") <- ...` matches the new
   message ("hvg123 subset does not exist."). It still checks that the
   call fails; only the wording it expects changes.

`make test` and `make coverage` (must not drop) before hand-off;
`make check-full` and `make bioccheck` before a PR.

### NEWS and version

- Version 1.23.29 → 1.23.30.
- NEWS entry under 1.23.30: `reducedDims<-` no longer masks
  SingleCellExperiment's, so `withDimnames =` and other arguments work
  again, with and without `subsetName =`; `subsetName` on an object
  that isn't an `ExperimentSubset` is now ignored, as for `altExps<-`;
  `reducedDimNames<-` with an unknown subset gives the same error as the
  other methods (#26).

## Out of scope

- Dropping the `reducedDimNames` and `reducedDimNames<-` generics. That
  needs SCE's generics to gain `...`; a draft upstream issue can follow
  separately if wanted.
- Any other generic in `R/allGenerics.R`.
