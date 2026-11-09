# RNA and ADT workflows in the PBMC article (#93)

- **Status:** Draft, for review
- **Date:** 2026-11-09
- **Branch:** `feature/cite-seq-article`, from `devel`
- **Lands:** one PR to `devel` (1.23.49); devel only, no release port
- **ADR:** none. No package code, class or dependency changes; only a
  pkgdown article, its website packages, the pkgdown workflow's cache
  key, a dev script's header, NEWS and the version.

## Intent

Since #58, a subset sees the main object's altExps, cut to its cells, and
stores only the altExp assays it changes. `dev/scripts/check-adt-workflow.R`
showed on real data that RNA and ADT (CITE-seq) analyses can share one
object. Issue #93 asks for that as user-facing documentation.

Decided in chat (2026-11-09):

- **Move the existing article to the CITE-seq data** (issue option 2),
  rather than add a second article. One walkthrough covers RNA and ADT;
  there is no second copy of the RNA steps to maintain, and the site
  downloads one data set.
- **Store Cell Ranger's raw matrix as it comes:** all 6,794,880 barcodes
  and 33,570 features. The only change is that `splitAltExps()` moves the
  32 antibodies to `altExp("ADT")`. Storing everything is the point of
  the package. Measured on 2026-11-09: the count matrix is 213 MB
  without names (16.4 million non-zero entries; an all-zero barcode
  costs 4 bytes); with barcode names and colData, R uses about 1.8 GB;
  reading takes 11 s.
- **ADT steps in the article:**
  1. set-up (ADT in the main object's altExp, seen by every subset
     without a copy);
  2. decontPro and CLR on `clean`;
  3. the RNA-labelled T cells checked with protein;
  4. decontPro again on the T cells alone, with
     `unlink = "ADT/decontPro"`.
- **Left out:** renaming the main ADT and dropping its isotype controls
  (carry along). They come back once subsets can be made from an altExp
  (#92); the idea is recorded on #92.

Out of scope: the short Bioconductor vignette (it stays small and
offline), package code, the deprecated scuttle calls inside singleCellTK's
`runCellQC()` (singleCellTK's to fix), and decontPro's "Pareto k"
warnings (a decontX question; an issue for the decontX repo can be drafted
separately).

## The article

Same file, `vignettes/articles/pbmc-workflow.Rmd`, and URL, so the
vignette's two links stay valid. New title and description that name RNA
and protein; the run time and download size in the intro are updated
from the render. The disclaimer (a one-object analysis, not a definitive
PBMC annotation) stays and covers the protein checks too. The global
chunk options stay `message = FALSE, warning = FALSE`. The code follows
the user-facing rules: one step per line, a comment per step,
`library()` plus plain calls, `set.seed()` before anything random.

Sections, with what changes:

1. **Intro.** "4,000 PBMCs" becomes the 5k PBMC CITE-seq data. Says the
   antibody counts live in `altExp("ADT")` and that every subset sees
   them, cut to its cells, without copying them. Experiment and direct
   styles as now.
2. **The data.** `tenx-3.1.0-5k_pbmc_protein_v3` raw and filtered
   matrices from DropletTestFiles, read with `read10xCounts()`, gene
   symbols as row names (`uniquifyFeatureNames()`). `splitAltExps()` by
   `rowData()$Type` with `ref = "Gene Expression"`, the altExp named
   `ADT`, and `_TotalSeqB` stripped from the antibody names. Prose: a
   zero-count barcode is a whitelist entry never seen in sequencing, not
   an empty droplet; it is kept because it costs almost nothing.
3. **One object; finding cells; QC; ambient RNA; high-quality cells.**
   As now (`filtered`, `emptydrops`, `barcoderanks`, `runCellQC()`,
   decontamination, `clean`), on the new data, with the thresholds
   re-checked against it (the script used 15% mitochondrial counts,
   the article 10%; keep 10% unless the new data's distribution clearly
   needs more, and say so). One short new point: `filtered` and `clean`
   see the ADT of their cells and store none of it.
   `decontX()` comes from the decontX package instead of celda, after
   checking its arguments and results match.
4. **Clustering, annotation (SingleR), UMAP.** As now.
5. **New: removing ambient protein and normalising it** (direct style).
   decontPro on `clean`'s ADT counts, with the RNA clusters as cell
   types and the ADT of 20,000 droplets that aren't cells (and have any
   ADT counts) as the ambient profile, sampled after `set.seed()`. CLR as
   a short visible function. Both stored as assays of `clean`'s ADT with
   `altExp(es, "ADT", subsetName = "clean") <- adt`; `counts` stays
   linked. The object size before and after shows growth of about one
   matrix.
6. **New: checking the cell labels with protein.** A plot of CD3, CD4 and
   CD8a CLR values by SingleR label (with `fig.alt`), and the median CD3
   CLR of T cells against other cells.
7. **Re-analysing the T cells.** As now. The T cells'
   SingleCellExperiment carries ADT with `counts`, `decontPro` and
   `clr`, all linked.
8. **New: decontPro on the T cells alone.** decontPro fitted to the
   T cells only, with their own clusters as cell types. Storing it
   without `unlink` is an error naming the linked assay (shown in a chunk
   with `error = TRUE`); `unlink = "ADT/decontPro"` stores it, and
   `counts` and `clr` stay linked. This replaces the prose-only
   paragraph about `unlink = TRUE`.
9. **The whole analysis in one object.** `subsetSummary()`, and the size
   of the whole object next to that of the counts, with a sentence that
   the 6.8 million barcode names account for most of the difference.
10. **How this maps to MAMS.** The table gains a row: the `ADT` altExp is
    a MAMS feature-observation matrix (FOM) with `modality` "protein";
    `decontPro` and `clr` are `processing` values.
11. **Session information.** As now.

## Other files

- `DESCRIPTION`: `Config/Needs/website` replaces celda with decontX;
  `Version: 1.23.49`. Imports and Suggests unchanged.
- `.github/workflows/pkgdown.yaml`: the devel job's cache key
  `article-data-devel-v1` becomes `article-data-devel-v2`, since a cache
  whose key is found is not saved again and the new data would otherwise
  be downloaded on every run. The release job's key is unchanged: the
  release site keeps building the current article from `master`.
- `dev/scripts/check-adt-workflow.R`: only the header's last line
  changes, to point to the article and to #92 for the carry-along checks
  it still covers.
- `NEWS.md`: a 1.23.49 entry in the file's format: the workflow article
  now uses the 10x 5k PBMC CITE-seq data and adds an ADT workflow
  (decontPro, CLR, protein checks of the T cells, `unlink` of an altExp
  assay) on the same object (#93).
- `_pkgdown.yml`: unchanged (same article).

## Verification

No package code changes, so no new testthat tests.

- `make article FILTER=pbmc-workflow` renders with no error. In the
  rendered HTML: CD3 higher in RNA-labelled T cells than in other cells;
  CD4 above CD8a in CD4+ T cells and the reverse in CD8+ T cells; the
  object grows by about the decontPro matrix; no warnings or messages
  leak (the SeuratObject class-cache messages, decontPro's Pareto k
  warnings); total run time and peak memory noted.
- `make lint` and `make site-check` pass.
- `make test` and `make coverage` before hand-off; `make check-full` and
  `make bioccheck` before the PR.

If run time or memory at 6.8 million barcodes is unreasonable (for
example `runDropletQC()` or `object.size()` far slower than on the 4k
data), stop and bring it back with numbers before changing the design.
