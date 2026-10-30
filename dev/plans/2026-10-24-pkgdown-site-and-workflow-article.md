# pkgdown site, subset-experiment helpers and an end-to-end workflow article

- **Status:** Draft for review
- **Date:** 2026-10-24
- **Branch:** this spec is on `feature/pkgdown-site-workflow`; each PR in
  Delivery gets its own branch from `devel`
- **Lands:** after the Bioconductor 3.24 release, as several PRs to `devel`
  (see Delivery), plus one site-and-docs-only update to `RELEASE_3_24`

## Intent

Show users, on a public website, how ExperimentSubset keeps a whole
single-cell analysis in one object: raw droplets, the filtered cells, QC,
clustering, annotation and a re-analysed cell type all live as a chain of
subsets of the same object. The audience is analysts deciding whether the
package fits their workflow. Success means a reader can follow every step
of a real 10x analysis, see each stage as a subset, and reuse the pattern
with tools that only know about SingleCellExperiment.

Getting there needs two package changes first, so the article reads cleanly:
assay names that can repeat across subsets (#2), and helpers that hand a
subset to other tools and store the results back.

## Scope

In this spec, in build order:

0. **#2, linked and subset-specific assays.**
1. **`getSubsetExperiment()` / `setSubsetExperiment()`.**
2. **pkgdown site and GitHub Pages deployment** (release at the root, devel
   at `/dev/`), with ADR 0002.
3. **The workflow article**, `vignettes/articles/pbmc-workflow.Rmd`, with
   MAMS-aligned names and a MAMS section, and links to the site from the
   short vignette.
4. **A one-time site-and-docs update to `RELEASE_3_24`**, so the release
   site at the root works before the next release.

Out of scope: #4 (`[` subsetting), #26 (masked generics), #50 (inherited
rowData/colData), the #36 refactor beyond the shared name resolver, and the
lint backlog (#37).

## 0. Linked and subset-specific assays (#2)

Revised with the maintainer on 2026-10-25: each subset gets its own assay
namespace, made of assays linked from its parent and assays of its own.

**Today.** A subset exposes one assay, its slice of `parentAssay`. Assays
stored in subsets share one flat namespace: `assay(x, name)` looks at main
assays, then subset names, then stored assays, and the last match wins
silently. Nothing stops two related subsets from storing the same name with
different meanings.

**Terms.**

- *Container*: the main object or a subset.
- A subset's *parent container* is where its `parentAssay` lives: the main
  object for a main assay, subset `P` for `parentAssay = "P"` or
  `"P/assay"`.
- A subset's *linked assays* are all assays visible in its parent container,
  sliced to the subset's rows and columns. They are views computed on
  access by walking the parent chain (as `.rootIndices()` does), not copies.
  For a subset of the main object, these are the main assays.
- A subset's *own assays* are the ones stored in it with `setSubsetAssay()`.
- A subset's *visible assays* are its own assays plus its linked assays not
  hidden by an own assay of the same name. Its *primary* assay, what
  `assay(x, "subset")` returns, is its visible assay named after what its
  `parentAssay` refers to (followed through parent subsets). Usually that is
  its slice of `parentAssay`, as today; if the subset, or a subset on the
  way, stored its own assay of that name with `unlink = TRUE`, the primary
  assay is that one, consistent with "hides the linked one for that subset
  and every subset made from it".

Example: `logcounts` computed on the main object is linked into `filtered`,
`clean` and `tcells` under the same name, because log-normalisation is
per-cell. `scaled` (z-scores) stored in `filtered` is linked into `clean`
and `tcells`; a `scaled` computed on the T cells is a different assay and
must be stored deliberately, as below.

**Reading.**

- `assay(x, "s")` returns subset `s`'s primary assay (unchanged).
- `assay(x, "s", subsetAssayName = "a")`,
  `getSubsetAssay(x, "s", subsetAssayName = "a")` and the qualified name
  `"s/a"` return assay `a` as visible in subset `s`: its own if it has one,
  otherwise the linked view. `subsetAssayName` mirrors `setSubsetAssay()` and
  `assay<-`. (`getSubsetAssay()`'s second argument is already named
  `subsetName`, hence `subsetAssayName`.)
- A **bare name** means a main assay, or (as today) a subset name. A bare
  name that is only an assay stored in a subset is an error that names the
  subsets storing it and shows the qualified form. No deprecation cycle:
  no Bioconductor or CRAN package depends on, imports or suggests
  ExperimentSubset (checked 2026-10-25), and the maintainer chose to make
  the change directly.
- `createSubset(parentAssay = )` accepts a main assay, a subset name, or
  `"subset/assay"`. A child records its parent qualified
  (`"filtered/decontaminated"`), so later assays in other subsets can't
  change what it points to. Existing objects whose subsets record a bare
  stored-assay parent still resolve it when it is stored in exactly one
  subset, and `createSubset()` and `setSubsetAssay()` first rewrite such
  parents to `"subset/assay"`, so a later subset or assay with that name
  can't change or break them.
- Because stored assays are named `"subset/assay"`, a new subset may share
  its name with an assay stored in a subset; it may not share the name of a
  main assay.

**Writing.** `setSubsetAssay(x, subsetName, inputMatrix, subsetAssayName,
unlink = FALSE)` and `assay(x, "s", subsetAssayName = "a", unlink = ) <- m`:

- A name that is new to the subset (not own, not linked) is stored.
- A name the subset already owns is replaced, as today.
- A name that is linked is an error unless `unlink = TRUE`:
  `'scaled' in subset 'tcells' is linked to 'filtered/scaled'. Set
  unlink = TRUE to store a separate 'scaled' computed on the cells of
  'tcells'.` With `unlink = TRUE` the subset stores its own assay, which
  hides the linked one for that subset and every subset made from it. The
  parent's assay is unchanged.
- Sibling subsets don't see each other's assays, so they can each store the
  same name (decontX writing `decontXcounts` into several subsets works).
- `assay(x, "s/a") <- m` stores into subset `s` with the same rules.
- New subset names and stored assay names may not contain `/`.

**Listing.** `subsetAssayNames(x)` returns the subset names followed by
every own assay as `"subset/assay"`, so every listed name can be passed
back to `assay()`; `show()` prints the same. A new optional argument,
`subsetAssayNames(x, subsetName = "s")`, lists the names visible in subset
`s` (own and linked). This changes `subsetAssayNames()`'s output for objects
with stored assays; NEWS says so.

**Implementation.** One resolver maps a reference to a main assay, a subset,
or (container, assay); one function lists a subset's visible assays with
where each comes from; one function computes a linked view by composing the
subset's indices up to the container that owns the assay. `.assay()`,
`.parentSubsetName()`, `.owningSubset()`, `.subsetParent()`,
`.toRootIndices()`, `.rows()`/`.cols()` and `createSubset`'s parent check
use them. No new slots and no class changes.

**Behaviour changes.** Storing a name that is linked (for example
`logcounts` in a subset when the main object has `logcounts`) now needs
`unlink = TRUE`. Bare stored-assay names are an error. Listing output is
qualified. NEWS covers all three.

**Tests.** Linked views through one and two levels, matching a slice of the
source; own assays hiding linked ones only with `unlink = TRUE`, and only
for that subset and its descendants; siblings storing the same name;
qualified and `subsetAssayName =` reads and writes; the error for a bare
stored name; a child's parent
recorded qualified and unaffected by later assays elsewhere; `/` rejected;
the new `subsetAssayNames()` output; existing objects with bare parents.

## 1. `getSubsetExperiment()` / `setSubsetExperiment()`

New exported generics with methods for all five `Subset*` classes (thin
wrappers in `R/Subset<Class>.R`, shared code in `R/ExperimentSubset.R`).

**`getSubsetExperiment(x, subsetName, parentData = TRUE)`** returns a
standalone object of the main object's class (SummarizedExperiment,
RangedSummarizedExperiment, SingleCellExperiment, TreeSummarizedExperiment
or SpatialExperiment) for the subset:

- Assays: every assay visible in the subset (section 0): linked views under
  their source names (so `clean`'s primary assay is `decontaminated`, not
  `counts`) and its own assays.
- rowData and colData: the subset's own columns, plus with
  `parentData = TRUE` the main object's columns tagged as parent data (the
  `fromParentData` mark from #46).
- reducedDims, altExps and metadata: the subset's.
- Row and column names: the subset's own.

**`setSubsetExperiment(x, subsetName, value, unlink = FALSE)`** stores
what a tool produced:

- `value` must have the subset's dimensions and the same row and column
  names, in the same order; otherwise it is an error. Tools that drop cells
  or genes are followed by `createSubset()` instead.
- Each assay in `value` is compared with what the subset sees. Linked
  assays that are unchanged are skipped, so nothing is copied. New names are
  stored and own assays replaced, through `setSubsetAssay()`. A linked assay
  that the tool changed (for example a tool rewriting `scaled`) is an error
  unless `unlink = TRUE`, which stores it as the subset's own.
- rowData and colData replace the subset's own, without the tagged parent
  columns. reducedDims, altExps and metadata replace the subset's.
- Replacement, not merge: what was removed from `value` is removed from the
  subset. `value` normally comes from `getSubsetExperiment()`, so this is
  the expected behaviour; the help page says so.

**Tests.** A round trip for each of the five classes; a new assay,
reducedDim and colData column stored; unchanged linked assays not copied; a
changed linked assay erroring without `unlink = TRUE` and stored with it;
wrong dimensions and reordered names each error; parent columns not
stored; a subset made from another subset.

## 2. pkgdown site and GitHub Pages

**Site.** `_pkgdown.yml` with Bootstrap 5,
`url: https://campbio.github.io/ExperimentSubset/`,
`development: mode: auto` (overridden by `PKGDOWN_DEV_MODE` in CI), and a
reference index grouped as: constructor and classes; creating and managing
subsets; getting and setting subset data; class-specific accessors (links,
spatial); methods that take `subsetName`. The site URL goes first in
DESCRIPTION's `URL`. `make site-check` must pass. If roxygen's skipped method
topics (#15) block it, #15 is fixed in this work; otherwise it is left.

Articles: the existing vignette is "Get started"; the workflow article is
listed under "Workflows".

**Dependencies.** The article's packages go in `Config/Needs/website`, not
Suggests: singleCellTK, celda, DropletUtils, DropletTestFiles, scrapper,
SingleR, celldex, scater. `.Rbuildignore` gains `^_pkgdown\.yml$`, `^docs$`,
`^pkgdown$` and `^vignettes/articles$`; `.gitignore` gains `docs/`. The
Bioconductor build never installs these packages or sees the article.

**Deployment, `.github/workflows/pkgdown.yaml`.** Two jobs deploy to the
`gh-pages` branch:

- `devel-site`: on push to `devel`, on `workflow_dispatch`, and build-only on
  pull requests to `devel` that touch `R/`, `man/`, `vignettes/`,
  `_pkgdown.yml` or DESCRIPTION. Runs in
  `bioconductor/bioconductor_docker:devel`, installs the package and
  `Config/Needs/website`, builds with `PKGDOWN_DEV_MODE=devel` into `/dev/`,
  and deploys with `target-folder: dev`. Until `master` has `_pkgdown.yml`
  (below), it also builds the root, so the site isn't empty.
- `release-site`: on `workflow_run` when "Sync stable branch" completes
  successfully, and on `workflow_dispatch`. The sync updates `master` with
  the bot token, and pushes made with it don't start workflows, but
  `workflow_run` does fire. It checks out `master`, skips if `master` has no
  `_pkgdown.yml` or if its SHA equals the last deployed release SHA (kept in
  a file on `gh-pages`), and otherwise builds with
  `PKGDOWN_DEV_MODE=release` into the root and deploys with
  `clean-exclude: dev`. Runs in `bioconductor/bioconductor_docker:RELEASE_3_24`;
  the tag is bumped each release, as in `R-CMD-check.yaml`.
- Security, following the review of the release diff
  (campbio/r-bioc-dev-standards#13): `permissions: contents: read` at the top
  level, `contents: write` only on the two deploying jobs, the deploy action
  pinned to a commit SHA, and no secrets other than `GITHUB_TOKEN`.
- `workflow_run` workflows only run from the default branch (`master`), so
  the release job starts working once the workflow file reaches `master`,
  through item 4 below.

**Manual step.** After the first deploy, the maintainer enables GitHub Pages
in the repository settings with source `gh-pages`.

**ADR 0002** records: adding a pkgdown site with GitHub Pages; release/devel
split and the `workflow_run` trigger; `Config/Needs/website` for
article-only dependencies; and alternatives considered (devel-only site;
building the release site on a schedule).

**AGENTS.md** replaces "no pkgdown site" with the site's location, the
article directory, and the `make site-check` / `make article` targets.

## 3. The workflow article and the short vignette

**Data.** 10x PBMC 4k (Cell Ranger 2.1.0) from DropletTestFiles: the raw
matrix (about 737k barcodes) and the filtered matrix, read with
`DropletUtils::read10xCounts()`.

**Object and names.** One object `es`, with names following MAMS
(<https://github.com/single-cell-mams/mams>) suggested categories:

```
counts                        main assay, all droplets (MAMS processing counts; obs/feature_subset full)
├── background                empty droplets used as decontX's ambient profile
└── filtered                  Cell Ranger cells (obs_subset filtered)
    │  decontaminated         decontX output (processing decontaminated)
    └── clean                 cells passing QC (obs_subset clean; parent filtered/decontaminated)
        │  lognormalized      own (processing lognormalized)
        │  scaled             own (processing scaled)
        │  reduction.pca, embedding.umap
        ├── clean.variable    top variable genes (feature_subset variable)
        └── tcells            T cells (parent clean)
            │  lognormalized  linked from clean: per-cell, so not recomputed
            │  scaled         own, stored with unlink = TRUE: z-scores of T cells only
            │  reduction.pca, embedding.umap
            └── tcells.variable
```

**Steps.** Each step says which of two styles it uses: *experiment*
(`getSubsetExperiment()`, run a tool, `setSubsetExperiment()`) or *direct*
(store a matrix with `setSubsetAssay()`, or set `colData`/`rowData`/a
reducedDim on the subset).

1. Import raw and filtered; `ExperimentSubset()` on raw; `createSubset(es,
   "filtered", cols = <filtered barcodes>)`.
2. Droplet QC on the main object with singleCellTK `runDropletQC()`
   (emptyDrops, barcode ranks); compare its calls with Cell Ranger's.
3. Cell QC on `filtered`, experiment style: singleCellTK `runCellQC()`
   (QC metrics, mitochondrial percentage, scDblFinder).
4. decontX, direct style: `celda::decontX()` on `assay(es, "filtered")` with
   the `background` subset; `setSubsetAssay(..., "decontaminated")`;
   contamination into the subset's colData.
5. Filter on detected genes, mitochondrial percentage, doublets and
   contamination; `createSubset(es, "clean", cols = keep,
   parentAssay = "filtered/decontaminated")`.
6. Clustering, direct style, with scrapper: normalise (`lognormalized`),
   choose variable genes (`clean.variable`), scale them (`scaled`), PCA
   (`reduction.pca`), SNN graph and clusters (colData), UMAP
   (`embedding.umap`).
7. Annotation, direct style: SingleR with celldex's Monaco reference, main
   labels into colData; marker genes flagged in rowData.
8. T cells, experiment style: `createSubset(es, "tcells", ...)` with parent
   `clean`. `lognormalized` is linked from `clean` and reused, since
   log-normalisation is per-cell. Variable genes (`tcells.variable`),
   `scaled` (z-scores across T cells only, so stored with `unlink = TRUE`,
   with a short explanation of linked vs own assays), PCA, clustering and
   UMAP on the T cells; SingleR fine labels.
9. `subsetSummary(es)`; UMAPs of `clean` and `tcells` with scater, from
   `getSubsetExperiment()`; a note on memory compared with separate objects;
   `sessionInfo()`.

**MAMS section.** "How this maps to MAMS": a table relating subsets to
`obs_subset`/`feature_subset`, assay and reducedDim names to `processing`,
`parentAssay` to `parent_id`, and subsets of subsets to
`parent_relationship`, with a link to the specification.

**Style.** Written for R novices per the standards: one operation per line,
a comment per step, `set.seed()` before anything random, a `library()` call
per package and no `pkg::` prefixes in shown code. The article runs in full
when the site builds; a failing step fails the build.

**Short vignette.** `vignettes/ExperimentSubset.Rmd` gets a "Learn more"
note near the top and at the end, with absolute links to the site and to the
workflow article, so they work from the Bioconductor landing page.

**Checks while implementing.** singleCellTK's exact arguments; emptyDrops'
runtime on about 737k droplets (raise `lower` if needed rather than skip the
step); whether `getSubsetExperiment()` output is accepted as-is by
`runCellQC()` and scrapper.

## 4. One-time update to `RELEASE_3_24`

After 3.24 is released and the site work is merged to `devel`, a
release-branch PR (base `RELEASE_3_24`) carries only site and documentation
files: `_pkgdown.yml` (without the workflow article in its article list),
`.github/workflows/pkgdown.yaml`, the `.Rbuildignore`/`.gitignore` entries,
DESCRIPTION's `URL`, and the short vignette's "Learn more" links. No R code.
The release version's z is bumped; the sync workflow updates `master`; the
maintainer pushes `RELEASE_3_24` to Bioconductor. Bioconductor's release
policy allows documentation changes but not new features, which is why #2,
the helpers and the article stay on `devel` until the next release. From
then on the release site at the root documents 1.24, and `/dev/` documents
devel with the article.

## Delivery

Separate PRs to `devel`, in order, each on its own branch from `devel` and
with tests, NEWS and a z bump:

1. `feature/subset-assay-scoping`: #2, linked and subset-specific assays.
2. `feature/subset-experiment`: `getSubsetExperiment()` /
   `setSubsetExperiment()`.
3. `feature/pkgdown-site`: site, deployment workflow, ADR 0002, AGENTS.md,
   `Config/Needs/website`.
4. `feature/workflow-article`: the article and the short vignette's links.

This spec is committed on `feature/pkgdown-site-workflow` and goes in with
PR 1.

Then the release-branch PR in item 4. Each PR follows the usual review: a
spec check, `/code-review`, `make test`, `make coverage`, `make check-full`,
`make bioccheck`, and for 3 and 4 `make site-check` and
`make article FILTER=pbmc-workflow`. `make article` needs the package
installed locally; if it isn't, the developer is asked to install it.

## Risks

- The article's runtime on CI (emptyDrops on about 737k droplets, SingleR,
  UMAP twice). If the site build is too slow, cache ExperimentHub downloads
  between runs, then reduce emptyDrops work.
- singleCellTK and celldex are large dependencies of the site build only;
  version drift there can break the site without affecting the package.
- `workflow_run` only fires from the default branch, so the release job does
  nothing until item 4 lands on `master`.
