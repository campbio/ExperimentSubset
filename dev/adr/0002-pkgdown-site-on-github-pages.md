# 0002. Publish a pkgdown site: release at the root, devel under /dev/

- **Status:** Accepted
- **Date:** 2026-10-25
- **Deciders:** Joshua D. Campbell (approved with the spec in
  `dev/plans/2026-10-24-pkgdown-site-and-workflow-article.md`)

## Context

The package's only long-form documentation is its Bioconductor vignette,
which must stay short and build quickly on the Bioconductor machines. A
full single-cell workflow (empty droplets, QC, decontamination, clustering
and annotation on PBMC data) shows how subsets fit together, but it
downloads data and needs heavy packages (singleCellTK, celda, SingleR,
celldex and others) that the package itself doesn't use. Adding them to
Suggests would make every Bioconductor build and every check install them.

Users read documentation for the release they have installed, while
development happens on `devel`, which can differ in its API for months.
Bioconductor's release policy allows documentation fixes on the release
branch but not new features.

`master` is a CI-maintained copy of the current release, updated by the
shared "Sync stable branch" workflow with a bot token. Pushes made with
that token don't start other workflows.

## Decision

We will publish a pkgdown site on GitHub Pages, deployed to the
`gh-pages` branch by `.github/workflows/pkgdown.yaml`:

- The release version is built from `master` at the site root. It runs on
  `workflow_run` after "Sync stable branch" succeeds, which fires even
  though the sync's push doesn't, and skips the build when `master` hasn't
  changed since the last deploy (the commit is kept in `release-sha.txt`
  on `gh-pages`). That check runs on a bare runner, so the release image
  is only pulled when there is something to build. A manual dispatch from
  `master` always rebuilds.
- The devel version is built from `devel` under `/dev/` on each push.
  Pull requests to `devel` that touch the site's sources build it without
  deploying. If the site has no home page yet, because there is no
  release site, the first devel deploy adds a one-line home page that
  sends readers to `/dev/`; the release site replaces it. The devel jobs
  never otherwise write the root.
- Each site is built in the matching Bioconductor container
  (`bioconductor_docker:devel` and `:RELEASE_X_Y`, bumped each release)
  with read-only permissions. Separate deploy jobs, which run no package
  code, get `contents: write`; each site has its own concurrency group, so
  a devel deploy can't cancel a pending release deploy. The deploy action
  is pinned to a commit SHA, and `GITHUB_TOKEN` is the only secret.
- Long articles live in `vignettes/articles/`, which `.Rbuildignore`
  excludes, so Bioconductor never builds them. Their packages go in
  `Config/Needs/website`, which only the site build installs.

## Consequences

- Readers get documentation that matches the release they installed,
  plus a preview of devel.
- Heavy article dependencies stay out of Suggests and out of the
  Bioconductor build.
- `workflow_run` only fires from the default branch, so the release job
  only works once the workflow file is on `master`. That needs a one-time,
  documentation-only pull request to the release branch (spec item 4), and
  each later release carries the file over.
- The release container tag must be bumped every release, next to the one
  in `R-CMD-check.yaml`.
- The articles can break when their packages change, without the package
  or its checks failing; only the site build shows it.
- After the first deploy, the maintainer enables GitHub Pages once, with
  source `gh-pages`.

## Alternatives considered

- **A devel-only site.** Simpler, one job, but readers on the release
  would see functions and arguments their version doesn't have.
- **Building the release site on a schedule.** Avoids `workflow_run`, but
  rebuilds when nothing changed, or lags behind a release by up to a day.
- **Article packages in Suggests.** Every Bioconductor build and check
  would install them, and the article would have to build on the
  Bioconductor machines.
