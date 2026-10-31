# NNNN. Drop stats from Suggests

- **Status:** Proposed
- **Date:** 2026-10-31
- **Deciders:** (pending maintainer approval)

## Context

DESCRIPTION lists `stats` in Suggests. The review of the package metadata
(issue #33) flagged it, because `stats` is a base package. PR #49 left it
in place because a DESCRIPTION dependency change needs an ADR.

`stats` ships with every R installation and is one of R's default
packages, so it is attached unless `R_DEFAULT_PACKAGES` excludes it.
`R CMD check` attaches it when it runs examples, tests and vignettes.

The package calls `stats` functions only in those places:
- `rnorm()` in the roxygen examples in `R/allGenerics.R`;
- `rpois()` and `kmeans()` in `tests/testthat/test-ExperimentSubset.R`;
- `rpois()` in the vignette.

Package code in `R/` calls no `stats` function, and NAMESPACE imports
nothing from it.

## Decision

We will remove `stats` from Suggests. Nothing else in DESCRIPTION changes.

## Consequences

- DESCRIPTION lists only packages that need installing, so readers of
  Suggests aren't pointed at a base package.
- Examples, tests and the vignette are unaffected, since `stats` is
  attached in all of them.
- If package code in `R/` ever calls a `stats` function, it must be added
  to Imports and imported in NAMESPACE (`@importFrom stats ...`), or
  `R CMD check` gives a NOTE about undefined globals. Listing it in
  Suggests would not have covered that case anyway.

## Alternatives considered

- **Keep it in Suggests.** Harmless to `R CMD check`, but it implies an
  optional dependency that doesn't exist, and it was flagged in review.
- **Move it to Imports.** Nothing in `R/` uses it, so this would add an
  import with no use.
