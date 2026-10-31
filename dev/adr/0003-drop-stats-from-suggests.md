# 0003. Drop stats from Suggests

- **Status:** Proposed
- **Date:** 2026-10-31
- **Deciders:** Joshua D. Campbell

## Context

DESCRIPTION lists `stats` in Suggests. The review of the package metadata
(issue #33) flagged it: `stats` is a base package that ships with every R
installation and is attached in every session, including the ones
`R CMD check` uses for examples, tests and vignettes. PR #49 left it in
place because a DESCRIPTION dependency change needs an ADR.

The package calls `stats` functions only in places where `stats` is
already attached: `rnorm()` in the roxygen examples in `R/allGenerics.R`,
`rpois()` in the test fixture and in the vignette. Package code in `R/`
calls no `stats` function, and NAMESPACE imports nothing from it.

## Decision

We will remove `stats` from Suggests. Nothing else in DESCRIPTION changes.

## Consequences

- DESCRIPTION lists only packages that need installing, so
  `devtools::install_dev_deps()` and readers of Suggests aren't pointed at
  a base package.
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
