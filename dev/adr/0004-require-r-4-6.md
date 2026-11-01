# 0004. Require R 4.6.0

- **Status:** Accepted
- **Date:** 2026-11-01
- **Deciders:** Joshua D. Campbell

## Context

DESCRIPTION has `Depends: R (>= 4.0.0)`. BiocCheck 1.48.1 gives the
NOTE "Update R version dependency from 4.0.0 to 4.6.0" (issue #52).
Bioconductor 3.24, the current devel, is built and checked with R 4.6
(bioconductor.org/config.yaml), as was 3.23. Bioconductor only installs a
release on the R version it was built for, so users on older R don't get
this version of the package from Bioconductor anyway.

The declared floor is also not tested. No CI job or Bioconductor builder
runs the package on R 4.0 to 4.5, and issue #52 notes that the code
likely needs at least R 4.2, since it relies on `if()` failing on a
condition of length greater than one, which became an error in R 4.2.

## Decision

We will change the dependency to `R (>= 4.6.0)`. Nothing else in
DESCRIPTION's dependencies changes.

## Consequences

- The declared R version is the one the package is built and checked
  with, and BiocCheck no longer gives the NOTE.
- Installing this version from GitHub on R older than 4.6.0 fails at
  install time with a clear message, instead of running on an untested R.
- Each time Bioconductor devel moves to a new R version, BiocCheck will
  suggest raising it again; that is a new decision, recorded in a later
  ADR that supersedes this one.

## Alternatives considered

- **Keep R 4.0.0.** It claims support for R versions nothing tests, and
  BiocCheck keeps giving the NOTE.
- **Raise it to R 4.2.0, the oldest R the code likely needs.** It
  would also be untested, since Bioconductor does not build on R 4.2 to
  4.5, and BiocCheck would still give the NOTE.
