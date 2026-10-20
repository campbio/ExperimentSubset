# Adopt r-bioc-dev-standards

## What and why

Set ExperimentSubset up like the other campbio packages that use
campbio/r-bioc-dev-standards (v1), following its ADOPTING.md, so every
Claude Code session loads the shared standards, the guardrail permissions
and the lint-on-edit hook, and the branches follow the devel / RELEASE_X_Y
/ stable-master model. Choices follow decontX where they apply.

## Changes

- Claude tooling, copied from the `v1` tag: `dev/hooks/load-standards.sh`,
  `dev/hooks/lint-changed.sh`, `Makefile` (no extra targets, empty
  `PEOPLE_ONLY`), `.worktreeinclude`, and `dev/claude-settings.json`
  allowing `git fetch campbio` (the shared remote here) and `git blame`.
- `.lintr` from the template with `indent = 2L`: the code is written with
  2 spaces (about 1,300 lines start with 2 spaces and 490 with 4).
- `AGENTS.md` (layout, object model, tests, setup) and `CLAUDE.md`;
  `dev/adr/` from the templates.
- `.Rbuildignore` and `.gitignore` entries from ADOPTING.md step 5.
- Stable-branch workflows in their own commit, so it can be cherry-picked
  onto RELEASE_3_23: `sync-stable.yaml` (`stable-branch: master`),
  `pr-base-devel.yaml`, `pull_request_template.md`.
- Replace the outdated `R-CMD-check.yaml` (it ran only on master, on
  ubuntu-20.04 with v1 actions) with decontX's layout: Linux checks in the
  Bioconductor devel and RELEASE_3_23 docker images on every push and PR
  to devel and RELEASE_*, and macOS and Windows on push only.

No change to package code or results. Every new file is excluded from
the build, so there is no version bump and no NEWS entry.

## Testing

- `make help` lists the standard targets.
- `make lint` runs. Record the lint count; fix nothing.
- `make test` passes; `make coverage` records a baseline.
- `make bioccheck` reports no new errors and nothing about `.claude`.
- After merge, a developer runs `make claude-setup` outside Claude Code,
  restarts `claude`, and checks that the standards load (ADOPTING.md
  step 8).

## After merge

Push RELEASE_3_23 from Bioconductor to GitHub, cherry-pick the
stable-workflows commit onto it, and push to trigger the first sync.
RELEASE_3_23's `.Rbuildignore` lacks `^\.github$` (it is only on devel),
so add it in a second commit there; otherwise a later release fix pushed
to Bioconductor would ship `.github/` in the tarball. Then,
in the GitHub UI, add a master ruleset (restrict deletions, block force
pushes) and check the sync run.
