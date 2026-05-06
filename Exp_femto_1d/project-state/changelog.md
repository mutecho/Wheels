# Changelog

## 2026-05-06

- added `scripts/run_exp_femto_1d.sh` for config-driven `build-cf`/`fit` runs
- documented the run helper in README and synchronized `project-state/`
- checked run-helper shell syntax, help output, and non-ROOT stage dispatch
- fixed direct run-helper execution by auto-entering `O2Physics/latest-master-o2`
  when the caller shell lacks a complete ROOT runtime
- fixed the zero-argument helper path under Bash `set -u` by avoiding empty
  array expansion during runtime re-entry
- detached internal `THnSparse` projections immediately from ROOT directories
  to remove repeated `TROOT::Append` replacement warnings during `build-cf`
- changed the run-helper default stage to `build-cf`; full `build-cf` plus
  `fit` now requires explicit `--stage all`

## 2026-04-20

- created standalone `Exp_femto_1d` CMake project structure
- added public config, logging, workflow, and CATS model interfaces
- implemented structured `build-cf` output with `meta/SliceCatalog`
- added fit-side `FitCatalog`, per-slice output directories, region summaries, and TSV writing
- added example TOML configs and test entry points
- established `project-state/` for status, decisions, tests, issues, work items, and handoff
- passed non-sandboxed O2Physics configure/build/`ctest --output-on-failure` with 4/4 tests green
