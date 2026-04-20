# Changelog

## 2026-04-20

- created standalone `Exp_femto_1d` CMake project structure
- added public config, logging, workflow, and CATS model interfaces
- implemented structured `build-cf` output with `meta/SliceCatalog`
- added fit-side `FitCatalog`, per-slice output directories, region summaries, and TSV writing
- added example TOML configs and test entry points
- established `project-state/` for status, decisions, tests, issues, work items, and handoff
- passed non-sandboxed O2Physics configure/build/`ctest --output-on-failure` with 4/4 tests green
