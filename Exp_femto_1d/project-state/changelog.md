# Changelog

## 2026-05-27

- added `[build].cf_rebin_factor` with default `1`; values above one merge
  adjacent k* bins in the SE/ME operands used for `CF1D`
- kept `SE_raw1d` and `ME_raw1d` in the original selected k* binning and
  persisted `cf_rebin_factor` in `meta/SliceCatalog` with legacy catalogs
  defaulting to `1`
- extended config/catalog/workflow smoke tests; local build, guarded `ctest`,
  and O2Physics ROOT executor `ctest --output-on-failure` passed all 4
  registered tests
- guarded ROOT histogram error-buffer initialization so projections, clones,
  and rebinned CF operands no longer emit duplicate `TH1D::Sumw2` warnings
  during build-cf
- re-ran local build, O2Physics ROOT executor `ctest --output-on-failure`, and
  direct `workflow_smoke_test`; all passed and the direct smoke output had no
  `TH1D::Sumw2` warnings

## 2026-05-18

- added `[build].split_mixed_event_by_phi` with default `false` to preserve the
  historical integrated MinBias mixed-event denominator behavior
- added opt-in split mode where mixed-event denominators follow the same
  event-plane region as each same-event slice
- persisted `split_mixed_event_by_phi` in `meta/SliceCatalog` with legacy
  catalogs defaulting to `false`
- added `skipped_zero_mixed_event_slices` accounting for split-mode denominator
  failures and printed it in the build-cf CLI summary
- extended config and workflow smoke tests; local build and O2Physics ROOT
  executor `ctest --output-on-failure` passed all 4 registered tests

## 2026-05-08

- added build-cf `cent_slices/<cent_id>/<region_name>/CFByMtCanvas` output
  canvases that overlay all available mT-bin CF histograms for each centrality
  and event-plane region
- extended ROOT-backed tests to verify the new canvas schema in both reopen
  and shared-output build modes
- passed local build, guarded local `ctest`, and O2Physics
  `ctest --output-on-failure`
- changed `CFByMtCanvas` to line-only markers-off output by default and added
  `build.cf_by_mt_show_markers` as the opt-in marker switch

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
