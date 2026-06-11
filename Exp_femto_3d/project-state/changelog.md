# Changelog

## 2026-06-10

- Upgraded build/fit progress rendering to the Eventgen-style line with stage
  label, 50-column bar, percent, activity frame, and ETA.
- Added a one-second heartbeat for enabled progress output so long ROOT work
  keeps showing an active frame and refreshed ETA between completed slices.
- Added `progress_render_test` and linked the core target with
  `Threads::Threads`; O2Physics ROOT executor configure/build/`ctest
  --output-on-failure` passed all four registered tests.
- Added `scripts/run_exp_femto_3d.sh` as a project-local OO run wrapper that
  defaults to `config/oo_build_and_fit.toml`, runs `build-cf` then `fit`, and
  supports `--stage`, `--model`, `--input-cf-root`, and `--binary` overrides.
- Fixed the wrapper's no-argument `alienv` re-entry path so macOS Bash 3.2
  `set -u` does not fail on an empty preserved-argument array.
- Verified the new run wrapper with `bash -n` and `--help`; the full
  `oo_build_and_fit.toml` real-data workflow was not rerun because it writes
  configured production outputs.
- Added standalone fit report ROOT output controlled by
  `[output].fit_report_directory` and `[output].fit_report_root_name`.
- Kept older configs valid by defaulting `fit_report_directory` to
  `output_directory` and `fit_report_root_name` to `fit_report.root`.
- Wrote report `meta/FitCatalog` so the ROOT report carries the same
  fit-result fields as the TSV summary.
- Mirrored legacy `summary/R2_vs_phi/...` graphs into the report ROOT file.
- Added `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`
  canvases using the Eventgen-style alpha/radius panel layout.
- Added `eps_vs_mt/<cent>/epsf_vs_mt` graphs and canvases from the fitted
  `Rside2(phi)` second-harmonic response.
- Added a fit-time guard that rejects report paths resolving to the same ROOT
  file as the CF or detailed fit output.
- Updated example and working TOML configs, README output contract, config
  parsing tests, and workflow smoke checks; O2Physics ROOT executor
  `ctest --output-on-failure` passed all three registered tests.

## 2026-05-18

- Added `[build].split_mixed_event_by_phi` with default `false` to preserve the
  historical integrated mixed-event denominator per centrality/mT group.
- Added opt-in per-phi ME projection so each SE phi slice can use a matching
  mixed-event denominator while the phi-integrated slice still uses the full
  phi range.
- Persisted `split_mixed_event_by_phi` in `meta/SliceCatalog`; legacy catalogs
  without the branch read it as `false`.
- Added split-mode `skipped_zero_mixed_event_slices` accounting and printed it
  in the build-cf CLI summary.
- Extended config and catalog roundtrip tests; O2Physics ROOT executor
  `ctest --output-on-failure` passed all three registered tests.

## 2026-04-19

- Synced `project-state/` to the current `3d_cf_from_exp` refactor work rather
  than leaving it at the earlier ROOT-runtime bootstrap state.
- Recorded the new build/fit progress-mode controls exposed through TOML.
- Recorded that build now persists `build_uses_symmetric_phi_range` into
  `meta/SliceCatalog`.
- Recorded that fit can either follow the input CF phi metadata or override it
  from stored `raw_phi_*` coordinates without rebuilding the CF file.
- Recorded backward-compatible inference for legacy `SliceCatalog` trees that do
  not yet carry the new branch.
- Recorded the authoritative `2026-04-19` non-sandboxed O2Physics `ctest` pass
  and the matching sandbox skip signature for context.

## 2026-04-11

- Bootstrapped `project-state/` for `Exp_femto_3d`.
- Adopted `project-state/` as the durable coordination ledger path for
  `Exp_femto_3d`.
- Confirmed that prior ROOT/THnSparse failures during agent execution were
  caused by sandboxed `alienv` bootstrap failure, not by project logic.
- Completed the previously incomplete ROOT-dependent tests in a clean
  non-sandboxed O2Physics environment.
- Added `ROOT_RUNTIME_AGENT_NOTE.md` as a durable diagnostic reference for
  future agents.
- Synced the ledger again to make the non-hidden `project-state/` convention
  explicit and record `project_state_sync_status`.
- Added the required Chinese `project-state/guide.md` after confirming the
  current skill explicitly requires it for bootstrap.
