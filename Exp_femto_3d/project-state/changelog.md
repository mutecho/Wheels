# Changelog

## 2026-08-13

- Resolved the remote/local conflict as a semantic union of three independent
  contracts: qn-aware ME denominators, configurable Levy fit parameters, and
  build-side mT/phi rebin.
- Composed phi and qn selections in the ME projection helper so split-ME mode
  follows the final merged phi interval and, when enabled, the current qn
  range; `phi_all` and qn-all retain full normal-axis coverage.
- Kept both qn policy and rebin provenance in SliceCatalog/FitCatalog/TSV and
  preserved legacy catalog defaults.
- Reconciled the Wenya config, restored the stable OO default runner, expanded
  combined qn/phi/mT coverage, and synchronized README, formula docs, and the
  project ledger.
- O2Physics ROOT executor build and full CTest returned `PRIMARY_OK`; all six
  registered tests passed in 32.06 seconds.

## 2026-08-12

- Added `[build.rebin.mt]` and `[build.rebin.phi]` with explicit enable
  switches, factor/range modes, `[[bins.phi]]`, pre-projection sparse-axis
  grouping, rebin-aware paths, and catalog/TSV provenance.
- Added preflight validation for edge alignment, divisibility, phi mapping
  seams, duplicate paths, and safe output preservation.
- Added `build_cf_rebin_test` and expanded config/catalog/workflow coverage.

## 2026-07-01

- Added optional `[build].split_mixed_event_by_qn`, preserved qn-integrated ME
  as the default, and recorded the policy in SliceCatalog with legacy fallback.
- Added optional `[fit.parameters.<name>]` controls for all main Levy
  parameters, including fixed lambda/alpha support shared by chi2 and PML.
- Added strict config validation and ROOT-backed fit/qn regression coverage.

## 2026-06-30

- Added `[build].split_same_event_by_qn` and `[[bins.qn]]` support so
  qn-specific SAME slices can be written while preserving the legacy qn-all
  slices and keeping MIXED qn integrated.
- Updated `SliceCatalog`, `FitCatalog`, `CoulombKernelCatalog`, and TSV output
  metadata with `qn_index`, `qn_low`, `qn_high`, `qn_label`, and
  `is_qn_integrated`; legacy catalogs default to `qn_all`.
- Updated `config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml` for
  `/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb/wenya/3Dfemto_cent_mt_q2_phi_LHC23_merge.root`
  to write qn-all plus qn1/qn2/qn3 outputs in the standard
  `Exp_femto_3d` `SliceCatalog`/`slices` structure.
- Added config-parse, catalog roundtrip, and smoke coverage for qn metadata and
  qn split output paths.
- Verified the Wenya production build/fit through the O2Physics ROOT executor:
  `build-cf` stored `1040` slices, `fit` fitted `468/468` selected slices, ROOT
  inspection found qn_all/qn1/qn2/qn3 counts as expected, and the TSV has `469`
  lines.
- Added `docs/数学物理公式流程说明.md` as the current math/physics workflow
  reference for the refactored 3D femtoscopy pipeline.
- Documented the sparse-axis contract, SE/ME normalization, CF construction,
  phi coordinate mapping, `SliceCatalog`, diag/full Levy fit formulas,
  Gamow/finite-source Coulomb branches, PML objective, fit/report outputs,
  `R2_vs_phi`, and `epsf_vs_mt` semantics with code-location pointers.
- Linked the new formula workflow document from `README.md`.
- Synced `project-state/` to record that this was a docs-only update; no
  analysis code or runtime artifacts were changed.

## 2026-06-23

- Changed `scripts/cmake.sh` back to CMake's incremental build path by
  default; `EXP_FEMTO_3D_CLEAN_FIRST=1` now means an explicit full clean
  rebuild request.
- Added active ROOT runtime drift detection: when `root-config --prefix` points
  at a different `ROOTConfig.cmake` than cached `ROOT_DIR`, the helper clears
  ROOT CMake cache entries, reconfigures against the active ROOT, and refreshes
  link rules so targets relink as needed.
- Kept the CATS binary guard and added an on-demand CATS link-rule refresh for
  the stale no-CATS binary case, without forcing every ordinary build to
  recompile.
- Verified through the O2Physics ROOT executor that the first run refreshed
  `ROOT/v6-36-10-alice1-local7` to `ROOT/v6-36-10-alice1-local8`, the second
  run performed no compile/link work, `bin/exp_femto_3d` still links CATS/GSL,
  and `ctest` passed `5/5`.

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
