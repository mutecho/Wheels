# Current Status

## Task Snapshot

- scope: implement finite-source Coulomb fit support from
  `docs/plan/fit_finite_coul.md`
- current conclusion: `fit` now supports explicit Coulomb modes
  `none`, `gamow`, and `finite_source`; finite-source mode uses optional local
  CATS tables keyed by centrality/mT, seeded from the corresponding `phi_all`
  slice, and writes the selected mode plus finite-source radius into fit
  metadata
- primary evidence:
  - `[fit].coulomb_mode` parses `none|gamow|finite_source`; legacy
    `[fit].use_coulomb` remains accepted only when it maps unambiguously to
    `none` or `gamow`
  - `[fit].finite_source_mode` parses `fixed_1d|iterative_1d` and is rejected
    unless `coulomb_mode = "finite_source"`
  - CMake enables CATS-backed finite-source support when the local CATS/GSL
    libraries are available, and a no-CATS build fails finite-source workflows
    explicitly instead of silently falling back
  - finite-source fitting builds a one-dimensional Coulomb kernel from CATS
    using `k* = 0.5 * q_inv`, seeds each centrality/mT group from `phi_all`,
    and supports one-pass fixed or one-pass iterative source-radius updates
  - fit summaries preserve `usesCoulomb` and add `coulombMode`,
    `finiteSourceMode`, and `finiteSourceRadiusFm`; detailed and report ROOT
    outputs include `meta/CoulombKernelCatalog`
  - `workflow_smoke_test` covers `fixed_1d`, `iterative_1d`, and the no-CATS
    finite-source failure path; `coulomb_kernel_validation_test` covers the
    q-to-k mapping, table interpolation/clamping, no-FSI unity reference, and
    finite-source radius ordering
  - `scripts/cmake.sh` now clean-builds the selected build tree by default and
    verifies that `bin/exp_femto_3d` links CATS whenever the current build rule
    expects CATS, preventing a stale no-CATS binary from surviving in the
    shared source-tree `bin/` output directory
  - `2026-06-21` O2Physics ROOT executor configure/build/`ctest
    --output-on-failure` returned `PRIMARY_OK` in both CATS-enabled and
    no-CATS build matrices; all five registered tests passed in both matrices
  - `2026-06-22` O2Physics ROOT executor `scripts/cmake.sh` returned
    `PRIMARY_OK`, relinked the default CATS-enabled targets, `otool -L` showed
    `libCATS` and GSL on `bin/exp_femto_3d`, and `ctest --output-on-failure`
    passed all five registered tests

## Previous Snapshot

- scope: add a project-local run script for the OO 3D build/fit workflow
- current conclusion: `scripts/run_exp_femto_3d.sh` is now the default
  operator entry for running `build-cf` followed by `fit` with
  `config/oo_build_and_fit.toml`
- primary evidence:
  - default config resolves to `Exp_femto_3d/config/oo_build_and_fit.toml`
  - default stage is `all`, with `--stage build-cf` and `--stage fit`
    available for partial reruns
  - fit stages accept `--model full|diag` and `--input-cf-root <path>`
    overrides before forwarding to the public CLI
  - the script re-enters `O2Physics/latest-master-o2` through `alienv` when the
    ROOT runtime is not already active
  - `2026-06-10` `bash -n scripts/run_exp_femto_3d.sh` and
    `scripts/run_exp_femto_3d.sh --help` passed
  - the no-argument `alienv` re-entry path now avoids empty-array expansion
    under macOS Bash 3.2 `set -u`
  - the full OO real-data build/fit was not rerun because the default config
    writes production ROOT/TSV/report outputs

## Earlier Snapshot

- scope: add a standalone fit report ROOT output for `fit` results
- current conclusion: `fit` writes an independently configured report ROOT file
  in addition to the legacy `fit_root_name` and TSV summary outputs
- primary evidence:
  - `[output].fit_report_directory` and `[output].fit_report_root_name` parse
    from TOML; `fit_report_directory` defaults to `output_directory` for older
    configs
  - report ROOT includes `meta/FitCatalog` with the same fit-result fields as
    the TSV summary
  - report ROOT mirrors `summary/R2_vs_phi/...`
  - report ROOT adds
    `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`
  - report ROOT adds `eps_vs_mt/<cent>/epsf_vs_mt` and
    `eps_vs_mt/<cent>/epsf_vs_mt_canvas`
  - `fit` rejects report paths that resolve to the same file as the CF ROOT or
    detailed fit ROOT output
  - `2026-06-10` O2Physics ROOT executor `ctest --output-on-failure` passed all
    three registered tests and the workflow smoke checks the new report objects

## Historical Snapshot

- scope: sync `project-state/` with the current `3d_cf_from_exp` refactor work
  around phi-mapping persistence/override, explicit progress-mode control, and
  mixed-event phi-splitting control
- current conclusion: the active worktree has advanced beyond the earlier
  ROOT-runtime-diagnosis state; the current implementation now treats phi
  mapping as durable file metadata and lets fit follow or override it
- primary evidence:
  - `[build].progress` and `[fit].progress` now parse `true`, `false`, and
    `"auto"`
  - build writes `build_uses_symmetric_phi_range` into `meta/SliceCatalog`
  - build writes `split_mixed_event_by_phi` into `meta/SliceCatalog`
  - fit follows input CF metadata by default and can override it via
    `[fit].map_pair_phi_to_symmetric_range`
  - legacy `SliceCatalog` trees without the new branch are still readable
    through raw/display phi inference
  - legacy `SliceCatalog` trees without `split_mixed_event_by_phi` read it as
    `false`
  - `2026-04-19` non-sandboxed O2Physics `ctest --output-on-failure` passed
    all registered tests
  - `2026-05-18` O2Physics ROOT executor `ctest --output-on-failure` passed all
    three registered tests after adding split-ME coverage

## Verification Status

- verification_status: locally verified for config parsing, Coulomb kernel
  behavior, and ROOT-backed smoke coverage
- project_state_sync_status: written

Reason:

- `2026-06-21` O2Physics ROOT executor configure/build/`ctest
  --output-on-failure` returned `PRIMARY_OK` for the default CATS-enabled build
- the CATS-enabled configure step reported
  `CATS finite-source Coulomb support enabled` using the local CATS install
- the CATS-enabled `ctest` run passed 5/5 tests:
  `coulomb_kernel_validation_test`, `config_parse_validation_test`,
  `progress_render_test`, `slice_catalog_roundtrip_test`, and
  `workflow_smoke_test`
- `2026-06-21` O2Physics ROOT executor configure/build/`ctest
  --output-on-failure` also returned `PRIMARY_OK` with
  `-DEXP_FEMTO_3D_ENABLE_CATS=OFF`
- the no-CATS `ctest` run passed the same 5/5 registered tests, including the
  explicit finite-source unavailable failure path
- `2026-06-22` `scripts/cmake.sh` returned `PRIMARY_OK` through the O2Physics
  ROOT executor, clean-rebuilt the default build, and produced a CATS-linked
  `bin/exp_femto_3d`
- `2026-06-22` `ctest --test-dir Exp_femto_3d/build --output-on-failure`
  returned `PRIMARY_OK`; all five registered tests passed after the build
  helper change
- `2026-06-22` `otool -L bin/exp_femto_3d | grep -E 'CATS|gsl'` showed
  `libCATS.dylib`, `libgsl`, and `libgslcblas`
- `git diff --check` passed after the finite-source implementation and ledger
  updates
- full real-data physics regression on production OO/PbPb inputs has not yet
  been rerun; current coverage is toy ROOT workflow smoke plus direct kernel
  behavior checks
- `2026-06-10` O2Physics ROOT executor `cmake -S ... -B ...`, `cmake --build
  ... -j4`, and `ctest --output-on-failure` all returned `PRIMARY_OK`
- `ctest` passed the earlier four-test matrix after Eventgen-style progress
  rendering was added
- `2026-06-10` `bash -n scripts/run_exp_femto_3d.sh` and
  `scripts/run_exp_femto_3d.sh --help` passed for the OO run script entry
- the full `oo_build_and_fit.toml` real-data build/fit remains intentionally
  unrereun during local smoke work because it writes configured production
  outputs under `/Users/allenzhou/ALICE/alidata/femtoep_res/OO` and
  `Exp_femto_3d/res`
- authoritative local test execution for the phi-mapping worktree passed in a
  clean non-sandboxed O2Physics environment on `2026-04-19`
- O2Physics ROOT executor `ctest --output-on-failure` passed on `2026-05-18`
  after the mixed-event phi-splitting switch was added
- sandboxed runs that fail during `alienv` bootstrap with `/dev/fd/... Operation
  not permitted` remain non-authoritative environment noise

## Active Constraints

- ROOT-dependent validation must be run from a fully entered O2Physics
  environment
- finite-source Coulomb fits require CATS support at configure/build time;
  builds without CATS reject `coulomb_mode = "finite_source"` during workflow
  execution
- CATS kernel tables are one-dimensional in `k*` and seeded per centrality/mT
  group from the selected `phi_all` slice; physics closure still requires a
  real-data regression on a known-good dataset
- `scripts/cmake.sh` defaults to a clean-first build because all build trees
  currently write executables into the shared source-tree `bin/` directory;
  use `EXP_FEMTO_3D_CLEAN_FIRST=0` only for a deliberate fast incremental loop
- `scripts/run_exp_femto_3d.sh` defaults to `config/oo_build_and_fit.toml` and
  executes both `build-cf` and `fit`, so running it without `--stage` writes the
  configured OO outputs
- progress remains controlled by `[build].progress` and `[fit].progress`; in
  `"auto"` mode the ETA line is shown only when `stderr` is attached to a TTY
- sandboxed tool runs that fail during `alienv` bootstrap are not reliable
  evidence for code-level ROOT regressions
- `build.split_mixed_event_by_phi = false` preserves the historical
  phi-integrated mixed-event denominator per centrality/mT group; `true` is the
  opt-in mode for denominators that follow each SE phi slice

## Ledger Convention

- `Exp_femto_3d` adopts `project-state/` as its coordination ledger path
- this project does not use the hidden `.project-state/` path

## Active Worktree Highlights

- `fit` exposes `fit.coulomb_mode = "none"|"gamow"|"finite_source"` and keeps
  legacy `fit.use_coulomb` compatibility for unambiguous none/Gamow configs
- `fit.finite_source_mode = "fixed_1d"|"iterative_1d"` controls whether the
  final CATS kernel uses the Gamow-seeded `phi_all` radius directly or performs
  one finite-source seed refit before rebuilding the final kernel
- detailed and report ROOT outputs write `meta/CoulombKernelCatalog`; TSV and
  `meta/FitCatalog` expose Coulomb mode and finite-source radius metadata
- CMake auto-detects local CATS/GSL and can be forced off with
  `-DEXP_FEMTO_3D_ENABLE_CATS=OFF` for explicit no-CATS validation
- `scripts/cmake.sh` is the default local build helper; it resolves the project
  root from the script path, clean-builds by default, and checks the final
  operator binary for CATS linkage when the current link rule includes CATS
- `scripts/run_exp_femto_3d.sh` provides the project-local OO run entry with
  runtime re-entry, stage selection, and fit overrides
- build and fit expose explicit progress-mode control plus Eventgen-style
  activity and ETA rendering
- `SliceCatalog` carries file-level build phi mapping metadata for downstream
  consumers
- fit can reinterpret stored `raw_phi_*` slices into either raw `[0, pi]` or
  symmetric `[-pi/2, pi/2]` summary coordinates without rebuilding CF files
- build can either reuse one integrated-ME denominator per centrality/mT group
  or project ME with the current SE phi range via
  `build.split_mixed_event_by_phi`
- legacy catalogs remain readable through mapping-state inference
- fit writes a standalone report ROOT file controlled by
  `output.fit_report_directory` and `output.fit_report_root_name`

## Coordination Ledger State

- `project-state/` is the active adopted coordination ledger for
  `Exp_femto_3d`
- this sync records the finite-source Coulomb implementation, CATS/no-CATS
  validation matrix, the hardened `scripts/cmake.sh` build helper, and the
  remaining real-data physics-regression gap
