# Current Status

## Task Snapshot

- scope: support all `fit_selection` slices in process-parallel profile-only
  execution and provide a strict 10-worker operator configuration
- current conclusion:
  `profile_only + process + fit_selection` is supported. The coordinator
  materializes the full ordered selection, assigns one complete
  `(centrality,mT,qn)` group per child, and each child is forced to an internal
  exact `listed` scope so it cannot repeat the global selection
- primary evidence:
  - added `config/oo_build_and_fit_6bins_profile_strict_parallel.toml`, preserving
    the baseline strict physics/scans while selecting all configured fit slices,
    using profile-only legacy TMinuit with `workers=10`
  - child process startup requires the marker, exact slice list, and temporary
    chunk output together; common numerical-library thread pools are fixed to 1
  - chunk reuse/merge now requires both readable `ProfilePoints` and
    `AttemptPoints` for every expected slice/scan
  - `ProfileExecution` preserves the parent scope and records selected
    slice/group counts plus configured/effective workers
- verification:
  `verified`: `2026-09-05` O2Physics ROOT executor returned `PRIMARY_OK`; build
  and full CTest passed 7/7 in 41.56 s. Two-group fit-selection process/resume,
  serial/process numerical equivalence, changed-selection digest rejection,
  missing-tree rejection, and production-output isolation passed. ROOT inspector
  returned `STATUS: OK`. Real strict estimate-only reported 84 slices, 12 groups,
  153,468 maximum attempts, and 10/10 workers without creating output. No real
  strict profile scan was started

## Previous Task Snapshot

- scope: resolve the remote/local merge and combine qn-aware ME splitting,
  configurable Levy parameters, and build-side mT/phi rebin
- current conclusion:
  the resolved workflow preserves all three feature contracts. mT/phi rebin is
  performed on sparse selection axes before q projection; ME projections can
  independently follow the merged phi interval and current qn interval; fit
  parameter overrides remain shared by chi2 and PML paths
- primary evidence:
  - `BuildCfConfig` and TOML parsing retain `split_mixed_event_by_qn`, both
    per-axis rebin specifications, `[[bins.phi]]`, and all ten
    `[fit.parameters.*]` overrides
  - `RunBuildCf()` applies the selected mT range to SE/ME and composes phi and
    qn denominator policies for all four split-switch combinations
  - `SliceCatalog`, `FitCatalog`, and TSV output retain qn policy plus mT/phi
    rebin provenance; legacy catalogs default to qn-integrated ME, legacy mT
    ranges, and native phi
  - the Wenya config keeps qn-all plus qn1/qn2/qn3 SAME slices, phi-matched ME,
    qn-integrated ME, explicit native phi, range-based mT grouping, and fixed
    `alpha = 2`
  - the operator script again defaults to `config/oo_build_and_fit.toml`; NeNe
    and six-bin OO runs remain available through `--config`
- verification:
  `2026-08-13` O2Physics ROOT executor build and full CTest returned
  `PRIMARY_OK`; all six tests passed in 32.06 seconds, including ROOT-backed
  combined qn/phi/mT coverage, catalog compatibility, fit-parameter smoke, and
  the complete workflow smoke test

## Previous Snapshot

- scope: add and validate configurable build-side mT/phi rebin support with an
  explicit on/off switch and visible rebin metadata
- current conclusion:
  `Exp_femto_3d` now accepts `[build.rebin.mt]` and `[build.rebin.phi]` with
  `enabled = false|true`, supports either factor or explicit range modes, keeps
  phi native when the new switch is omitted, and records the final rebin mode in
  `SliceCatalog` / `FitCatalog` / TSV metadata
- primary evidence:
  - `build.rebin.phi` now allows overlapping explicit `[[bins.phi]]` ranges and
    still rejects exact duplicates
  - build and fit catalog entries carry `mt_rebin_enabled`, `mt_rebin_mode`,
    `phi_rebin_enabled`, and `phi_rebin_mode`
  - `build-cf` resolves all target selections before projection, merges native
    bins in factor mode, and keeps `phi_all` on the full native phi span
  - invalid axis edges, factor divisibility, symmetric-phi seam crossings,
    duplicate paths, and dynamic fit mT selections fail before existing CF/fit
    ROOT outputs are reset
  - `config_parse_validation_test` covers explicit disabled, factor, range,
    and invalid rebin contracts
  - `build_cf_rebin_test` covers factor rebin, metadata persistence, path
    uniqueness, CF=1 proportional SE/ME closure, finite stored errors/Sumw2,
    q-axis preservation, and axis/seam failure modes
  - CMake build ran through the O2Physics ROOT executor with `PRIMARY_OK`
  - `ctest --test-dir build --output-on-failure` reported `100% tests passed`
    for all six registered tests in 16.61 seconds with executor status
    `PRIMARY_OK`
- verification: locally verified through the O2Physics ROOT executor with
  passing config, ROOT-backed, and full `ctest` coverage

## Earlier Snapshot

- scope: add the required math/physics formula workflow document for
  `Exp_femto_3d`
- current conclusion: `docs/数学物理公式流程说明.md` now documents the current
  sparse-input, CF-building, Levy/Coulomb-fitting, and summary-output flow with
  formulas, abstract operations, and code-location pointers
- primary evidence:
  - the document is grounded in `README.md`, `include/exp_femto_3d/Types.h`,
    `src/Config.cpp`, `src/Workflow.cpp`, `CMakeLists.txt`, public TOML
    configs, and the existing ROOT-backed tests
  - it covers the 7-axis `THnSparseF` contract, SE/ME visible-range
    normalization, `CF3D = SE_norm / ME_norm`, phi-integrated versus per-phi
    ME denominator semantics, raw/display phi mapping, and `SliceCatalog`
  - it records the diag/full Levy model, optional \(q^2\) baseline, PML
    objective, Gamow mode, CATS-backed 1D finite-source Coulomb mode, and
    `fixed_1d`/`iterative_1d` kernel-preparation flow
  - it records output meanings for `FitCatalog`, `CoulombKernelCatalog`,
    `R2_vs_phi`, source-parameter report canvases, and `epsf_vs_mt`
  - `README.md` now links the new formula workflow document
- verification: docs-only update; verified by source/code-reference review and
  Markdown inspection, with no runtime analysis code changed

## Older Snapshot

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
  - `scripts/cmake.sh` now defaults to CMake's incremental build path, refreshes
    stale ROOT cache entries only when the active ROOT runtime changes, and
    verifies that `bin/exp_femto_3d` links CATS whenever the current build rule
    expects CATS
  - `2026-06-21` O2Physics ROOT executor configure/build/`ctest
    --output-on-failure` returned `PRIMARY_OK` in both CATS-enabled and
    no-CATS build matrices; all five registered tests passed in both matrices
  - `2026-06-22` O2Physics ROOT executor `scripts/cmake.sh` returned
    `PRIMARY_OK`, relinked the default CATS-enabled targets, `otool -L` showed
    `libCATS` and GSL on `bin/exp_femto_3d`, and `ctest --output-on-failure`
    passed all five registered tests

## Older Run-Script Snapshot

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

## Older Historical Snapshot

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

## Oldest Historical Snapshot

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

- verification_status: partially verified for the current profile-likelihood
  change; ROOT-independent scan semantics and static checks pass, while the
  ROOT executor build/CTest/schema/canvas closeout is blocked. Earlier rebin,
  catalog, Wenya, and Coulomb evidence below remains historical and valid for
  those earlier snapshots
- project_state_sync_status: written

Reason:

- `2026-09-03` standalone `ProfileLikelihoodDriverTest` compiled and passed,
  and `git diff --check` passed after the final profile fixes
- the required O2Physics executor sandbox route failed during environment entry;
  the explicitly authorized escalated call was then rejected by the approval
  backend with HTTP 403, so ROOT runtime claims were not inferred from source
  inspection
- `2026-08-12` O2Physics ROOT executor build plus
  `ctest --test-dir build --output-on-failure` returned `PRIMARY_OK`; all six
  registered tests passed in 16.61 seconds
- the synthetic rebin test verifies SE/ME integral conservation, proportional
  `CF3D = 1`, finite non-negative errors with Sumw2, unchanged q axes, final
  physical metadata ranges, dynamic mT fit selection, and output preservation
  on build/fit validation failures
- `2026-06-30` added Wenya PbPb LHC23 qn-all plus qn1/qn2/qn3 production
  support and verified it with O2Physics ROOT executor build, `ctest`,
  production `build-cf`, production `fit`, ROOT catalog/report inspection, and
  TSV line-count checks
- the production Wenya `build-cf` wrote
  `/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb/wenya/EP_dependence_CF_wenya_lhc23_merge_qn_split_plus_integrated.root`
  with `1040` `SliceCatalog` entries and no skipped slices
- the production Wenya `fit` wrote the detailed fit ROOT, TSV summary, and
  report ROOT for `468` selected slices with no missing objects
- `2026-06-30` docs-only review created `docs/数学物理公式流程说明.md`, linked it
  from `README.md`, and synced `project-state/`; no C++/TOML runtime behavior
  changed in this pass
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
- `2026-06-23` `scripts/cmake.sh` returned `PRIMARY_OK` through the O2Physics
  ROOT executor, refreshed stale `ROOT/v6-36-10-alice1-local7` CMake cache
  entries to active `ROOT/v6-36-10-alice1-local8`, and a second helper run
  performed no compile or link work
- `2026-06-23` `ctest --test-dir build --output-on-failure` returned
  `PRIMARY_OK`; all five registered tests passed after the incremental helper
  change
- `2026-06-23` `otool` checks showed `bin/exp_femto_3d` has ROOT
  `LC_RPATH` at `ROOT/v6-36-10-alice1-local8/lib` and still links
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
- `scripts/cmake.sh` defaults to CMake's incremental dependency checks; use
  `EXP_FEMTO_3D_CLEAN_FIRST=1` only for a deliberately full clean rebuild
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

- `docs/数学物理公式流程说明.md` is the source-of-truth formula workflow document
  for the current `build-cf` and `fit` computation chain; update it when sparse
  axes, CF normalization, fit formulas, Coulomb kernels, output semantics, or
  cited code ownership changes
- `config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml` is the dedicated
  Wenya PbPb LHC23 production config; it preserves qn-all output and adds
  qn1/qn2/qn3 SAME-side slices while keeping the standard
  `SliceCatalog`/`slices` structure
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
  root from the script path, refreshes stale ROOT cache/link rules only when
  the active ROOT runtime changes, and checks the final operator binary for
  CATS linkage when the current link rule includes CATS
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
- this sync records the profile-likelihood implementation and its incomplete
  ROOT-runtime closeout, together with the earlier Wenya production,
  finite-source Coulomb, CATS/no-CATS, build-helper, and real-data regression
  history

## 2026-09-03 Documentation Closeout

- docs_reviewed:
  `README.md`, `config/examples/exp_femto_3d.example.toml`,
  `docs/数学物理公式流程说明.md`, and all adopted `project-state/` ledger files
- docs_written:
  `README.md`, public example configuration, formula workflow document,
  `current-status.md`, `handoff.md`, `tests.md`, `changelog.md`, `guide.md`,
  `decisions.md`, `issues.md`, and `work-items.md`
- docs_intentionally_unchanged: production OO/PbPb TOML configurations and
  `ROOT_RUNTIME_AGENT_NOTE.md`; this pass does not enable profiles in production
- docs_stale_candidates: none identified
- docs_missing: none
- formula_workflow_doc: updated with profile/slice definitions, reference
  convention, minimum validity, frozen finite-source kernel, ROOT schema, and
  implementation/test pointers
- project_state_sync_status: written
- closeout_write_status: small_patch_and_docs
- project-state closeout: checked

## 2026-09-04 Acceleration Closeout

- docs_reviewed: README, public example TOML, formula workflow document, staged
  configs/runner, and all adopted project-state ledgers
- docs_written: README, public example, formula workflow document,
  current-status, handoff, tests, changelog, guide, decisions, issues, work-items
- docs_intentionally_unchanged: preserved strict OO likelihood configuration,
  existing user-modified `run_exp_femto_3d.sh`, and unrelated OO configs
- docs_stale_candidates: none after resolving the prior executor-blocked issue
- docs_missing: real-data benchmark evidence remains an open work item, not a
  documentation omission
- formula_workflow_doc: updated for process ownership, checkpoint/atomic merge,
  HESSE semantics, and ordered PML bin caching
- project_state_sync_status: written
- closeout_write_status: implementation_and_docs
- project-state closeout: checked
