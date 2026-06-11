# Current Status

## Task Snapshot

- scope: upgrade build/fit progress rendering to match the Eventgen-style ETA
  display
- current conclusion: `build-cf` and `fit` progress now show stage label,
  percent, an activity frame, and ETA, with a one-second heartbeat while ROOT
  work is still running
- primary evidence:
  - `FormatProgressLine` renders the label, 50-column bar, integer percent,
    activity frame, and ETA from a value-only snapshot
  - `ProgressReporter` starts a heartbeat when progress is enabled and stops it
    before finish/abort cleanup
  - `Logger` progress rendering is mutex-protected so warning/info output closes
    the progress line before writing log messages
  - `progress_render_test` covers half-complete ETA, not-started unknown ETA,
    and over-complete clamping to 100%
  - `2026-06-10` O2Physics ROOT executor configure/build/`ctest
    --output-on-failure` returned `PRIMARY_OK`; all four registered tests
    passed

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

- verification_status: locally verified for progress formatting, config parsing,
  and ROOT-backed smoke coverage
- project_state_sync_status: written

Reason:

- `2026-06-10` O2Physics ROOT executor `cmake -S ... -B ...`, `cmake --build
  ... -j4`, and `ctest --output-on-failure` all returned `PRIMARY_OK`
- `ctest` now runs four tests: `config_parse_validation_test`,
  `progress_render_test`, `slice_catalog_roundtrip_test`, and
  `workflow_smoke_test`; all passed
- `2026-06-10` `bash -n scripts/run_exp_femto_3d.sh` and
  `scripts/run_exp_femto_3d.sh --help` passed for the new script entry
- `2026-06-10` `bash -u` branch checks passed for both empty and non-empty
  preserved-argument cases after the no-argument re-entry fix
- the full `oo_build_and_fit.toml` real-data build/fit was intentionally not
  rerun during script creation because it would write configured production
  outputs under `/Users/allenzhou/ALICE/alidata/femtoep_res/OO` and
  `Exp_femto_3d/res`
- O2Physics ROOT executor `ctest --output-on-failure` passed on `2026-06-10`
  after adding the standalone fit report ROOT file and report-object smoke
  checks
- authoritative local test execution for the current worktree passed in a clean
  non-sandboxed O2Physics environment on `2026-04-19`
- O2Physics ROOT executor `ctest --output-on-failure` passed on `2026-05-18`
  after the mixed-event phi-splitting switch was added
- the same date's sandboxed run still produced `/dev/fd/... Operation not
  permitted`, so ROOT-guarded skips remain non-authoritative environment noise
- full real-data equivalence validation against the legacy macro has not yet
  been rerun after the current phi-mapping update
- the new run script is verified as an operator wrapper, but the default OO
  `all` stage has not been executed in this update

## Active Constraints

- ROOT-dependent validation must be run from a fully entered O2Physics
  environment
- `scripts/run_exp_femto_3d.sh` defaults to `config/oo_build_and_fit.toml` and
  executes both `build-cf` and `fit`, so running it without `--stage` writes the
  configured OO outputs
- progress remains controlled by `[build].progress` and `[fit].progress`; in
  `"auto"` mode the ETA line is shown only when `stderr` is attached to a TTY
- sandboxed tool runs that fail during `alienv` bootstrap are not reliable
  evidence for code-level ROOT regressions
- physics-level closure still requires a real-data regression on a known-good
  dataset
- `build.split_mixed_event_by_phi = false` preserves the historical
  phi-integrated mixed-event denominator per centrality/mT group; `true` is the
  opt-in mode for denominators that follow each SE phi slice

## Ledger Convention

- `Exp_femto_3d` adopts `project-state/` as its coordination ledger path
- this project does not use the hidden `.project-state/` path

## Active Worktree Highlights

- `scripts/run_exp_femto_3d.sh` provides the project-local OO run entry with
  runtime re-entry, stage selection, and fit overrides
- build and fit now expose explicit progress-mode control plus Eventgen-style
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

- `project-state/` is now the active adopted coordination ledger for
  `Exp_femto_3d`
- this sync updates the ledger from the earlier ROOT-runtime bootstrap state to
  the current phi-mapping/progress-control development state, the OO run script
  entry, and the Eventgen-style progress rendering update
