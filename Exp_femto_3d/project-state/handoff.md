# Handoff

## Latest Durable Handoff

- completed:
  - upgraded `build-cf` and `fit` progress rendering to show stage label,
    percent, activity frame, and ETA
  - added a one-second heartbeat so enabled progress output keeps refreshing
    during long ROOT operations between completed slices
  - added `progress_render_test` and linked `exp_femto_3d_core` with
    `Threads::Threads`
  - reran O2Physics ROOT executor configure/build/`ctest --output-on-failure`
    on `2026-06-10`; all four registered tests passed with `PRIMARY_OK`
  - added `scripts/run_exp_femto_3d.sh` as the project-local OO run entry
    defaulting to `config/oo_build_and_fit.toml`
  - the script runs `build-cf` followed by `fit` by default, supports
    `--stage build-cf|fit|all`, and forwards fit-only `--model full|diag` plus
    `--input-cf-root` overrides
  - the script re-enters `O2Physics/latest-master-o2` through `alienv` when the
    ROOT runtime is not already active
  - verified the new script with `bash -n` and `--help` on `2026-06-10`; the
    full OO real-data build/fit was not rerun because it writes production
    output paths from the TOML config
  - fixed the no-argument `alienv` re-entry path so macOS Bash 3.2 `set -u`
    does not treat an empty preserved-argument array as unbound
  - added standalone fit report ROOT output controlled by
    `[output].fit_report_directory` and `[output].fit_report_root_name`
  - report ROOT now includes `meta/FitCatalog`, mirrored
    `summary/R2_vs_phi/...`, Eventgen-style
    `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`, and
    `eps_vs_mt/<cent>/epsf_vs_mt(_canvas)`
  - updated active/example TOML configs and README output contract for the new
    report path
  - reran O2Physics ROOT executor `ctest --output-on-failure` on `2026-06-10`;
    all three registered tests passed
  - added explicit TOML progress-mode parsing for build and fit
  - persisted build-side phi mapping state into `meta/SliceCatalog`
  - taught fit to follow input CF phi metadata by default or override it from
    stored `raw_phi_*` coordinates
  - added backward-compatible inference for legacy `SliceCatalog` trees that do
    not yet contain `build_uses_symmetric_phi_range`
  - added `build.split_mixed_event_by_phi` with default `false` and an opt-in
    mode where ME denominators follow each SE phi slice
  - persisted `split_mixed_event_by_phi` in `meta/SliceCatalog` with legacy
    default `false`
  - extended config, catalog roundtrip, and workflow smoke coverage for the new
    phi/progress and split-ME semantics
  - reran `ctest --output-on-failure` in a non-sandboxed O2Physics environment
    on `2026-04-19`; all three registered tests passed
  - reran O2Physics ROOT executor `ctest --output-on-failure` on `2026-05-18`;
    all three registered tests passed

## Next Recommended Owner Action

- run the real-data regression comparison between the refactored executable and
  the legacy macro on a known-good dataset
- include the standalone fit report ROOT file in that regression by checking
  the `source_parameters` overview canvases and `eps_vs_mt` graphs on a real
  multi-mT dataset
- cover both phi conventions during that regression: follow-input mapping and
  explicit fit-side override
- cover both ME denominator modes during that regression:
  `build.split_mixed_event_by_phi = false` and `true`
- keep treating sandbox-only `alienv` failures as environment noise unless a
  non-sandboxed O2Physics rerun reproduces them

## Suggested Next Commands

```bash
/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh

/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh \
  --stage fit \
  --input-cf-root /path/to/existing_cf.root
```

For regression coverage:

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build &&
  ctest --output-on-failure
'
```

Then run a build-cf / fit comparison on a previously validated real input set,
once with the input CF mapping followed as-is and once with an explicit fit-side
phi override.
