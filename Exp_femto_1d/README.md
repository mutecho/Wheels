# Exp_femto_1d

Refactored `C++17 + CMake` version of the legacy `get_cf_from_exp.cpp` macro,
plus a first-pass 1D `pi-pi` CATS fit workflow.

## Scope

This project keeps the original experimental input structure and event-plane
slice definitions, but moves the logic into a maintainable two-stage workflow:

1. `build-cf`
   - reads `THnSparseF` same-event and mixed-event inputs
   - builds `MinBias`, `InPlane`, and `OutOfPlane` 1D slices
   - writes an explicit `meta/SliceCatalog` plus structured slice directories
2. `fit`
   - reads `meta/SliceCatalog`
   - fits each selected slice with a fixed `baseline polynomial x CATS`
     `pi-pi` model
   - writes `meta/FitCatalog`, per-slice fit artifacts, region summaries, and
     `fit_summary.tsv`

Version `v1` is intentionally narrow:

- like-sign `pi-pi` only
- one fixed model: `cats-pipi-poly`
- no legacy flat-name compatibility layer
- no source-family switching
- fixed event-plane regions: `MinBias / InPlane / OutOfPlane`

## Layout

- `include/exp_femto_1d/`: public types, config, logging, workflow, and CATS model interfaces
- `src/`: config parsing, logging, build/fit workflow, and CATS implementation
- `app/`: CLI entry point
- `config/examples/`: minimal and commented TOML examples
- `scripts/`: project-local run helpers for config-driven workflows
- `tests/`: config-only tests plus ROOT-guarded workflow/CATS smoke tests
- `legacy/`: archived macro reference kept out of the CMake build
- `project-state/`: coordination ledger for status, decisions, tests, and follow-up work

## Build

This project expects:

- ROOT from the ALICE/O2Physics runtime
- `toml++` from Homebrew
- installed CATS headers and `libCATS.dylib`

The CMake build prefers `CATS_ROOT` and otherwise falls back to
`/Users/allenzhou/Research_software/CATS/install`.

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_1d \
        -B /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/build &&
  cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/build
'
```

Configured build trees live under `build/`, while executable targets are
emitted into `bin/`.

## Run

The project-local helper runs the config-driven workflow from any current
directory. By default, it runs only `build-cf` with
`config/pbpb_build_and_fit.toml`; use `--stage all` when the slower fit stage
should run immediately afterward. Direct shell invocations automatically
re-enter `O2Physics/latest-master-o2` with `alienv` when the current shell does
not already provide a complete ROOT runtime:

```bash
/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/scripts/run_exp_femto_1d.sh
```

Common variants:

```bash
./scripts/run_exp_femto_1d.sh --config config/pbpb_build_and_fit.toml
./scripts/run_exp_femto_1d.sh --stage all --config config/pbpb_build_and_fit.toml
./scripts/run_exp_femto_1d.sh --stage build-cf --config config/pbpb_build_and_fit.toml
./scripts/run_exp_femto_1d.sh --stage fit --config config/pbpb_build_and_fit.toml
./scripts/run_exp_femto_1d.sh --stage fit --config config/pbpb_build_and_fit.toml --input-cf-root /path/to/cf.root
```

The helper preserves the underlying CLI contract:

```bash
./bin/exp_femto_1d build-cf --config config/pbpb_build_and_fit.toml
./bin/exp_femto_1d fit --config config/pbpb_build_and_fit.toml
./bin/exp_femto_1d fit --config config/pbpb_build_and_fit.toml --input-cf-root /path/to/cf.root
```

On success, the CLI prints compact stage summaries:

- `build-cf`: `stored_slices`, `skipped_zero_mixed_event_groups`,
  `skipped_zero_mixed_event_slices`, `skipped_zero_same_event_slices`
- `fit`: `catalog_slices`, `selected_slices`, `fitted_slices`, `skipped_missing_objects`, `skipped_failed_fits`

On argument or runtime errors, the CLI prints `[error] <message>`, then usage,
and exits non-zero.

## Config Schema

Top-level TOML tables:

- `[input]`
- `[output]`
- `[build]`
- `[fit]`
- `[[bins.centrality]]`
- `[[bins.mt]]`
- `[[fit_selection.centrality]]`
- `[[fit_selection.mt]]`

Key semantics:

- ROOT input and stored histograms use `GeV/c`
- `PiPiCatsModel` converts the fit-side `k*` to `MeV/c` before calling CATS
- missing `fit_selection` falls back to the full build bin lists
- output file extensions are normalized to `.root` and `.tsv`
- `build.split_mixed_event_by_phi` defaults to `false`; `false` keeps one
  MinBias mixed-event denominator per centrality/mT group, while `true`
  projects mixed-event denominators with the same event-plane region as each
  same-event slice
- `build.cf_by_mt_show_markers` defaults to `false`; set it to `true` to draw
  markers on `CFByMtCanvas` trend lines
- build/fit progress accepts `true`, `false`, `"auto"`, `"enabled"`, or `"disabled"`

See:

- [config/examples/minimal.toml](/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/examples/minimal.toml)
- [config/examples/exp_femto_1d.example.toml](/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/examples/exp_femto_1d.example.toml)
- [config/pbpb_build_and_fit.toml](/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/pbpb_build_and_fit.toml)

## Output Contract

`build-cf` writes:

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw1d`
- `slices/<slice_id>/ME_raw1d`
- `slices/<slice_id>/CF1D`
- `cent_slices/<cent_id>/<region_name>/CFByMtCanvas`

`meta/SliceCatalog` records `split_mixed_event_by_phi` for every slice. Older
CF files without this branch are read as `false`.

`fit` writes:

- `meta/FitCatalog`
- `fits/<slice_id>/DataCF`
- `fits/<slice_id>/FitFunction`
- `fits/<slice_id>/Baseline`
- `fits/<slice_id>/PureFemto`
- `fits/<slice_id>/FitCanvas`
- `summary/by_region/<group_id>/SourceSize`
- `summary/by_region/<group_id>/Norm`
- `summary/by_region/<group_id>/P1`
- `summary/by_region/<group_id>/P2`
- `fit_summary.tsv`

## Test Notes

`ctest` always runs the config-only test. ROOT-backed tests are guarded by
`tests/run_root_guard.sh`, which skips with code `77` when the local ROOT
runtime cannot safely create, persist, and project a toy `THnSparseF`.

Because the authoritative ROOT runtime comes from O2Physics, final build and
workflow validation should be run from a fully entered non-sandboxed O2
environment. See
[ROOT_RUNTIME_AGENT_NOTE.md](/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/ROOT_RUNTIME_AGENT_NOTE.md).
