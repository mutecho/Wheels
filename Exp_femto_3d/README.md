# Exp_femto_3d

Refactored C++17/CMake version of the original `3d_cf_from_exp.cpp` ROOT macro.

## Layout

- `include/exp_femto_3d/`: public types and workflow interfaces
- `src/`: configuration parsing, logging, CF building, and Levy fitting
- `app/`: CLI entry point
- `docs/`: human-facing workflow notes, including the math/physics formula flow
- `config/examples/`: example TOML configurations
- `tests/`: config and workflow verification helpers
- `legacy/`: archived macro reference kept out of the CMake build

## Documentation

- [`docs/数学物理公式流程说明.md`](docs/数学物理公式流程说明.md): current
  math/physics workflow reference from sparse inputs through CF building,
  Levy/Coulomb fitting, and summary outputs

## Build

This project expects ROOT from the ALICE/O2 environment and `toml++` from
Homebrew.

```bash
cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
      -B /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build
cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build
```

Configured build trees still live under `build/`, but all executable targets are
emitted into `bin/`. `compile_commands.json` is exported in the configured build
directory as `build/compile_commands.json`.

## Run

```bash
./bin/exp_femto_3d build-cf --config config/pbpb_build_and_fit.toml
./bin/exp_femto_3d fit --config config/pbpb_build_and_fit.toml
```

Optional fit overrides:

- `--model full|diag`
- `--input-cf-root /absolute/or/relative/path.root`

Coulomb fit modes:

- `[fit].coulomb_mode = "none" | "gamow" | "finite_source"`
- `[fit].finite_source_mode = "fixed_1d" | "iterative_1d"` is only valid with
  `coulomb_mode = "finite_source"`
- older `[fit].use_coulomb = false|true` configs remain accepted and map to
  `none` and `gamow`; when both fields are present they must agree
- finite-source mode uses local CATS support when CMake finds
  `/Users/allenzhou/Research_software/CATS/install` and Homebrew GSL; builds
  without CATS still support `none` and `gamow`

Independent fit-report output:

- `[output].fit_report_directory`
- `[output].fit_report_root_name`
- these fields control the standalone ROOT report file produced by `fit`; when
  `fit_report_directory` is omitted, it defaults to `output_directory`

TOML progress control:

- `[build].progress`
- `[fit].progress`
- accepted values: `true`, `false`, or `"auto"` (also `"enabled"` / `"disabled"`)
- default when omitted: `"auto"`; the bar only appears when `stderr` is attached to a TTY
- enabled progress lines include the stage label, percentage, an activity frame,
  and an ETA computed from completed slices; a heartbeat refreshes the line once
  per second during long ROOT operations

Phi mapping control:

- `[build].map_pair_phi_to_symmetric_range` controls how build writes `display_phi_*`
  into `meta/SliceCatalog`
- `[fit].map_pair_phi_to_symmetric_range` is optional; when omitted, fit follows the
  input CF metadata from `SliceCatalog`
- when `[fit].map_pair_phi_to_symmetric_range` is set explicitly, fit reinterprets the
  stored slices from `raw_phi_*`, so the summary `R2_vs_phi` graphs can use either the
  raw `[0, pi]` range or the symmetric `[-pi/2, pi/2]` range without rebuilding the CF
- if the current config's `[build]` mapping flag disagrees with the input CF metadata,
  fit warns and trusts the input CF file

Mixed-event denominator control:

- `[build].split_mixed_event_by_phi` defaults to `false`
- `false` preserves the historical behavior: one mixed-event denominator is
  projected over all phi bins for each centrality/mT group, then reused by the
  phi-integrated and per-phi SE slices
- `true` projects the mixed-event denominator with the same phi range as the SE
  slice; the phi-integrated slice still uses the full phi range
- this switch is independent of `map_pair_phi_to_symmetric_range`, which only
  controls stored/displayed phi coordinates for downstream summaries

## Output Contract

- `build-cf` writes `meta/SliceCatalog` plus `slices/<slice_id>/...`; the
  catalog records both `build_uses_symmetric_phi_range` and
  `split_mixed_event_by_phi`, with older catalogs defaulting the latter to
  `false`
- `fit` reads `meta/SliceCatalog` and writes `meta/FitCatalog`,
  `meta/CoulombKernelCatalog`, `fits/<slice_id>/...`, `summary/R2_vs_phi/...`,
  and `fit_summary.tsv`
- `fit` also writes the standalone report ROOT file configured by
  `fit_report_directory` and `fit_report_root_name`; this report mirrors
  `meta/FitCatalog` and `summary/R2_vs_phi/...`, adds
  `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`, and adds
  `eps_vs_mt/<cent>/epsf_vs_mt(_canvas)` summaries

## Test Notes

`ctest` always runs config-only tests. ROOT integration tests are guarded by a
runtime probe and are skipped when the local ROOT runtime cannot safely create,
persist, and project a `THnSparseF`.

For agent-oriented diagnosis of sandboxed ROOT/O2 failures, see
[`ROOT_RUNTIME_AGENT_NOTE.md`](./ROOT_RUNTIME_AGENT_NOTE.md).
