# Exp_femto_3d

Refactored C++17/CMake version of the original `3d_cf_from_exp.cpp` ROOT macro.

## Layout

- `include/exp_femto_3d/`: public types and workflow interfaces
- `src/`: configuration parsing, logging, CF building, and Levy fitting
- `app/`: CLI entry point
- `config/examples/`: example TOML configurations
- `tests/`: config and workflow verification helpers
- `legacy/`: archived macro reference kept out of the CMake build

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
./bin/exp_femto_3d build-cf --config config/examples/pbpb_build_and_fit.toml
./bin/exp_femto_3d fit --config config/examples/pbpb_build_and_fit.toml
```

Optional fit overrides:

- `--model full|diag`
- `--input-cf-root /absolute/or/relative/path.root`

TOML progress control:

- `[build].progress`
- `[fit].progress`
- accepted values: `true`, `false`, or `"auto"` (also `"enabled"` / `"disabled"`)
- default when omitted: `"auto"`; the bar only appears when `stderr` is attached to a TTY

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

## Output Contract

- `build-cf` writes `meta/SliceCatalog` plus `slices/<slice_id>/...`
- `fit` reads `meta/SliceCatalog` and writes `meta/FitCatalog`,
  `fits/<slice_id>/...`, `summary/R2_vs_phi/...`, and `fit_summary.tsv`

## Test Notes

`ctest` always runs config-only tests. ROOT integration tests are guarded by a
runtime probe and are skipped when the local ROOT runtime cannot safely create,
persist, and project a `THnSparseF`.

For agent-oriented diagnosis of sandboxed ROOT/O2 failures, see
[`ROOT_RUNTIME_AGENT_NOTE.md`](./ROOT_RUNTIME_AGENT_NOTE.md).
