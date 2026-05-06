# Test Ledger

## Implemented Tests

- `config_parse_validation_test`
  - covers required fields, extension normalization, progress parsing, invalid
    fit limits, duplicate bins, invalid `fit_selection`, and example config parsing
- `slice_catalog_roundtrip_test`
  - builds a toy 4D sparse input and verifies `SliceCatalog` metadata plus
    structured slice paths
- `workflow_smoke_test`
  - runs toy `build-cf` and `fit`, then checks `SliceCatalog`, `FitCatalog`,
    per-slice objects, summaries, and TSV headers
- `cats_fit_smoke_test`
  - fits a synthetic `baseline x CATS` histogram and checks `p0` plus source-size recovery
- `root_runtime_probe`
  - confirms that the local ROOT runtime can create, store, reopen, and project
    a toy `THnSparseF`

## Current Validation State

- `2026-04-20` non-sandboxed O2Physics configure/build/`ctest --output-on-failure`
  passed all 4 registered tests
- `2026-05-06` `bash -n scripts/run_exp_femto_1d.sh`,
  `scripts/run_exp_femto_1d.sh --help`, and a non-ROOT dispatch check with
  `--binary /bin/echo` passed for the project-local run helper
- `2026-05-06` reproduced direct shell failure as missing O2Physics ROOT
  runtime (`ROOT_DYN_PATH`/PCM failures), then confirmed direct
  `scripts/run_exp_femto_1d.sh --stage build-cf` auto-enters
  `O2Physics/latest-master-o2` and writes 75 slices
- `2026-05-06` `/bin/echo` dispatch check confirmed `--stage fit`,
  `--config`, and `--input-cf-root` survive the auto-reexec path
- `2026-05-06` `/bin/echo` dispatch check confirmed default helper execution
  runs only `build-cf`; explicit `--stage all` dispatches `build-cf` before
  `fit`
- `2026-05-06` fixed and smoke-checked the zero-argument path under `set -u`
  so Bash 3.2 does not expand an empty `original_args` array
- `2026-05-06` `cmake --build Exp_femto_1d/build` passed after detaching
  internal `THnSparse` projections immediately from ROOT directories
- `2026-05-06` direct `scripts/run_exp_femto_1d.sh --stage build-cf` completed
  with 75 stored slices and no `TROOT::Append` replacement warnings
- `2026-05-06` non-sandboxed
  `alienv setenv O2Physics/latest-master-o2 -c sh -lc 'cd .../Exp_femto_1d/build && ctest --output-on-failure'`
  passed all 4 registered tests
- `2026-05-06` changed the run helper default to `build-cf` so direct script
  execution returns after CF construction instead of silently entering the long
  fit stage
- `2026-05-06` real zero-argument
  `scripts/run_exp_femto_1d.sh` completed `build-cf` with 75 stored slices and
  exited with status 0
- real-data regression against the legacy macro is still pending

## Required Follow-up Validation

- run direct `scripts/run_exp_femto_1d.sh --stage all --config
  config/pbpb_build_and_fit.toml` after confirming fit runtime is acceptable
- compare new `build-cf` against `legacy/get_cf_from_exp.cpp` on a known-good dataset
- inspect at least one `MinBias` and one EP-differential real-data fit canvas
