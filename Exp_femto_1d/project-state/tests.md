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
- real-data regression against the legacy macro is still pending

## Required Follow-up Validation

- run full configure/build/ctest under O2Physics
- compare new `build-cf` against `legacy/get_cf_from_exp.cpp` on a known-good dataset
- inspect at least one `MinBias` and one EP-differential real-data fit canvas
