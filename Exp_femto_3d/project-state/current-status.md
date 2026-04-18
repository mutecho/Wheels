# Current Status

## Task Snapshot

- scope: sync `project-state/` with the current `3d_cf_from_exp` refactor work
  around phi-mapping persistence/override and explicit progress-mode control
- current conclusion: the active worktree has advanced beyond the earlier
  ROOT-runtime-diagnosis state; the current implementation now treats phi
  mapping as durable file metadata and lets fit follow or override it
- primary evidence:
  - `[build].progress` and `[fit].progress` now parse `true`, `false`, and
    `"auto"`
  - build writes `build_uses_symmetric_phi_range` into `meta/SliceCatalog`
  - fit follows input CF metadata by default and can override it via
    `[fit].map_pair_phi_to_symmetric_range`
  - legacy `SliceCatalog` trees without the new branch are still readable
    through raw/display phi inference
  - `2026-04-19` non-sandboxed O2Physics `ctest --output-on-failure` passed
    all registered tests

## Verification Status

- verification_status: locally verified for config parsing and ROOT-backed
  smoke coverage
- project_state_sync_status: written

Reason:

- authoritative local test execution for the current worktree passed in a clean
  non-sandboxed O2Physics environment on `2026-04-19`
- the same date's sandboxed run still produced `/dev/fd/... Operation not
  permitted`, so ROOT-guarded skips remain non-authoritative environment noise
- full real-data equivalence validation against the legacy macro has not yet
  been rerun after the current phi-mapping update

## Active Constraints

- ROOT-dependent validation must be run from a fully entered O2Physics
  environment
- sandboxed tool runs that fail during `alienv` bootstrap are not reliable
  evidence for code-level ROOT regressions
- physics-level closure still requires a real-data regression on a known-good
  dataset

## Ledger Convention

- `Exp_femto_3d` adopts `project-state/` as its coordination ledger path
- this project does not use the hidden `.project-state/` path

## Active Worktree Highlights

- build and fit now expose explicit progress-mode control
- `SliceCatalog` carries file-level build phi mapping metadata for downstream
  consumers
- fit can reinterpret stored `raw_phi_*` slices into either raw `[0, pi]` or
  symmetric `[-pi/2, pi/2]` summary coordinates without rebuilding CF files
- legacy catalogs remain readable through mapping-state inference

## Coordination Ledger State

- `project-state/` is now the active adopted coordination ledger for
  `Exp_femto_3d`
- this sync updates the ledger from the earlier ROOT-runtime bootstrap state to
  the current phi-mapping/progress-control development state
