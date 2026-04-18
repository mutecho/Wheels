# Tests

## T-001: Non-Sandboxed O2Physics CTest Run After Phi-Mapping Update

- date: 2026-04-19
- environment: `alienv setenv O2Physics/latest-master-o2 -c sh -lc`
- command:

```bash
cd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build &&
ctest --output-on-failure
```

- result: passed
- evidence:
  - `config_parse_validation_test`: passed
  - `slice_catalog_roundtrip_test`: passed
  - `workflow_smoke_test`: passed
- significance: verifies progress-mode parsing, build-side phi-mapping
  persistence, legacy `SliceCatalog` compatibility, and fit-side phi remapping
  overrides in a real ROOT/O2 runtime

## T-002: Sandboxed O2Physics CTest Run On The Same Worktree

- date: 2026-04-19
- environment: sandboxed tool execution
- command: same `ctest --output-on-failure` entry as above
- result: partially executed; config-only test passed, ROOT-backed tests were
  skipped by the runtime guard
- signature:
  - `/dev/fd/... Operation not permitted`
  - incomplete `alienv` bootstrap
- interpretation: still an environment-entry limitation, not authoritative
  evidence of a code regression

## T-003: Direct ROOT Runtime Sanity Check

- date: 2026-04-11
- environment: non-sandboxed O2Physics runtime
- command type: ROOT snippet creating `THnSparseF`
- result: passed
- significance: confirms the local ROOT installation is usable when the O2
  environment is entered correctly
