# Handoff

## Latest Completed Work

- scaffolded `Exp_femto_1d` as a standalone CMake package
- added public `Types/Config/Logging/Workflow/CatsModel` interfaces
- implemented structured `build-cf` output with `SliceCatalog`
- added fit-side `FitCatalog`, per-slice fit artifacts, summary histograms, and TSV writing
- created example TOML configs, README, ROOT runtime note, and `project-state/`
- passed non-sandboxed O2Physics configure/build/`ctest --output-on-failure` on `2026-04-20`

## Current Risk

- real-data equivalence to the legacy macro is still open
- fit stability on real slices still needs inspection

## Recommended Next Action

1. Execute one real-data `build-cf` regression against the legacy macro.
2. Run one real-data `fit` on a `MinBias` slice and one EP-differential slice.
3. Inspect `FitCanvas`, `Baseline`, and `PureFemto` outputs.

## Recommended Commands

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/bin/exp_femto_1d build-cf \
    --config /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/pbpb_build_and_fit.toml &&
  /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/bin/exp_femto_1d fit \
    --config /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/pbpb_build_and_fit.toml
'
```
