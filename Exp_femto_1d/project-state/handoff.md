# Handoff

## Latest Completed Work

- scaffolded `Exp_femto_1d` as a standalone CMake package
- added public `Types/Config/Logging/Workflow/CatsModel` interfaces
- implemented structured `build-cf` output with `SliceCatalog`
- added fit-side `FitCatalog`, per-slice fit artifacts, summary histograms, and TSV writing
- created example TOML configs, README, ROOT runtime note, and `project-state/`
- passed non-sandboxed O2Physics configure/build/`ctest --output-on-failure` on `2026-04-20`
- added `scripts/run_exp_femto_1d.sh` as the config-driven run helper on `2026-05-06`
- fixed direct run-helper execution to re-enter `O2Physics/latest-master-o2`
  automatically when the caller has not already entered a complete ROOT runtime
- added `build-cf` `cent_slices/<cent_id>/<region_name>/CFByMtCanvas` overlays
  on `2026-05-08` so each centrality and event-plane region shows all mT-bin
  CF curves together
- changed `CFByMtCanvas` to line-only output by default; use
  `build.cf_by_mt_show_markers = true` only when marker symbols are desired
- added `build.split_mixed_event_by_phi` on `2026-05-18`; default `false`
  preserves the integrated MinBias mixed-event denominator, while `true`
  projects mixed events in the same event-plane region as each SE slice
- persisted the ME-splitting choice in `meta/SliceCatalog` and passed local
  build plus O2Physics ROOT executor `ctest --output-on-failure`

## Current Risk

- real-data equivalence to the legacy macro is still open
- real-data comparison of integrated-ME and split-ME modes is still open
- fit stability on real slices still needs inspection

## Recommended Next Action

1. Execute one real-data `build-cf` regression against the legacy macro with
   `build.split_mixed_event_by_phi = false`.
2. Inspect the new `cent_slices/*/*/CFByMtCanvas` overlays for at least one
   centrality bin.
3. Repeat the same dataset with `build.split_mixed_event_by_phi = true` and
   compare the EP-differential denominators.
4. Run one real-data `fit` on a `MinBias` slice and one EP-differential slice.
5. Inspect `FitCanvas`, `Baseline`, and `PureFemto` outputs.

## Recommended Commands

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/scripts/run_exp_femto_1d.sh \
    --stage all \
    --config /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/config/pbpb_build_and_fit.toml
'
```
