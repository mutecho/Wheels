# Handoff

## Latest Durable Handoff

- completed:
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
- cover both phi conventions during that regression: follow-input mapping and
  explicit fit-side override
- cover both ME denominator modes during that regression:
  `build.split_mixed_event_by_phi = false` and `true`
- keep treating sandbox-only `alienv` failures as environment noise unless a
  non-sandboxed O2Physics rerun reproduces them

## Suggested Next Commands

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build &&
  ctest --output-on-failure
'
```

Then run a build-cf / fit comparison on a previously validated real input set,
once with the input CF mapping followed as-is and once with an explicit fit-side
phi override.
