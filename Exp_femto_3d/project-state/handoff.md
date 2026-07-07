# Handoff

## Latest Durable Handoff

- completed:
  - added optional `[build].split_mixed_event_by_qn`; default `false` preserves
    existing qn-specific SAME slices with qn-integrated ME denominators
  - required `split_same_event_by_qn = true` before
    `split_mixed_event_by_qn = true` is accepted, so ME qn splitting only
    applies when qn-specific SE slices are written
  - updated `RunBuildCf()` so ME projections apply the current qn selection only
    when the new switch is enabled; qn-all slices remain integrated over the
    full qn axis
  - persisted `split_mixed_event_by_qn` in `meta/SliceCatalog`, with legacy
    catalog reads defaulting the missing branch to `false`
  - added config validation coverage and extended `slice_catalog_roundtrip_test`
    so toy qn1 `ME_raw3d` under the new switch is compared against the old
    qn-integrated denominator
  - documented the new switch in README, example TOMLs,
    `docs/数学物理公式流程说明.md`, and `project-state/`
  - restored
    `config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml` to explicit
    `split_same_event_by_qn = true` because the current test contract and qn
    production baseline expect that setting enabled
  - ran `scripts/cmake.sh`; it configured CATS support and rebuilt changed
    targets successfully
  - ran `ctest --test-dir build --output-on-failure`; it reported `100% tests
    passed`, with ROOT-backed tests skipped by the local guard
  - ran `git diff --check`; it passed
  - attempted direct O2Physics ROOT executor validation of
    `bin/slice_catalog_roundtrip_test`; sandbox execution returned
    `STATUS: ESCALATION_REQUIRED`, and the required escalated rerun was
    auto-rejected by the platform usage limit

## Next Recommended Owner Action

- to enable qn-matched ME denominators in a qn-split build, use:

```toml
[build]
split_same_event_by_qn = true
split_mixed_event_by_qn = true
```

- rerun `bin/slice_catalog_roundtrip_test` through a usable O2Physics ROOT
  runtime to execute the new qn1 `ME_raw3d` integral assertion

## Previous Durable Handoff

- completed:
  - added optional `[fit.parameters.<name>]` TOML subtables for the main 3D
    Levy fit parameters: `norm`, `lambda`, `rout2`, `rside2`, `rlong2`,
    `routside2`, `routlong2`, `rsidelong2`, `alpha`, and `baseline_q2`
  - each supported parameter can override `initial`, `min`, and `max`; only
    `lambda` and `alpha` support `fixed_value`
  - kept all legacy defaults when the new subtables are omitted, including the
    existing q2-baseline default bounds derived from `fit_q_max`
  - rejected ambiguous or ignored configurations: unknown names/fields,
    non-finite values, one-sided or invalid limits, unsupported fixed values,
    out-of-bounds fixed values, `lambda` overrides with
    `use_core_halo_lambda=false`, and `baseline_q2` overrides with
    `use_q2_baseline=false`
  - routed both diag/full `TF3` builders and the PML `TMinuit` setup through the
    same effective parameter configuration so fixed `lambda`/`alpha` work in
    both chi2 and PML paths
  - updated README, example TOMLs, the formula workflow document, tests, and
    `project-state/`
  - ran `scripts/cmake.sh`; it configured CATS support and rebuilt changed
    targets successfully
  - ran `ctest --test-dir build --output-on-failure`; it reported `100% tests
    passed`, with ROOT-backed tests skipped by the local guard
  - reran ROOT-backed tests in a clean non-sandboxed O2Physics runtime:
    `bin/slice_catalog_roundtrip_test` and `bin/workflow_smoke_test` both
    passed

## Previous Recommended Owner Action

- use `[fit.parameters.<name>]` only for intentional fit-control changes; leave
  it omitted to preserve the old production defaults
- for fixed source-shape scans, prefer:

```toml
[fit.parameters.lambda]
fixed_value = 0.65

[fit.parameters.alpha]
fixed_value = 1.20
```

- keep running ROOT-backed validation through a clean O2Physics runtime when
  changing fit math, parameter limits, or ROOT output semantics

## Previous Durable Handoff

- completed:
  - extended
    `config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml` for the Wenya
    PbPb LHC23 input so it preserves qn-all and adds qn1/qn2/qn3 SAME-side
    slices
  - kept the original `Exp_femto_3d` output structure:
    `meta/SliceCatalog`, `slices/<slice_id>/...`, detailed fit ROOT, TSV
    summary, `summary/R2_vs_phi`, and standalone report ROOT
  - added `[build].split_same_event_by_qn` and `[[bins.qn]]`; qn-specific
    group ids append `__qn1`, `__qn2`, or `__qn3`, while qn-all keeps the
    historical cent/mT group id
  - kept MIXED qn integrated; qn splitting applies to SAME only
  - used the existing PbPb fit settings: full model, legacy Gamow Coulomb via
    `use_coulomb = true`, core-halo lambda, q2 baseline, PML, and
    `fit_q_max = 0.15`
  - added config-parse coverage, qn catalog roundtrip coverage, and qn metadata
    smoke coverage for the new output schema
  - ran `scripts/cmake.sh`; it configured CATS support and rebuilt all targets
  - ran `ctest --test-dir build --output-on-failure`; it reported `100% tests
    passed`, with ROOT-backed tests then run directly through the O2Physics
    executor
  - ran `bin/slice_catalog_roundtrip_test` and `bin/workflow_smoke_test`
    through the O2Physics executor; both returned `PRIMARY_OK`
  - ran production `build-cf`; it returned `PRIMARY_OK` and stored `1040`
    slices with no skipped groups or slices
  - ran production `fit`; it returned `PRIMARY_OK` and fitted `468/468`
    selected slices with no missing objects or raw histograms
  - inspected produced ROOT files: `SliceCatalog=1040`, qn_all/qn1/qn2/qn3
    each `260` build slices, `FitCatalog=468`, qn_all/qn1/qn2/qn3 each `117`
    fitted slices, report/catalog/summary objects present, TSV line count
    `469`

## Previous Recommended Owner Action

- use the updated config directly for this Wenya qn-all plus qn1/qn2/qn3
  production path:

```bash
bin/exp_femto_3d build-cf --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml
bin/exp_femto_3d fit --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml
```

- when future work changes sparse axes, CF normalization, phi mapping,
  `SliceCatalog`, Levy/Coulomb formulas, fit metadata, or summary-output
  semantics, update `docs/数学物理公式流程说明.md` in the same pass

## Older Formula Handoff

- completed:
  - added `docs/数学物理公式流程说明.md` as the current formula workflow
    reference for `Exp_femto_3d`
  - documented the sparse-axis contract, build-side SE/ME normalization,
    `CF3D` construction, phi coordinate mapping, `SliceCatalog`, diag/full
    Levy formulas, optional PML objective, Gamow and finite-source Coulomb
    branches, fit output catalogs, `R2_vs_phi`, report canvases, and
    `epsf_vs_mt`
  - linked the new document from `README.md`
  - synced `project-state/current-status.md` and `project-state/changelog.md`
    for this docs-only update
  - verified by reading the current implementation and docs; no analysis code,
    configs, build files, or runtime outputs were changed

## Older Finite-Source Handoff

- completed:
  - implemented the `docs/plan/fit_finite_coul.md` finite-source Coulomb fit
    path
  - added explicit `fit.coulomb_mode = "none"|"gamow"|"finite_source"` and
    `fit.finite_source_mode = "fixed_1d"|"iterative_1d"` parsing
  - kept legacy `fit.use_coulomb` compatibility only for unambiguous none/Gamow
    configs and reject conflicting mixed legacy/new settings
  - added optional CATS/GSL detection to CMake, with
    `EXP_FEMTO_3D_ENABLE_CATS=OFF` available for no-CATS validation
  - wired finite-source fitting through CATS-backed one-dimensional kernel
    tables keyed by centrality/mT and seeded from the corresponding `phi_all`
    slice
  - implemented fixed and one-pass iterative source-radius flows before final
    selected-slice fitting and artifact writing
  - added fail-fast handling for invalid in-table CATS kernel values; only
    above-table high-k evaluation falls back to unity
  - preserved `usesCoulomb` and added Coulomb mode, finite-source mode, and
    finite-source radius metadata in TSV and `meta/FitCatalog`
  - added `meta/CoulombKernelCatalog` to both detailed fit ROOT and standalone
    report ROOT outputs
  - updated README and example configs for the new public mode names
  - added/expanded config parsing, workflow smoke, and kernel validation tests
  - reran O2Physics ROOT executor configure/build/`ctest --output-on-failure`
    on `2026-06-21` for the default CATS-enabled build; all five registered
    tests passed with `PRIMARY_OK`
  - reran the same matrix on `2026-06-21` with
    `-DEXP_FEMTO_3D_ENABLE_CATS=OFF`; all five registered tests passed with
    `PRIMARY_OK`
  - hardened `scripts/cmake.sh` so the default local build stays incremental,
    refreshes stale ROOT CMake cache/link rules when the active ROOT runtime
    changes, and verifies CATS linkage when the generated link rule expects
    CATS
  - reran `scripts/cmake.sh` through the O2Physics ROOT executor on
    `2026-06-23`; it returned `PRIMARY_OK`, refreshed cached
    `ROOT/v6-36-10-alice1-local7` links to active
    `ROOT/v6-36-10-alice1-local8`, and a second run performed no compile/link
    work
  - `otool` checks on `bin/exp_femto_3d` showed the active ROOT `LC_RPATH`,
    `libCATS`, and GSL
  - reran `ctest --test-dir build --output-on-failure` on `2026-06-23`; all
    five registered tests passed with `PRIMARY_OK`
  - reran `git diff --check`; it passed

## Older Recommended Owner Action

- run a real-data regression on a known-good OO/PbPb input set with
  `fit.coulomb_mode = "finite_source"`
- inspect `meta/CoulombKernelCatalog`, TSV `finiteSourceRadiusFm`, and report
  ROOT source-parameter canvases for representative centrality/mT groups
- compare the finite-source results against the existing Gamow baseline before
  treating the new mode as production physics default
- cover both `fit.finite_source_mode = "fixed_1d"` and `"iterative_1d"` during
  that real-data regression
- continue to keep no-CATS builds in the smoke matrix so finite-source requests
  fail explicitly on machines without CATS
- use `scripts/cmake.sh` for the default local build before running
  finite-source configs; set `EXP_FEMTO_3D_CLEAN_FIRST=1` only when a full
  clean rebuild is deliberately needed
- keep treating sandbox-only `alienv` failures as environment noise unless a
  non-sandboxed O2Physics rerun reproduces them

## Previous Suggested Next Commands

```bash
/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh \
  --stage fit \
  --input-cf-root /path/to/existing_cf.root
```

Set the fit mode in the TOML before the finite-source regression:

```toml
[fit]
coulomb_mode = "finite_source"
finite_source_mode = "fixed_1d"      # or "iterative_1d"
```

For local validation after edits:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command '/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/cmake.sh && ctest --test-dir /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build --output-on-failure'
```

Then run a build-cf / fit comparison on a previously validated real input set,
once against the Gamow baseline and once against the finite-source mode.
