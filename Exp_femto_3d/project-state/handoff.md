# Handoff

## Latest Durable Handoff

- completed on 2026-09-04:
  - added `profile_only`, process workers, explicit worker count, HESSE policy,
    checkpoints, resume digests, and `fit --profile-estimate-only`
  - made `profile_only` bypass detailed-fit ROOT, report ROOT, TSV, and
    `FitCatalog` writes entirely
  - used `posix_spawn` rather than post-ROOT `fork`; workers own complete
    Coulomb groups and final merge is catalog-validated and atomic
  - added ordered `PMLBinCache` plus a once-per-slice old/new objective check
  - added attempt FCN/timing/HESSE/error-validity metadata and process execution
    metadata
  - added scout/focused configs and `scripts/run_exp_femto_3d_PROFILE.sh`
- verified:
  - O2Physics executor build and full CTest passed 7/7 with `PRIMARY_OK`
  - workflow smoke covers output isolation, process chunk merge, matching
    resume, digest mismatch rejection, and preservation of the last final ROOT
  - real scout estimate-only passed: 14 slices, 1,386 maximum attempts, no
    likelihood-slice evaluations
  - ROOT inspector returned `PRIMARY_OK` / `STATUS: OK`; toy canvas was viewed
- intentionally pending:
  - real OO scan and 1/2/4/6-worker scaling/RSS/swap benchmark
  - Minuit2 thread backend; configuration is recognized but runtime execution is
    rejected until the specified numerical A/B gate passes

## Previous Durable Handoff

- completed:
  - added default-off `[fit.profile_likelihood]` support inside the existing
    `fit` command for multiple serial 1D/2D scans over listed exact slices or
    all `fit_selection` slices
  - extracted the original PML calculation into the shared
    `EvaluatePMLObjective()` path and consolidated nominal/profile Minuit setup
    in `RunPMLMinimization()` without changing the statistic or physics model
  - added nominal and bidirectional-neighbor starts, complete attempt retention,
    explicit point statuses, one regular interior refinement pass, and a
    per-slice global reference recomputation
  - froze finite-source profile evaluation to the group kernel produced by the
    nominal workflow and recorded its metadata in the profile catalog
  - added the independent profile ROOT schema and diagnostic 1D/2D displays;
    production fits are never replaced by a lower profile minimum
  - added config, pure-driver, and ROOT-guarded smoke coverage, and synchronized
    README, public example configuration, formula documentation, and this ledger
  - preserved the user's local `scripts/run_exp_femto_3d.sh` edit and the two
    untracked OO configurations
- verified:
  - standalone ROOT-independent profile driver compiled and passed
  - `git diff --check` passed
- blocked verification:
  - the required O2Physics ROOT executor cannot enter the runtime in the
    sandbox (`/dev/fd/... Operation not permitted`)
  - the user authorized escalation, but the approval backend rejected the
    escalated executor call with HTTP 403; therefore the full ROOT build/CTest,
    direct skipped-test reruns, ROOT schema inspection, and toy canvas visual QA
    remain pending
  - no real OO data run or Type I–VI physical attribution was performed

## Previous Durable Handoff

- completed:
  - resolved the pull conflict by composing remote qn-ME splitting and Levy
    parameter controls with local phi/mT rebin; no feature side was discarded
  - retained explicit `[build.rebin.phi]` / `[build.rebin.mt]` switches,
    factor/range modes, legacy configuration behavior, and pre-projection
    sparse-axis grouping
  - retained optional `[build].split_mixed_event_by_qn`, its dependency on
    SAME qn splitting, and qn-all semantics
  - retained all ten `[fit.parameters.<name>]` controls and fixed
    lambda/alpha behavior in chi2 and PML fits
  - made the ME projection compose the current rebin phi interval and qn
    selection, including the integrated-phi and qn-all cases
  - persisted qn denominator policy together with mT/phi rebin metadata across
    `SliceCatalog`, `FitCatalog`, and TSV output, with legacy defaults
  - kept `phi_all` on the full native phi span and used output indices in
    rebin-aware `group_id` / `slice_id` paths to avoid collisions
  - reconciled the Wenya config and restored
    `config/oo_build_and_fit.toml` as the stable no-argument runner default
  - updated README, the formula workflow document, tests, and the project-state
    ledger for the combined contract
  - ran `cmake --build build -j4` and the full
    `ctest --test-dir build --output-on-failure` through the O2Physics ROOT
    executor; all six registered tests passed in 32.06 seconds with
    `PRIMARY_OK`

## Recommended Owner Action

- inspect the current real scout estimate without creating output:

```bash
scripts/run_exp_femto_3d_PROFILE.sh --tier scout --profile-estimate-only
```

- then run the scout with `workers=2`; keep its output diagnostic-only
- benchmark workers 1/2/4/6 before raising concurrency, watching RSS and swap
- promote only anomalous slices to focused-1D/focused-2D and finally strict
- do not enable the Minuit2 thread backend until its A/B gate is completed

## Previous Recommended Owner Action

- use the new `[build.rebin.mt]` / `[build.rebin.phi]` contract in production
  configs as needed:

```bash
bin/exp_femto_3d build-cf --config <config.toml>
bin/exp_femto_3d fit --config <config.toml>
```

- the remaining open validation item is a real-data native/factor/ranges
  comparison against the legacy workflow; the current acceptance evidence uses
  synthetic ROOT sparses plus the existing smoke suite

- when future work changes sparse axes, CF normalization, phi mapping, rebin
  mode semantics, `SliceCatalog`, Levy/Coulomb formulas, fit metadata, or
  summary-output semantics, update `docs/数学物理公式流程说明.md` in the same
  pass

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
