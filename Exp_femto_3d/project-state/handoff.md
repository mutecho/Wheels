# Handoff

## Latest Durable Handoff

- configuration follow-up completed on 2026-09-07:
  - strict-parallel now has four 1D scans plus all six 2D pairs of
    `lambda`, `rout2`, `rside2`, and `rlong2`; each 2D scan uses an
    unrefined `21 x 21` grid with inherited hard bounds
  - removed the stale README statement that 2D display work remained deferred;
    synchronized runner help and the configuration coverage assertions
  - targeted config build/CTest passed 1/1 with O2Physics `PRIMARY_OK`;
    no real profile job or full-suite rerun was performed for this follow-up
  - adding scans changes the existing digest: earlier five-scan chunks are
    incompatible even if their display contract is already v3; use a fresh
    run ID if that old checkpoint generation exists

- completed on 2026-09-07 (issue 2):
  - added ROOT-independent `ProfileDisplay2D` geometry for clipped TH2 edges,
    stable coarse/refined best selection, and four-corner-valid marching-squares
  - replaced automatic NaN-contaminated contours with disconnected explicit
    segments and added gray invalid masks plus a fixed six-category status panel
  - added full-range profile canvases and independent 2D fixed-nuisance slice
    matrices/canvases, with persisted reference/count/boundary/contour annotations
  - advanced checkpoint compatibility to `display-contract-v3`; matching v3
    resume passes and v2 sidecars are rejected
  - O2Physics ROOT `PRIMARY_OK`; full CTest passed 8/8 in 47.52 s and six
    closed/reopened toy canvases passed visual QA
  - final code, physics-contract, and plot reviews found no unresolved blocker
- operator-owned next steps:
  - manually run the intended `_r2max64` OO profile configuration; do not reuse
    a same-name v2 chunk (choose a new `_displayv3` run ID if needed)
  - assess new-range real-data convergence separately; toy QA validates display
    semantics, not production minimizer quality

- completed on 2026-09-07:
  - set the four OO profile tiers to hard diagonal-radius bounds
    `[0.01,64.0] fm^2`, while preserving scout's explicit `[0.01,20.0]`
    diagnostic scan subrange and the code-wide legacy `[0.01,400]` defaults
  - versioned all four profile ROOT names and checkpoint `run_id`s with
    `_r2max64`
  - unified every persisted 1D graph and `Canvas_1D` on the resolved scan
    range; empty valid-point sets now retain the axis and show an explicit
    diagnostic note without fabricated points
  - upgraded checkpoint reuse to hash effective parameter domains/fixed state,
    fit/minimizer settings, and resolved scan ranges
- verified:
  - O2Physics ROOT executor build and full CTest passed 7/7 in 41.11 s with
    `PRIMARY_OK`
  - ROOT toy coverage reads catalogs, graphs, canvas frames, nuisance
    trajectories, invalid-point trees, and inherited-bound checkpoint mismatch
  - `git diff --check` passed
- intentionally pending:
  - the operator will run the real OO profile jobs manually
  - convergence quality in the new `[0.01,64]` production range remains to be
    assessed from that user-run output

## Previous Durable Handoff

- completed on 2026-09-05:
  - enabled `profile_only + process` for parent-side materialized
    `slice_scope="fit_selection"`; child workers are narrowed to their exact
    assigned group and cannot expand the source TOML selection again
  - added strict all-selection config
    `config/oo_build_and_fit_6bins_profile_strict_parallel.toml` and runner tier
    `scripts/run_exp_femto_3d_PROFILE.sh --tier strict-parallel`
  - fixed common child numerical-library thread counts to 1 and added original
    scope, selected slice/group counts, and configured/effective workers to
    `meta/ProfileExecution`
  - strengthened chunk validation so every expected scan must contain readable,
    non-empty `ProfilePoints` and `AttemptPoints`
- verified:
  - O2Physics executor build and full CTest passed 7/7 in 41.56 s with
    `PRIMARY_OK`
  - two-group toy covers complete fit-selection expansion, one-group-only child
    chunks, exact final catalog cardinality, matching resume, changed-selection
    digest rejection, both missing-tree cases, serial/process numerical
    equivalence, and production sentinels
  - ROOT inspector returned `PRIMARY_OK` / `STATUS: OK` for the 10-slice merged
    process output
  - real strict-parallel estimate-only passed: 84 slices, 12 groups, 10/10
    configured/effective workers, 153,468 maximum attempts, and no output file
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

## Superseded Issue-1 Owner Action

- rerun the desired OO profile tier manually; the new output/checkpoint
  generation is identified by `_r2max64`, so the prior generation is preserved
- interpret a best-fit radius at `64 fm^2` only as an optimum on the imposed
  constraint boundary, not as a localized unconstrained minimum
- the former instruction to defer 2D interpretation was completed by the
  issue-2 handoff above; use the v3 canvases and annotations for new outputs

## Previous Recommended Owner Action

- inspect the all-selection strict contract without creating output:

```bash
scripts/run_exp_femto_3d_PROFILE.sh --tier strict-parallel --profile-estimate-only
```

- the strict-parallel estimate is very expensive; normally run scout/focused
  first and reserve strict-parallel for selected operational needs
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
