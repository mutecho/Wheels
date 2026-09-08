# Tests

## T-023: Complete Strict-Parallel Pair Coverage

- date: 2026-09-07
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d --command 'cmake --build build --target config_parse_validation_test -j 2 && ctest --test-dir build -R config_parse_validation_test --output-on-failure'
```

- result: targeted build and CTest passed 1/1 in 1.25 s; executor `PRIMARY_OK`
  after retrying the identical command with escalation following a sandbox
  environment-entry failure
- coverage: all six unordered 2D pairs, exact parameter axis order, `21 x 21`
  grids, disabled 2D refinement, inherited bounds, unchanged four refined 1D
  scans, and eight process workers in the real strict-parallel TOML
- scope: config-only follow-up; no real OO profile job or full-suite rerun.
  The earlier full-suite evidence remains recorded below for the implementation.

## T-022: Two-Dimensional Profile And Slice Display Contract

- date: 2026-09-07
- command: complete CMake configure/build and `ctest --test-dir build
  --output-on-failure` through the shared O2Physics ROOT executor
- result: passed 8/8 in 47.52 s; executor status `PRIMARY_OK`
- coverage:
  - ROOT-independent midpoint/clipped edges, true-coordinate contours,
    invalid-corner holes, disconnected regions, constant fields, saddle and
    vertex-equality handling, duplicate suppression, and stable best tie-breaks
  - reopened TH2 ranges/titles/stats, exact `ix/iy` mapping, retained NaN bins,
    independently valid slice bins, four canvases, discrete status legend, and
    conditional best markers
  - sparse nonconverged profile points beside valid fixed-nuisance slice points,
    a finite slice surface below all 1/2/4 thresholds, an all-invalid 2D scan,
    and an unavailable common reference with no fabricated markers or colors
  - strict `write_likelihood_slice=false` object suppression, matching v3
    resume, and rejection of an injected display-contract-v2 checkpoint sidecar
- visual evidence: six canvases exported from closed/reopened toy ROOT files
  under `/tmp/exp_femto_3d_profile_2d_qa/`; no real OO profile job was run

## T-021: OO Radius Bounds, 1D Display Range, And Checkpoint Contract

- date: 2026-09-07
- command:
  `ctest --test-dir build --output-on-failure` through the shared O2Physics ROOT
  executor after a complete `cmake --build build -j 2`
- result: passed 7/7 in 41.11 s; executor status `PRIMARY_OK`
- coverage:
  - parses all four OO profile tiers and checks diagonal-radius hard bounds
    `[0.01,64.0]`, default lambda behavior, scout `[0.01,20.0]` subranges, and
    exact `_r2max64` output/checkpoint identifiers
  - rejects an explicit scan upper bound of `65.0` under the new hard domain;
    an unmodified legacy fixture still reports `[0.01,400.0]`
  - reads back `ProfileLikelihoodCatalog`, `ProfileParameterCatalog`, every
    present 1D graph, all named nuisance trajectories, and the canvas frame to
    verify one resolved X range
  - covers a fully PSD-invalid 1D scan with an empty `Profile1D`, retained
    `ProfilePoints`/`AttemptPoints` failure causes, no synthetic points, an
    explicit no-valid-points annotation, and an unchanged out-of-range nominal
    coordinate
  - changes only an inherited radius upper bound and confirms old checkpoint
    rejection plus preservation of the last complete profile ROOT output
  - preserves existing 2D grid-center, PSD-invalid, and `profile <= slice`
    regression checks without asserting 2D display semantics
- scope note: no real OO profile job was run; production rerun is operator-owned

## T-020: All-Selection Strict Parallel Profile-Only

- date: 2026-09-05
- environment: shared O2Physics ROOT executor, primary module
- result: passed, executor status `PRIMARY_OK`
- evidence:
  - configured build succeeded and complete CTest passed 7/7 in 41.56 seconds
  - config validation accepts process profile-only `fit_selection` without
    `slice_ids`; invalid process/alongside and fit-selection-plus-IDs contracts
    remain rejected
  - two-group ROOT smoke selected all 10 toy slices, produced two disjoint
    five-slice chunks, merged exactly 10 slice/scan catalog rows, and recorded
    parent scope plus 2/2 configured/effective workers
  - serial and two-process outputs matched coordinates, attempt order, winner,
    status, reference source/objective, raw objective, and delta
  - exact resume succeeded; changed `fit_selection`, missing `ProfilePoints`,
    and missing `AttemptPoints` each rejected old chunks while preserving the
    last complete final ROOT
  - production fit ROOT, report ROOT, and TSV sentinels were unchanged
  - ROOT inspector returned `STATUS: OK`; merged output contained 10
    `ProfilePoints`, 10 `AttemptPoints`, compact named nuisance graphs, and one
    `ProfileExecution` tree
  - real `strict-parallel --profile-estimate-only` returned 84 slices, 12
    groups, 10/10 workers, 153,468 maximum attempts, and 51,156 likelihood-slice
    evaluations; the profile output state remained absent
- scope note: no real strict profile calculation or physics attribution was run

## T-019: Compact Profile ROOT Display Contract

- date: 2026-09-05
- environment: shared O2Physics ROOT executor, primary module
- result: passed, executor status `PRIMARY_OK`
- evidence:
  - build succeeded and complete CTest passed 7/7 in 38.04 seconds
  - ROOT smoke confirms named `Nuisance_<parameter>` trajectories remain
  - duplicate `Nuisance_p<N>` aliases and redundant QA objects are absent
  - process resume test injects a legacy alias into a checkpoint and confirms
    the final merged ROOT prunes it while retaining the numerical trees
  - `git diff --check` passed

## T-018: Profile Acceleration, Process Resume, And ROOT Closeout

- date: 2026-09-04
- environment: shared O2Physics ROOT executor, primary module
- result: passed, executor status `PRIMARY_OK`
- evidence:
  - configured build succeeded
  - final complete CTest passed 7/7 in 35.49 seconds
  - `workflow_smoke_test` passed with profile-only production sentinels,
    `posix_spawn` chunk generation/merge, checkpoint reuse, mismatch rejection,
    timing/HESSE schema checks, PSD-invalid classification, and
    `profile <= slice`
  - ordered cache evaluator passed the built-in old/new objective tolerance
  - ROOT inspector reported `STATUS: OK` for the merged process output: five
    trees, thirty graphs, two canvases, and the expected profile directories
  - real scout `--profile-estimate-only` returned 14 slices, 1,386 maximum
    attempts, zero slice evaluations, and did not create an output
  - `bash -n scripts/run_exp_femto_3d_PROFILE.sh` and `git diff --check` passed
- scope note: no expensive real OO scan or worker-scaling benchmark was run

## T-017: PML Profile-Likelihood Diagnostic Mode

- date: 2026-09-03
- intended environment: O2Physics ROOT executor with the configured CATS build
- ROOT-independent command:

```bash
clang++ -std=c++17 -Wall -Wextra -Wpedantic \
  -Iinclude -Isrc \
  tests/ProfileLikelihoodDriverTest.cpp src/ProfileLikelihood.cpp \
  -o /tmp/exp_femto_3d_profile_driver_test
/tmp/exp_femto_3d_profile_driver_test
python3 -c 'import tomllib; tomllib.load(open("config/examples/exp_femto_3d.example.toml", "rb"))'
git diff --check
```

- result: passed
- covered evidence:
  - deterministic 1D/2D grids fix the intended target coordinates
  - nuisance re-minimization produces a true profile and satisfies the scripted
    `profile <= slice` case
  - nominal, forward-neighbor, and reverse-neighbor attempts are retained and
    the lowest valid branch wins
  - non-converged/error/domain/objective failures remain represented
  - interior coarse minima trigger one regular refinement pass
  - a lower valid profile point triggers per-slice reference recomputation
  - source inspection confirms the pure driver has no ROOT dependency
  - the updated public example remains valid TOML
- additional coverage written but not independently executed at closeout:
  - config tests for default-off behavior, PML prerequisite,
    listed/fit-selection scopes, multiple scans, safe IDs, dimensions, ranges,
    point counts, contours, and refinement pairing
  - separate-file schema and numerical tree fields
  - exact 1D/2D coordinates and attempt rows
  - likelihood slice comparison, nuisance output, failure/status mask,
    nominal/profile/boundary markers, and output-collision ordering
  - an unreachable full-model off-diagonal scan is expected to produce
    `model_domain_invalid` rows; the compact display contract does not duplicate
    that category as a standalone graph
  - CLI `--model` target revalidation and fixed/unknown target rejection
- blocked command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command 'cmake --build build -j4 && ctest --test-dir build --output-on-failure'
```

- blocked result:
  - sandbox execution returned `STATUS: ESCALATION_REQUIRED` with
    `/dev/fd/... Operation not permitted`
  - after explicit user authorization, the escalated call was rejected by the
    approval backend with HTTP 403
  - no alternate ROOT installation was used, and no real OO data was run
- significance:
  implementation and ROOT-independent scan semantics are verified; ROOT ABI,
  file schema execution, and canvas visual QA remain an explicit open item

## T-016: Resolved Merge Compatibility Matrix

- date: 2026-08-13
- environment: O2Physics ROOT executor, CATS finite-source support enabled
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'cmake --build build -j4 && ctest --test-dir build --output-on-failure'
```

- result: passed, executor status `PRIMARY_OK`
- evidence:
  - all affected production and test targets compiled successfully
  - full CTest passed 6/6 in 32.06 seconds
  - `config_parse_validation_test` covers rebin modes, qn-ME dependencies,
    and all supported Levy parameter override contracts
  - `slice_catalog_roundtrip_test` covers qn denominator metadata, rebin
    metadata, and legacy defaults
  - `build_cf_rebin_test` covers combined qn/phi/mT selections, SE/ME
    closure, integral behavior, q-axis preservation, paths, physical ranges,
    seam/alignment/divisibility failures, and safe output ordering
  - `workflow_smoke_test` covers FitCatalog propagation and fixed
    lambda/alpha behavior in chi2 and PML
- significance:
  validates the resolved semantic union rather than only either pre-merge side

## T-015: Optional ME Qn Split And Levy Parameter Overrides

- date: 2026-07-01
- result: passed in the original feature validation; superseded for merge
  confidence by T-016
- evidence:
  - config validation covered the qn split dependency and all ten fit parameter
    tables
  - ROOT-backed catalog/workflow checks covered qn-specific ME narrowing and
    fixed lambda/alpha in chi2 and PML

## T-014: Configurable Build-Side mT/phi Rebin Contract

- date: 2026-08-12
- result: passed in the original feature validation; superseded for merge
  confidence by T-016
- evidence:
  - O2Physics build and full six-test CTest returned `PRIMARY_OK`
  - synthetic ROOT coverage checked factor/ranges modes, CF closure, metadata,
    unique paths, q-axis preservation, invalid-boundary failures, and legacy
    catalog behavior


## T-013: Wenya PbPb LHC23 qn-All Plus qn1/qn2/qn3 Production Build/Fit

- date: 2026-06-30
- environment: O2Physics ROOT executor, `PRIMARY_OK` for ROOT-backed commands
- command:

```bash
./scripts/cmake.sh
ctest --test-dir build --output-on-failure

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/slice_catalog_roundtrip_test'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/workflow_smoke_test'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/exp_femto_3d build-cf --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/exp_femto_3d fit --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml'
```

- result: passed
- evidence:
  - `scripts/cmake.sh` configured with CATS finite-source support enabled and
    rebuilt all targets successfully
  - `ctest --test-dir build --output-on-failure` reported `100% tests passed`;
    the two ROOT-backed tests were skipped by the ctest environment rule, then
    were run directly through the O2Physics executor
  - `bin/slice_catalog_roundtrip_test` returned `PRIMARY_OK` and covered qn
    metadata defaults, legacy catalog compatibility, and qn-all/qn1/qn2/qn3
    toy build output
  - `bin/workflow_smoke_test` returned `PRIMARY_OK` and covered qn metadata in
    `FitCatalog` and `CoulombKernelCatalog`
  - production `build-cf` returned `PRIMARY_OK` and printed
    `stored_slices=1040`, `skipped_zero_me_groups=0`,
    `skipped_zero_mixed_event_slices=0`, and `skipped_zero_se_slices=0`
  - production `fit` returned `PRIMARY_OK` and printed `fitted_slices=468`,
    `selected_slices=468`, `missing_objects=0`, and
    `missing_raw_histograms=0`
  - ROOT inspection returned `PRIMARY_OK`: `meta/SliceCatalog` has `1040`
    entries; qn_all/qn1/qn2/qn3 each have `260` build slices and `20`
    phi-all build slices; qn-all paths retain the legacy no-qn suffix
  - detailed fit `meta/FitCatalog` has `468` entries; qn_all/qn1/qn2/qn3 each
    have `117` fitted slices and `9` phi-all fitted slices; representative
    qn-all and qn1/qn2/qn3 `summary/R2_vs_phi` objects exist
  - report ROOT contains representative qn-all and qn1 `summary/R2_vs_phi`,
    `source_parameters`, and `eps_vs_mt` objects
  - `fit_summary_wenya_lhc23_merge_qn_split_plus_integrated.tsv` has `469`
    lines, matching one header plus `468` fitted rows
- significance: validates that the Wenya LHC23 config now preserves the
  original qn-integrated result while adding reference-style qn1/qn2/qn3 CF
  and fit outputs in the standard `Exp_femto_3d` structure

## T-012: Wenya PbPb LHC23 qn-Integrated Production Build/Fit

- date: 2026-06-30
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bash scripts/cmake.sh'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'ctest --test-dir build --output-on-failure'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/exp_femto_3d build-cf --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'bin/exp_femto_3d fit --config config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml'
```

- result: passed
- evidence:
  - `scripts/cmake.sh` returned `PRIMARY_OK`; the active build kept CATS support
    enabled and rebuilt `config_parse_validation_test`
  - `ctest --test-dir build --output-on-failure` returned `PRIMARY_OK` and
    passed all five registered tests:
    `coulomb_kernel_validation_test`, `config_parse_validation_test`,
    `progress_render_test`, `slice_catalog_roundtrip_test`, and
    `workflow_smoke_test`
  - `config_parse_validation_test` now parses
    `config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml` and checks its
    target input path, PbPb-style bin selections, qn-integrated output names,
    raw phi mapping, integrated-ME denominator, full fit model, Gamow legacy
    Coulomb setting, q2 baseline, PML, and follow-input fit phi metadata
  - production `build-cf` returned `PRIMARY_OK` and printed
    `stored_slices=260`, `skipped_zero_me_groups=0`,
    `skipped_zero_mixed_event_slices=0`, and `skipped_zero_se_slices=0`
  - production `fit` returned `PRIMARY_OK` and printed `fitted_slices=117`,
    `selected_slices=117`, `missing_objects=0`, and
    `missing_raw_histograms=0`
  - ROOT inspection returned `PRIMARY_OK`: `meta/SliceCatalog` has `260`
    entries, no `qn_index` branch, no `slice_id` containing `qn`, and `20`
    phi-integrated slices; detailed fit `meta/FitCatalog` has `117` entries
  - detailed and report ROOT outputs contain representative `summary/R2_vs_phi`
    objects, and report ROOT contains `meta/FitCatalog` plus a representative
    `source_parameters` overview canvas
  - `fit_summary_wenya_lhc23_merge_qn_integrated.tsv` has `118` lines,
    matching one header plus `117` fitted rows
- significance: validates that the Wenya LHC23 input can be processed through
  the original `Exp_femto_3d` catalog/slice structure with qn integrated, not
  through the flat reference-macro object layout

## T-011: Incremental Build Helper Refreshes ROOT Links Only When Needed

- date: 2026-06-23
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash -n /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/scripts/cmake.sh

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command 'bash scripts/cmake.sh'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command 'bash scripts/cmake.sh'

rg -n "ROOT_DIR|local8|local7" build/CMakeCache.txt build/CMakeFiles/exp_femto_3d.dir/link.txt
otool -l bin/exp_femto_3d | rg -A2 "LC_RPATH|path .*ROOT|path .*CATS"
otool -L bin/exp_femto_3d | rg "CATS|gsl"

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command 'ctest --test-dir build --output-on-failure'
```

- result: passed
- evidence:
  - `bash -n scripts/cmake.sh` passed
  - first `scripts/cmake.sh` run detected cached
    `ROOT/v6-36-10-alice1-local7` and refreshed to active
    `ROOT/v6-36-10-alice1-local8`
  - second `scripts/cmake.sh` run returned `PRIMARY_OK` with only
    `Built target ...` lines, confirming no compile or link work when inputs
    are unchanged
  - `build/CMakeFiles/exp_femto_3d.dir/link.txt` and the binary `LC_RPATH`
    both point at `ROOT/v6-36-10-alice1-local8/lib`
  - `otool -L bin/exp_femto_3d` showed `libCATS.dylib`, `libgsl`, and
    `libgslcblas`
  - `ctest --test-dir build --output-on-failure` returned `PRIMARY_OK` and
    passed all five registered tests
- significance: keeps stale ROOT/CATS link protection while avoiding the
  previous unconditional clean rebuild on every helper invocation

## T-010: Clean-First Build Helper Rebuilds CATS Operator Binary

- date: 2026-06-22
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash -n /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/cmake.sh

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command '/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/cmake.sh'

otool -L /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/bin/exp_femto_3d | grep -E 'CATS|gsl'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'ctest --test-dir /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build --output-on-failure'

git diff --check
```

- result: passed
- evidence:
  - `bash -n scripts/cmake.sh` passed
  - `scripts/cmake.sh` returned `PRIMARY_OK`, configured with
    `CATS finite-source Coulomb support enabled`, and relinked
    `bin/exp_femto_3d`
  - `otool -L bin/exp_femto_3d` showed `libCATS.dylib`, `libgsl`, and
    `libgslcblas`
  - `ctest --test-dir build --output-on-failure` returned `PRIMARY_OK` and
    passed all five registered tests:
    `coulomb_kernel_validation_test`, `config_parse_validation_test`,
    `progress_render_test`, `slice_catalog_roundtrip_test`, and
    `workflow_smoke_test`
  - `git diff --check` returned exit code 0
- significance: verifies the build helper no longer trusts a stale source-tree
  `bin/exp_femto_3d` that may have been written by another build tree

## T-009: Finite-Source Coulomb Mode And Kernel Validation

- date: 2026-06-21
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_3d -B /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build && cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build -j4 && ctest --test-dir /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build --output-on-failure'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_3d -B /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build_no_cats -DEXP_FEMTO_3D_ENABLE_CATS=OFF && cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build_no_cats -j4 && ctest --test-dir /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build_no_cats --output-on-failure'

git diff --check
```

- result: passed
- evidence:
  - default configure reported `CATS finite-source Coulomb support enabled`
    with `/Users/allenzhou/Research_software/CATS/install/lib/libCATS.dylib`
  - CATS-enabled `ctest` passed 5/5 tests:
    `coulomb_kernel_validation_test`, `config_parse_validation_test`,
    `progress_render_test`, `slice_catalog_roundtrip_test`, and
    `workflow_smoke_test`
  - no-CATS configure/build/`ctest` also passed 5/5 tests and exercised the
    explicit finite-source unavailable failure path
  - `coulomb_kernel_validation_test` covers Coulomb mode predicates,
    `q_inv -> k*` conversion, finite-source table low-k clamping,
    interpolation, high-k unity fallback, no-FSI CATS unity reference, and
    source-radius ordering against point-source Gamow behavior
  - `workflow_smoke_test` covers finite-source `fixed_1d` and `iterative_1d`
    toy workflows plus detailed/report `meta/CoulombKernelCatalog` output
  - `git diff --check` returned exit code 0
- significance: verifies the plan's public mode contract, CATS integration
  boundary, kernel behavior, and ROOT output metadata without requiring a
  production real-data write
- limitation: production OO/PbPb real-data regression has not yet been rerun
  after this implementation

## T-008: Eventgen-Style Progress Rendering

- date: 2026-06-10
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_3d -B /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
  --command 'cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build -j4'

bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build \
  --command 'ctest --output-on-failure'
```

- result: passed
- evidence:
  - CMake configure found `Threads`
  - build completed all targets, including `progress_render_test`
  - `ctest` passed 4/4 tests:
    `config_parse_validation_test`, `progress_render_test`,
    `slice_catalog_roundtrip_test`, and `workflow_smoke_test`
- significance: verifies the new progress formatter prints percent/activity/ETA
  fields and that the thread-linked progress reporter still builds and coexists
  with ROOT-backed workflow tests

## T-007: No-Argument Runtime Re-Entry Nounset Check

- date: 2026-06-10
- environment: local shell
- command:

```bash
bash -n /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh
/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh --help
bash -uc 'original_arg_count=$#; original_args=("$@"); if (( original_arg_count > 0 )); then printf "arg:%s\n" "${original_args[@]}"; else printf "no-args\n"; fi' _
bash -uc 'original_arg_count=$#; original_args=("$@"); if (( original_arg_count > 0 )); then printf "arg:%s\n" "${original_args[@]}"; else printf "no-args\n"; fi' _ --stage fit
```

- result: passed
- evidence:
  - shell syntax validation returned exit code 0
  - help output returned exit code 0
  - empty-argument branch printed `no-args` without triggering `set -u`
  - non-empty preserved-argument branch printed `arg:--stage` and `arg:fit`
- significance: verifies the no-argument wrapper path no longer expands an
  empty argument array while preserving argument forwarding for explicit stage
  runs
- limitation: did not run the default `all` stage because
  `oo_build_and_fit.toml` writes configured production outputs

## T-006: Run Script Syntax And Help Check

- date: 2026-06-10
- environment: local shell
- command:

```bash
bash -n /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh
/Users/allenzhou/Research_software/Code_base/Exp_femto_3d/scripts/run_exp_femto_3d.sh --help
```

- result: passed
- evidence:
  - shell syntax validation returned exit code 0
  - help output showed the default `config/oo_build_and_fit.toml` run contract,
    `--stage all|build-cf|fit`, `--model full|diag`, `--input-cf-root`, and
    `--binary`
- significance: verifies the new operator wrapper can be parsed and documents
  the intended CLI surface before any ROOT-backed real-data run
- limitation: did not run the default `all` stage because
  `oo_build_and_fit.toml` writes configured production outputs

## T-005: O2Physics CTest Run After Standalone Fit Report Output

- date: 2026-06-10
- environment: O2Physics ROOT executor, `PRIMARY_OK`
- command:

```bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build \
  --command 'ctest --output-on-failure'
```

- result: passed
- evidence:
  - `config_parse_validation_test`: passed
  - `slice_catalog_roundtrip_test`: passed
  - `workflow_smoke_test`: passed
- significance: verifies parsing for `output.fit_report_directory` and
  `output.fit_report_root_name`, backward-compatible report-directory defaulting,
  path-collision rejection for report files that would overwrite existing ROOT
  artifacts,
  and ROOT smoke coverage for `meta/FitCatalog`, mirrored `summary/R2_vs_phi`,
  `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`, and
  `eps_vs_mt/<cent>/epsf_vs_mt(_canvas)` in the standalone report file

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

## T-004: O2Physics CTest Run After ME Phi-Splitting Switch

- date: 2026-05-18
- environment: O2Physics ROOT executor
- command:

```bash
ctest --output-on-failure
```

- result: passed
- evidence:
  - `config_parse_validation_test`: passed
  - `slice_catalog_roundtrip_test`: passed
  - `workflow_smoke_test`: passed
- significance: verifies `build.split_mixed_event_by_phi` parsing, catalog
  metadata persistence, legacy catalog default `false`, and phi-dependent
  `ME_raw3d` behavior when split-ME mode is enabled
