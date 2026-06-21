# Tests

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
