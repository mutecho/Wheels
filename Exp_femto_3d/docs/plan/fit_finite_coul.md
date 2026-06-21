# Exp_femto_3d Finite-Source Coulomb Kernel Implementation Plan

## Context

Repo: `/Users/allenzhou/Research_software/Code_base/Exp_femto_3d`.

Current state:
- Fit Coulomb is a boolean `use_coulomb`.
- `use_coulomb=true` uses point-source Gamow in `Workflow.cpp`.
- Fit options flow through `LevyFitOptions`.
- `fit` writes detailed ROOT, TSV, and standalone report ROOT.
- Current ROOT validation must run through:
  ```bash
  bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
    --cwd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
    --command '...'
  ```
- Local CATS exists at `/Users/allenzhou/Research_software/CATS`, with headers/libs under `install/`; examples for ππ Gaussian source are in `work/fit_pipi_test.cpp` and `work/theory_pipi.cpp`.

Goal:
Add Coulomb fit modes `none`, `gamow`, `finite_source`, with finite-source supporting `fixed_1d` and `iterative_1d`. This round only implements 1D spherical finite-source Coulomb; 3D anisotropic kernel is explicitly out of scope.

## Public Interface

Add config fields under `[fit]`:

```toml
coulomb_mode = "none"            # none | gamow | finite_source
finite_source_mode = "fixed_1d"  # fixed_1d | iterative_1d, only used by finite_source
```

Backward compatibility:
- If `coulomb_mode` is absent:
  - `use_coulomb = false` maps to `none`.
  - `use_coulomb = true` maps to `gamow`.
- If both are present and agree, accept.
- If both are present and conflict, throw `ConfigError`.
- Keep `use_coulomb` accepted in existing configs and tests.

Add output metadata:
- TSV and `meta/FitCatalog` keep existing `usesCoulomb`.
- Add fields:
  - `coulombMode`: string, `none|gamow|finite_source`
  - `finiteSourceMode`: string, empty or `fixed_1d|iterative_1d`
  - `finiteSourceRadiusFm`: double, NaN unless finite-source was used
- Add `meta/CoulombKernelCatalog` to fit ROOT and report ROOT with per cent/mT group:
  - `group_id`, centrality/mT bounds
  - `finiteSourceMode`
  - seed `R_eff`
  - final `R_eff`
  - CATS enabled flag
  - kernel momentum range and bin count

## Implementation Changes

Types/config:
- In `Types.h`, replace `bool use_coulomb` in `LevyFitOptions` with:
  ```cpp
  enum class CoulombMode { kNone, kGamow, kFiniteSource };
  enum class FiniteSourceMode { kFixed1D, kIterative1D };
  ```
  Keep a helper like `UsesCoulomb(mode)` for output compatibility.
- Add `ToString/Parse` helpers in `Config.cpp` and declarations in `Config.h`.
- Validate:
  - `finite_source_mode` is only meaningful with `coulomb_mode="finite_source"`.
  - finite-source requires CATS support at runtime/build-time; otherwise fail clearly before starting fit.

CMake:
- Add optional CATS detection:
  - include path: `/Users/allenzhou/Research_software/CATS/install/include`
  - lib: `/Users/allenzhou/Research_software/CATS/install/lib/libCATS.dylib`
  - GSL: Homebrew `/opt/homebrew/opt/gsl`
- If found, define `EXP_FEMTO_3D_HAS_CATS=1` and link CATS/GSL.
- If not found, build still succeeds for `none/gamow`.

Kernel layer:
- Add a small internal kernel abstraction in `Workflow.cpp` or a new focused source/header:
  - `EvaluateCoulombFactor(q_out, q_side, q_long, options, kernel_context)`
  - `Gamow` calls the existing implementation.
  - `None` returns 1.
  - `FiniteSource` looks up a CATS-generated spline/table.
- The finite-source kernel must be Coulomb-only:
  - CATS source: spherical Gaussian, same form as local CATS examples.
  - Use `SetQ1Q2(1)`, ππ reduced mass, and `SetQuantumStatistics(0)` to avoid double-counting QS.
  - The existing Bowler model keeps its own `1 + exp(-levy_exponent)` QS term.
- Momentum convention:
  - Fit q is in GeV/c.
  - `q_inv ≈ sqrt(qout² + qside² + qlong²)`.
  - CATS momentum is `k*` in MeV/c:
    ```text
    kstar_MeV = 0.5 * q_inv_GeV * 1000
    ```
  - Clamp below first table bin and above last table bin; high q should approach or use 1 if outside support.

Finite-source radius:
- For each cent/mT group:
  ```text
  R_eff = sqrt((Rout2 + Rside2 + Rlong2) / 3)
  ```
- Reject non-finite or non-positive diagonal radii with a clear error.

Fit flow:
- Split current `FitAndWriteSingleSlice` into:
  - `FitSingleSlice(...)` returns `TF3`/`LevyFitResult` without writing.
  - `WriteSingleSliceFitArtifacts(...)` writes final CF, fit function, projections, canvases.
- Seed/iteration fits must not write ROOT outputs.
- For finite-source:
  - First run `phi_all` for each cent/mT group using `gamow` to seed `R_eff`.
  - `fixed_1d`: build CATS kernel from seed `R_eff`; fit all selected slices with this kernel.
  - `iterative_1d`: build kernel from seed `R_eff`, refit `phi_all` once with finite-source, rebuild kernel from updated `R_eff`, then fit all selected slices.
- Keep final output only from the final selected mode.

Model parameters:
- Replace fixed `UseCoulomb` TF3 parameter with `CoulombModeCode`.
- Do not let Minuit vary the Coulomb mode; it remains fixed.
- Ensure PML FCN and non-PML `TH3::Fit` paths both use the same kernel evaluator.

Logging:
- At fit start, log `coulombMode` and `finiteSourceMode`.
- For finite-source, log each group’s seed/final `R_eff` and kernel source.
- If CATS unavailable and finite-source requested, fail before fitting any slice.

Docs/project-state:
- Update README config section.
- Update `project-state/current-status.md`, `project-state/tests.md`, and `project-state/handoff.md` after implementation and validation.

## Test Plan

Config tests:
- Existing `use_coulomb=false` config parses as `none`.
- Existing `use_coulomb=true` config parses as `gamow`.
- `coulomb_mode="none"|"gamow"|"finite_source"` parses.
- Invalid mode throws.
- Conflicting `use_coulomb` and `coulomb_mode` throws.
- `finite_source_mode` invalid throws.

Unit/smoke kernel tests:
- If CATS is available, generate a finite-source ππ Coulomb-only kernel for a known radius.
- Verify `SetQuantumStatistics(0)` behavior by confirming no QS-only enhancement is baked into the kernel.
- Verify finite-source small-radius trend approaches Gamow more closely than a large-radius kernel.
- Verify q to k* conversion and bounds behavior.

Workflow tests:
- Extend `workflow_smoke_test` for `none` and `gamow`.
- Add finite-source smoke only when `EXP_FEMTO_3D_HAS_CATS=1`.
- Check output ROOT contains `meta/FitCatalog`, `meta/CoulombKernelCatalog`, and final slice artifacts.
- Check TSV includes new Coulomb mode columns.

Validation commands:
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

Real-data validation:
- Do not overwrite current production files.
- Create separate output names for:
  - `none`
  - `gamow`
  - `finite_source_fixed_1d`
  - `finite_source_iterative_1d`
- Start with a limited cent/mT selection.
- Compare:
  - low-q projection pulls for X/Y/Z
  - `alpha` boundary count
  - `lambda` and `R²` φ fluctuations
  - `R²_vs_phi` secondary fit pull, not only visual error size

## Assumptions And Defaults

- CATS is optional at build time, required only for `finite_source`.
- `fixed_1d` and `iterative_1d` both seed from `phi_all` Gamow fit.
- `iterative_1d` means exactly one finite-source update round in this version; multi-round convergence is deferred.
- Finite-source kernel is spherical 1D Gaussian only.
- 3D anisotropic finite-source kernel and CRAB spline input are deferred.
- The implementation should stay local to `Exp_femto_3d`; do not modify CATS itself.
