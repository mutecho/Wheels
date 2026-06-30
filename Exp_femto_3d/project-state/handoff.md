# Handoff

## Latest Durable Handoff

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

## Next Recommended Owner Action

- when future work changes sparse axes, CF normalization, phi mapping,
  `SliceCatalog`, Levy/Coulomb formulas, fit metadata, or summary-output
  semantics, update `docs/数学物理公式流程说明.md` in the same pass
- no runtime validation is required for this documentation-only handoff, but
  the existing O2Physics ROOT executor validation path remains the authoritative
  route for behavior changes

## Previous Durable Handoff

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

## Previous Recommended Owner Action

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
