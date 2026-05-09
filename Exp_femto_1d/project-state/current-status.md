# Current Status

## Task Snapshot

- scope: execute the `Exp_femto_1d/PLAN.md` engineering refactor and establish
  `project-state/`
- current conclusion: the project has been scaffolded as an independent CMake
  package with public config/logging/workflow/model interfaces, structured CF
  output, per-centrality CF-by-mT overlay canvases, fit workflow, example configs,
  test entry points, and a project-local config-driven run helper under `scripts/`

## Verification Status

- verification_status: locally verified
- project_state_sync_status: written

Reason:

- config, workflow, and CATS smoke tests are implemented in-tree
- `2026-04-20` non-sandboxed O2Physics configure/build/`ctest --output-on-failure`
  passed all registered tests
- `2026-05-06` script-level validation checked
  `scripts/run_exp_femto_1d.sh` syntax, help output, and stage dispatch
- `2026-05-06` direct script execution without a pre-entered ROOT shell was
  fixed to re-enter `O2Physics/latest-master-o2` before launching the binary;
  `--stage build-cf` completed with 75 stored slices
- `2026-05-06` internal `THnSparse` projection ownership was fixed so repeated
  EP interval projections no longer emit `TROOT::Append` replacement warnings;
  O2Physics `ctest --output-on-failure` passed all 4 registered tests
- `2026-05-08` `build-cf` started writing
  `cent_slices/<cent_id>/<region_name>/CFByMtCanvas`; local build passed and
  O2Physics `ctest --output-on-failure` passed all 4 registered tests
- `2026-05-08` `CFByMtCanvas` markers became opt-in through
  `build.cf_by_mt_show_markers`; default output is line-only for clearer trends
- real-data equivalence to the legacy macro has not yet been run

## Active Constraints

- ROOT-backed validation requires a full O2Physics environment; the run helper
  now enters the default `O2Physics/latest-master-o2` module automatically when
  `alienv` is available
- direct helper execution defaults to `build-cf`; full build-and-fit execution
  requires explicit `--stage all`
- the run helper still requires a built `bin/exp_femto_1d` executable
- real-data closure is still pending for both `build-cf` and `fit`

## Ledger Convention

- `Exp_femto_1d` adopts `project-state/` as the repository-local coordination ledger
- this project does not use `.project-state/`
