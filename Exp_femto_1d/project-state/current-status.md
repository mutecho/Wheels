# Current Status

## Task Snapshot

- scope: execute the `Exp_femto_1d/PLAN.md` engineering refactor and establish
  `project-state/`
- current conclusion: the project has been scaffolded as an independent CMake
  package with public config/logging/workflow/model interfaces, structured CF
  output, fit workflow, example configs, and test entry points

## Verification Status

- verification_status: locally verified
- project_state_sync_status: written

Reason:

- config, workflow, and CATS smoke tests are implemented in-tree
- `2026-04-20` non-sandboxed O2Physics configure/build/`ctest --output-on-failure`
  passed all registered tests
- real-data equivalence to the legacy macro has not yet been run

## Active Constraints

- ROOT-backed validation requires a full O2Physics environment
- real-data closure is still pending for both `build-cf` and `fit`

## Ledger Convention

- `Exp_femto_1d` adopts `project-state/` as the repository-local coordination ledger
- this project does not use `.project-state/`
