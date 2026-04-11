# Current Status

## Task Snapshot

- scope: diagnose ROOT invocation failures seen during agent-driven validation,
  complete the previously skipped ROOT-dependent tests, and bootstrap durable
  project-state documentation
- current conclusion: the observed ROOT failure mode was caused by sandboxed
  `alienv` environment entry, not by `Exp_femto_3d` project logic
- primary evidence:
  - sandboxed `alienv` failed with `/dev/fd/... Operation not permitted`
  - the same ROOT checks and project tests pass in a clean non-sandboxed
    O2Physics environment

## Verification Status

- verification_status: partially verified
- project_state_sync_status: written

Reason:

- authoritative local test execution for the current unit/integration suite has
  passed in a clean O2Physics environment
- full real-data equivalence validation against the legacy macro has not yet
  been rerun during this documentation task

## Active Constraints

- ROOT-dependent validation must be run from a fully entered O2Physics
  environment
- sandboxed tool runs that fail during `alienv` bootstrap are not reliable
  evidence for code-level ROOT regressions

## Ledger Convention

- `Exp_femto_3d` adopts `project-state/` as its coordination ledger path
- this project does not use the hidden `.project-state/` path

## Why Project-State Was Missing Previously

- `Exp_femto_3d` did not yet adopt a `project-state/` directory
- the previous closeout therefore had no adopted project-state ledger to sync
- this task explicitly requested project-state updates, so the ledger is now
  bootstrapped

## Coordination Ledger State

- `project-state/` is now the active adopted coordination ledger for
  `Exp_femto_3d`
- this task completed the required parent-owned writeback after collecting
  validation evidence
- `project-state/guide.md` has now been added to satisfy the current bootstrap
  requirement for a human-facing Chinese overview
