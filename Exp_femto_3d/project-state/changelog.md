# Changelog

## 2026-04-11

- Bootstrapped `project-state/` for `Exp_femto_3d`.
- Adopted `project-state/` as the durable coordination ledger path for
  `Exp_femto_3d`.
- Confirmed that prior ROOT/THnSparse failures during agent execution were
  caused by sandboxed `alienv` bootstrap failure, not by project logic.
- Completed the previously incomplete ROOT-dependent tests in a clean
  non-sandboxed O2Physics environment.
- Added `ROOT_RUNTIME_AGENT_NOTE.md` as a durable diagnostic reference for
  future agents.
- Synced the ledger again to make the non-hidden `project-state/` convention
  explicit and record `project_state_sync_status`.
- Added the required Chinese `project-state/guide.md` after confirming the
  current skill explicitly requires it for bootstrap.
