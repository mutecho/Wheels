# Changelog

## 2026-04-19

- Synced `project-state/` to the current `3d_cf_from_exp` refactor work rather
  than leaving it at the earlier ROOT-runtime bootstrap state.
- Recorded the new build/fit progress-mode controls exposed through TOML.
- Recorded that build now persists `build_uses_symmetric_phi_range` into
  `meta/SliceCatalog`.
- Recorded that fit can either follow the input CF phi metadata or override it
  from stored `raw_phi_*` coordinates without rebuilding the CF file.
- Recorded backward-compatible inference for legacy `SliceCatalog` trees that do
  not yet carry the new branch.
- Recorded the authoritative `2026-04-19` non-sandboxed O2Physics `ctest` pass
  and the matching sandbox skip signature for context.

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
