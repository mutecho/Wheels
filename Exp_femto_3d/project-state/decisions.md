# Decisions

## DEC-001: Treat Sandboxed Alienv Bootstrap Failure As Environment-Entry Failure

- date: 2026-04-11
- decision: when `alienv` emits `/dev/fd/... Operation not permitted` during an
  agent-run ROOT command, classify the result as environment-entry failure
  before investigating project code
- rationale:
  - the same ROOT checks pass in a clean non-sandboxed O2Physics environment
  - the same built project tests pass in that environment
  - the failure signature is therefore execution-context dependent

## DEC-002: Preserve A Durable Agent Diagnostic Note For ROOT Runtime Triage

- date: 2026-04-11
- decision: keep `ROOT_RUNTIME_AGENT_NOTE.md` in the project root as the
  canonical agent-facing explanation of this failure mode
- rationale: future agent runs need a local reference that prevents repeated
  misdiagnosis of sandboxed ROOT/O2 entry failures

## DEC-003: Use `project-state/` As The Coordination Ledger Path

- date: 2026-04-11
- decision: `Exp_femto_3d` uses `project-state/` rather than `.project-state/`
  as its adopted coordination ledger path
- rationale:
  - this matches the project's current local convention
  - future agent closeout should write to the explicit non-hidden ledger path
    instead of inferring the hidden default

## DEC-004: Keep Integrated ME As Default And Make Phi-Following ME Opt-In

- date: 2026-05-18
- decision: `build.split_mixed_event_by_phi = false` remains the default so each
  centrality/mT group uses one mixed-event denominator integrated over phi;
  `true` is an explicit opt-in for denominators that follow each SE phi slice
- rationale:
  - preserving the default keeps existing CF production numerically backward
    compatible
  - the opt-in mode makes phi-differential denominator checks possible without
    changing `map_pair_phi_to_symmetric_range`, which only controls coordinate
    mapping for display and fit summaries

## DEC-005: Keep Fit Detail Output And Fit Report Output Separate

- date: 2026-06-10
- decision: keep per-slice fit objects in `output.fit_root_name`, and write the
  summary/report presentation to a separate ROOT file controlled by
  `output.fit_report_directory` and `output.fit_report_root_name`
- rationale:
  - the detailed fit file remains compatible with the existing
    `fits/<slice_id>/...` object layout
  - report consumers can open a smaller summary file organized by centrality and
    mT without depending on the full per-slice fit object tree
  - the report path can be redirected independently from the CF/fit output
    directory while older configs keep working through defaults

## DEC-006: Keep Qn-Integrated ME As Default And Make Qn-Following ME Opt-In

- date: 2026-07-01
- decision: `build.split_mixed_event_by_qn = false` remains the default so
  qn-specific SAME slices continue to use qn-integrated mixed-event
  denominators unless explicitly requested
- rationale:
  - preserving the default keeps the existing qn-split production CFs
    numerically backward compatible
  - requiring `build.split_same_event_by_qn = true` prevents a configuration
    from claiming qn-specific ME denominators when no qn-specific SE slices are
    written
  - the opt-in mode supports direct checks of qn-matched SE/ME denominators
    without changing qn labels, group IDs, or fit selection semantics
