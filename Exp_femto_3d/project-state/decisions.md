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
