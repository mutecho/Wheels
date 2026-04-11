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
