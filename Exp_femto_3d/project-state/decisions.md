# Decisions

## DEC-011: Parallelize Complete Groups In Processes, Keep Scans Sequential

- date: 2026-09-04
- decision:
  legacy TMinuit profiles may run concurrently only as isolated process tasks
  owning complete `(centrality,mT,qn)` groups; each slice/scan remains serial,
  and `profile_only` is required for the process backend
- rationale:
  - TMinuit, ROOT, PML, and Coulomb callback state are process-local rather than
    thread-safe in the current implementation
  - group ownership builds and freezes one finite-source kernel in one worker
  - `posix_spawn` avoids inheriting ROOT locks, which a post-ROOT `fork` can do
  - chunk validation and atomic final replacement preserve the last complete
    diagnostic output across interruption or worker failure
- deferred alternative:
  Minuit2/thread execution stays gated until grid minima, status classes, deltas,
  and 1/2/4-worker scaling pass the documented A/B matrix

## DEC-010: Keep Profile Likelihood Separate, Serial, And Diagnostic-Only

- date: 2026-09-03
- decision:
  add profile likelihood as a default-off mode of the existing `fit` command;
  write it to a separate ROOT file, keep scans serial, compute deltas from the
  lowest valid nominal/profile point per slice, and never replace production
  fit results with a profile minimum
- rationale:
  - the production PML statistic, physics model, catalog success semantics, and
    ordinary output files remain backward compatible when the mode is disabled
  - serial evaluation avoids races in the current global PML/Coulomb context
  - numerical trees can preserve failed attempts and exact coordinates without
    forcing invalid points into regular histogram bins
  - diagnostic contours must not be presented as confidence intervals
- finite-source constraint:
  every profile point uses the final group Coulomb kernel prepared by the
  nominal workflow; scanning a radius does not regenerate the kernel

## DEC-006: Keep Build-Side Rebin Explicit And Default Phi To Native

- date: 2026-08-12
- decision: require explicit `[build.rebin.mt]` / `[build.rebin.phi]` tables to
  enable build-side rebin; when the phi table is omitted, phi stays on native
  bins and the workflow records that fact as `phi_rebin = native`
- rationale:
  - the visible `enabled` switch makes build behavior obvious in config reviews
  - native phi remains the default for old configs, so the new contract does not
    silently change existing production behavior
  - the explicit switch keeps future rebin changes in one place instead of
    encoding them indirectly through bin names alone
  - old configs without the new tables keep native phi and legacy
    `[[bins.mt]]` semantics

## DEC-007: Rebin Sparse Selection Axes Before q-Space Projection

- date: 2026-08-12
- decision: implement mT/phi rebin as closed normal-bin selections on the input
  sparse axes before projecting q_out/q_side/q_long; do not call ROOT `Rebin`
  on the output `TH3D`
- rationale:
  - this preserves the existing q-axis bins and fit coordinates
  - `split_mixed_event_by_phi = true` can apply the identical merged phi
    selection to SE and ME, while `phi_all` remains the complete normal range
  - explicit ranges may overlap for integral selections, but exact duplicates,
    non-aligned edges, non-divisible factors, and symmetric-map intervals that
    cross `pi/2` are rejected

## DEC-008: Validate Before Resetting Existing ROOT Outputs

- date: 2026-08-12
- decision: resolve input axes, rebin selections, phi seams, output identifiers,
  SliceCatalog availability, and explicit dynamic mT fit selections before
  resetting the corresponding CF or fit ROOT output
- rationale:
  - a configuration error must not destroy the last valid analysis artifact
  - sentinel regression tests make this failure-ordering contract durable

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

## DEC-009: Keep Qn-Integrated ME As Default And Make Qn-Following ME Opt-In

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
