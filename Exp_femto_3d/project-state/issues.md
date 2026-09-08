# Issues

## ISSUE-003: Two-Parameter Profile Display Meaning Needs Separate Review

- status: resolved 2026-09-07
- severity: medium
- description:
  the current `rout2_lambda` 2D heatmaps/canvases have unresolved axis and
  interpretation concerns reported by the operator
- resolution:
  implemented the accepted issue-2 contract: exact resolved TH2 domains,
  independent profile/slice validity, invalid masks, four-corner-only contours,
  threshold/full-range canvases, stable best/nominal markers, and persisted
  annotations. Deterministic geometry, ROOT readback, and toy visual QA passed
- remaining follow-up:
  production convergence on a user-run `_r2max64` OO result is not inferred
  from the toy display validation

## ISSUE-002: Profile ROOT Runtime And Canvas QA Blocked

- status: resolved 2026-09-04
- severity: medium
- description:
  the profile implementation, config coverage, and ROOT-independent scan driver
  are present, but the closeout O2Physics ROOT executor could not be entered;
  sandbox execution failed on `/dev/fd`, and the explicitly authorized
  escalation was rejected by the approval backend with HTTP 403
- impact:
  the current pass cannot claim an independently executed ROOT schema smoke,
  full CTest result, or manual toy-canvas visual inspection
- resolution:
  executor access was restored; the primary O2Physics module built successfully,
  CTest passed 7/7, the merged profile ROOT passed structural inspection, and a
  representative toy canvas was visually checked

## ISSUE-001: Real-Data Regression Evidence Still Missing

- status: open
- severity: medium
- description: the local config and ROOT-backed smoke coverage now passes in a
  clean O2Physics environment for the current phi-mapping/progress update, but
  real-data equivalence against the legacy macro has not yet been revalidated
  after this change set
- impact: project behaviour is locally well-covered, but physics-equivalence on
  production-like input is not yet fully closed
