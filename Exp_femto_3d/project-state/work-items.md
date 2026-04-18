# Work Items

## WI-001: Real-Data Equivalence Validation

- status: open
- owner: next engineering pass
- goal: compare refactored `build-cf` and `fit` outputs against the validated
  legacy macro on a known-good real dataset after the current phi-mapping
  control changes
- success condition:
  - raw counts match within expected floating tolerance
  - CF/projection outputs remain physics-equivalent
  - both phi conventions are checked: follow-input mapping and explicit
    fit-side override
  - key fit parameters and `R2_vs_phi` trends remain acceptable

