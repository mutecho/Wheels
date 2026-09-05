# Work Items

## WI-002: Complete Profile ROOT Closeout Validation

- status: complete 2026-09-04
- owner: completed in current engineering pass
- goal: validate the new detailed PML profile mode in the required ROOT runtime
- success condition:
  - configured build succeeds through the shared O2Physics ROOT executor
  - complete CTest passes; any code-77 ROOT tests pass when run directly through
    the same executor
  - toy profile ROOT contains the documented catalogs, point/attempt trees,
    exact coordinates, status masks, nuisance trajectories, and markers
  - valid profiled values obey the toy `profile <= slice` check
  - a full-model non-PSD scan point is classified as
    `model_domain_invalid` and is not displayed as a valid likelihood value
  - representative 1D/2D canvases pass manual visual inspection

## WI-003: Benchmark Real OO Profile Scaling

- status: open
- owner: operator/performance follow-up
- goal: run the documented nominal, 11-point 1D, 9x9 2D, strict 1D, then strict
  2D benchmark sequence on the 85 MB CF with isolated output names
- success condition:
  - record `/usr/bin/time -l`, RSS, I/O, FCN, MIGRAD, and HESSE fractions
  - compare 1/2/4/6 process workers and require at least 1.6x at two workers
    without increased swap pressure before raising the recommended worker count
  - keep real-data physics interpretation separate from the scout tier

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
  - standalone fit report ROOT file contains the expected source-parameter
    overview canvases and `eps_vs_mt` graphs for the real selected cent/mT bins
