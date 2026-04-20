# Open Issues

## ISSUE-001 Real-data equivalence not closed

- status: open
- summary: the new structured `build-cf` output has not yet been numerically
  compared against the legacy macro on a known-good real dataset
- impact: physics-level equivalence remains unproven
- next action: run the new `build-cf` and the legacy macro on the same dataset,
  then compare `SE_raw1d`, `ME_raw1d`, and `CF1D` slice-by-slice

## ISSUE-002 Real-data fit inspection still pending

- status: open
- summary: real-data fit stability and canvas readability have not yet been checked
- impact: build/test correctness is established, but physics-facing fit quality is not yet closed
- next action: run one `MinBias` and one EP-differential real-data fit and inspect `FitCanvas`, baseline, and pure-femto separation
