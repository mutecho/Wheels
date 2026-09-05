# Exp_femto_3d

Refactored C++17/CMake version of the original `3d_cf_from_exp.cpp` ROOT macro.

## Layout

- `include/exp_femto_3d/`: public types and workflow interfaces
- `src/`: configuration parsing, logging, CF building, and Levy fitting
- `app/`: CLI entry point
- `docs/`: human-facing workflow notes, including the math/physics formula flow
- `config/examples/`: example TOML configurations
- `tests/`: config and workflow verification helpers
- `legacy/`: archived macro reference kept out of the CMake build

## Documentation

- [`docs/数学物理公式流程说明.md`](docs/数学物理公式流程说明.md): current
  math/physics workflow reference from sparse inputs through CF building,
  Levy/Coulomb fitting, and summary outputs

## Build

This project expects ROOT from the ALICE/O2 environment and `toml++` from
Homebrew.

```bash
cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_3d \
      -B /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build
cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build
```

Configured build trees still live under `build/`, but all executable targets are
emitted into `bin/`. `compile_commands.json` is exported in the configured build
directory as `build/compile_commands.json`.

## Run

```bash
./bin/exp_femto_3d build-cf --config config/pbpb_build_and_fit.toml
./bin/exp_femto_3d fit --config config/pbpb_build_and_fit.toml
```

`scripts/run_exp_femto_3d.sh` defaults to `config/oo_build_and_fit.toml` and
`--stage all`; either value can be overridden with the documented CLI flags.

Optional fit overrides:

- `--model full|diag`
- `--input-cf-root /absolute/or/relative/path.root`

Coulomb fit modes:

- `[fit].coulomb_mode = "none" | "gamow" | "finite_source"`
- `[fit].finite_source_mode = "fixed_1d" | "iterative_1d"` is only valid with
  `coulomb_mode = "finite_source"`
- older `[fit].use_coulomb = false|true` configs remain accepted and map to
  `none` and `gamow`; when both fields are present they must agree
- finite-source mode uses local CATS support when CMake finds
  `/Users/allenzhou/Research_software/CATS/install` and Homebrew GSL; builds
  without CATS still support `none` and `gamow`

Main Levy fit-parameter overrides:

- optional `[fit.parameters.<name>]` subtables can override `initial`, `min`,
  and `max`; omitted fields keep the legacy defaults
- supported names: `norm`, `lambda`, `rout2`, `rside2`, `rlong2`,
  `routside2`, `routlong2`, `rsidelong2`, `alpha`, and `baseline_q2`
- `lambda` and `alpha` additionally support `fixed_value`; other parameters
  reject `fixed_value`
- `[fit.parameters.lambda]` requires `use_core_halo_lambda = true`; otherwise
  lambda is fixed by the core-halo switch
- `[fit.parameters.baseline_q2]` requires `use_q2_baseline = true`

```toml
[fit.parameters.lambda]
initial = 0.5
min = 0.0
max = 1.0
fixed_value = 0.65

[fit.parameters.alpha]
initial = 1.5
min = 0.5
max = 2.0
```

Independent fit-report output:

- `[output].fit_report_directory`
- `[output].fit_report_root_name`
- these fields control the standalone ROOT report file produced by `fit`; when
  `fit_report_directory` is omitted, it defaults to `output_directory`

PML profile-likelihood diagnostics:

- `[fit.profile_likelihood].enabled = false` is the default; when disabled,
  `fit` does not create or reset the profile ROOT file
- the mode requires `[fit].use_pml = true` and runs inside the existing `fit`
  command without changing the production fit result, `FitCatalog`, or TSV
- one run may contain multiple 1D/2D scans over exact `slice_ids`, or over all
  slices selected by `fit_selection`
- every grid point fixes the named scan parameter(s) and re-minimizes all
  remaining free nuisance parameters; `Slice1D` is written separately for the
  optional fixed-nuisance comparison
- finite-source scans reuse the group Coulomb kernel prepared by the nominal
  workflow; the kernel is never rebuilt as a scanned radius changes
- `[output].profile_root_name` controls the independent diagnostic file and
  must not resolve to the CF, detailed-fit, or report ROOT path
- profile contours are labelled `diagnostic only`; they are not confidence
  intervals and are not labelled as 68%/95% CL regions
- `execution_mode = "profile_only"` requires explicitly listed slices and never
  creates or modifies the detailed-fit ROOT, report ROOT, TSV, or `FitCatalog`
- `parallel_backend = "process"` assigns complete `(centrality,mT,qn)` Coulomb
  groups to isolated legacy-TMinuit workers; scans within one slice remain serial
- process checkpoints are written per completed group and reused only when the
  input-CF content and complete scan/model contract digest match
- `hesse_strategy = "none"` skips HESSE for the private nominal fit and every
  profile attempt; stored errors are NaN with `parameter_errors_valid=false`
- the Minuit2/thread backend is parsed as an experimental contract but execution
  is deliberately rejected until numerical A/B validation is completed

```toml
[output]
profile_root_name = "profile_likelihood.root"

[fit.profile_likelihood]
enabled = true
execution_mode = "profile_only"       # "alongside_fit" | "profile_only"
parallel_backend = "process"          # "serial" | "process" | experimental "thread"
workers = 2
minimizer_backend = "legacy_tminuit"
hesse_strategy = "none"               # "all_attempts" | "none"
slice_scope = "listed" # or "fit_selection"
slice_ids = ["exact_slice_id"] # remove this field for fit_selection
retry_strategy = "reference_and_bidirectional_neighbors"
write_likelihood_slice = true
contour_levels = [1.0, 2.0, 4.0]

[fit.profile_likelihood.checkpoint]
enabled = true
resume = true
run_id = "oo_6phi_scout_v1"
directory = "profile_likelihood_scout.work"

[[fit.profile_likelihood.scans]]
id = "rout2"
parameters = ["rout2"]
points = [41]
min = [0.01]
max = [100.0]
refine = true
refinement_points = [21]

[[fit.profile_likelihood.scans]]
id = "rout2_lambda"
parameters = ["rout2", "lambda"]
points = [21, 21]
refine = false
```

Finite scan bounds may be omitted to inherit the active fit bounds. Parameters
without finite default bounds, including full-model off-diagonal radii, require
explicit `min`/`max`. Targets must exist in the active `diag`/`full` model and
must be free in the nominal fit. `refine = true` requires a matching
`refinement_points` entry for every scan axis.

Use `fit --profile-estimate-only` to perform config, catalog, effective-model,
slice, range, grid, and output-collision validation and print the upper bounds
on Minuit attempts and fixed-nuisance evaluations without creating any output.
The staged runner `scripts/run_exp_femto_3d_PROFILE.sh` provides `scout`,
`focused-1d`, `focused-2d`, and preserved `strict` tiers. Start with:

```bash
scripts/run_exp_femto_3d_PROFILE.sh --tier scout --profile-estimate-only
scripts/run_exp_femto_3d_PROFILE.sh --tier scout
```

TOML progress control:

- `[build].progress`
- `[fit].progress`
- accepted values: `true`, `false`, or `"auto"` (also `"enabled"` / `"disabled"`)
- default when omitted: `"auto"`; the bar only appears when `stderr` is attached to a TTY
- enabled progress lines include the stage label, percentage, an activity frame,
  and an ETA computed from completed slices; a heartbeat refreshes the line once
  per second during long ROOT operations

Axis rebin control:

- `[build.rebin.mt]` and `[build.rebin.phi]` make the selected axis explicit and
  show at a glance whether rebin is enabled
- `enabled = false` keeps the native input-axis bins and forbids `factor`,
  `min`/`max`, or `[[bins.*]]` output ranges
- `enabled = true` requires exactly one of:
  - `factor = <integer>` for grouped native bins, with an optional aligned
    `min`/`max` window
  - `[[bins.mt]]` or `[[bins.phi]]` for explicit output ranges that align with
    the input sparse axis
- factor mode merges contiguous normal bins before projection; the workflow does
  not call ROOT histogram `Rebin` on the final `TH3D`
- when symmetric phi mapping is enabled, every non-integrated phi output range
  must stay on one side of the `pi/2` mapping seam; this applies to native,
  factor, and explicit-range modes
- old configs stay valid: phi defaults to native bins when no explicit rebin
  table is present, and mT keeps the legacy `[[bins.mt]]` contract

For example, phi can use a factor while mT uses explicit physical ranges:

```toml
[build.rebin.phi]
enabled = true
factor = 2
min = 0.0
max = 3.141592653589793

[build.rebin.mt]
enabled = true

[[bins.mt]]
min = 0.2
max = 0.6
label = "mt_020_060"
```

Phi mapping control:

- `[build].map_pair_phi_to_symmetric_range` controls how build writes `display_phi_*`
  into `meta/SliceCatalog`
- `[fit].map_pair_phi_to_symmetric_range` is optional; when omitted, fit follows the
  input CF metadata from `SliceCatalog`
- when `[fit].map_pair_phi_to_symmetric_range` is set explicitly, fit reinterprets the
  stored slices from `raw_phi_*`, so the summary `R2_vs_phi` graphs can use either the
  raw `[0, pi]` range or the symmetric `[-pi/2, pi/2]` range without rebuilding the CF
- if the current config's `[build]` mapping flag disagrees with the input CF metadata,
  fit warns and trusts the input CF file
- when mT groups are generated dynamically by native/factor rebin, omitted
  `[[fit_selection.mt]]` means all catalog mT groups; explicit
  `[[fit_selection.mt]]` entries are validated against `meta/SliceCatalog`
  before the fit output ROOT file is reset

Mixed-event denominator control:

- `[build].split_mixed_event_by_phi` defaults to `false`
- `false` preserves the historical behavior: one mixed-event denominator is
  projected over all phi bins for each centrality/mT group, then reused by the
  phi-integrated and per-phi SE slices
- `true` projects the mixed-event denominator with the same phi range as the SE
  slice; the phi-integrated slice still uses the full phi range
- this switch is independent of `map_pair_phi_to_symmetric_range`, which only
  controls stored/displayed phi coordinates for downstream summaries
- `[build].split_same_event_by_qn` defaults to `false`
- `true` requires `[[bins.qn]]` entries and writes the legacy qn-integrated
  slices plus one additional SAME-side qn slice set per configured qn label
- qn-specific groups append the qn label to `group_id`; the qn-integrated
  group keeps the historical cent/mT `group_id` so existing paths remain valid
- `[build].split_mixed_event_by_qn` defaults to `false`, preserving the current
  behavior where qn-specific SAME slices still use a qn-integrated ME
  denominator
- `true` requires `split_same_event_by_qn = true`; ME is then projected with the
  same qn interval as the current qn-specific SAME slice, while the qn-all slice
  remains integrated over the full qn axis

## Output Contract

- `build-cf` writes `meta/SliceCatalog` plus `slices/<slice_id>/...`; the
  catalog records `build_uses_symmetric_phi_range`, both
  `split_mixed_event_by_phi` and `split_mixed_event_by_qn`, the four mT/phi
  rebin fields, and qn metadata
  `qn_index/qn_low/qn_high/qn_label/is_qn_integrated`; older catalogs default
  to qn-integrated ME, legacy-range mT, native phi, and `qn_all`
- `fit` reads `meta/SliceCatalog` and writes `meta/FitCatalog`,
  `meta/CoulombKernelCatalog`, `fits/<slice_id>/...`, `summary/R2_vs_phi/...`,
  and `fit_summary.tsv`; `FitCatalog`, `CoulombKernelCatalog`, and TSV rows
  carry the same qn metadata, and `FitCatalog`/TSV rows carry the same rebin
  metadata and final physical phi ranges
- `fit` also writes the standalone report ROOT file configured by
  `fit_report_directory` and `fit_report_root_name`; this report mirrors
  `meta/FitCatalog` and `summary/R2_vs_phi/...`, adds
  `source_parameters/<cent>/<mt>/source_parameters_overview_canvas`, and adds
  `eps_vs_mt/<cent>/epsf_vs_mt(_canvas)` summaries. qn-specific report
  summaries are written below an extra `<qn_label>` directory, while qn-all
  summaries keep the historical paths.

## Test Notes

`ctest` always runs config-only tests. ROOT integration tests are guarded by a
runtime probe and are skipped when the local ROOT runtime cannot safely create,
persist, and project a `THnSparseF`.

For agent-oriented diagnosis of sandboxed ROOT/O2 failures, see
[`ROOT_RUNTIME_AGENT_NOTE.md`](./ROOT_RUNTIME_AGENT_NOTE.md).
