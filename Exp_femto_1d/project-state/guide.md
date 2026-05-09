# Exp_femto_1d Project Guide

## Purpose

`Exp_femto_1d` is the engineering refactor of the legacy
`legacy/get_cf_from_exp.cpp` macro into a maintainable `CMake + C++17` project
with an explicit two-stage workflow and a narrow `pi-pi` CATS fit model.

## Workflow

1. `build-cf`
   - reads 4D `THnSparseF` input with axes `k* / mT / centrality / event-plane`
   - builds `MinBias`, `InPlane`, and `OutOfPlane` slices
   - stores `SE_raw1d`, `ME_raw1d`, `CF1D`, and `meta/SliceCatalog`
2. `fit`
   - reads `SliceCatalog`
   - fits selected slices with `baseline polynomial x CATS`
   - writes `FitCatalog`, per-slice fit objects, region summaries, and a TSV summary

## Configuration

The project is driven by TOML:

- `[input]`
- `[output]`
- `[build]`
- `[fit]`
- `[[bins.centrality]]`
- `[[bins.mt]]`
- optional `[[fit_selection.*]]`

If `fit_selection` is omitted, fit follows the full build bin grid.
`build.cf_by_mt_show_markers` defaults to `false`; enable it only when marker
symbols are useful on `CFByMtCanvas` trend overlays.

## Run Helper

`scripts/run_exp_femto_1d.sh` is the project-local entry point for the
config-driven run workflow. It defaults to
`config/pbpb_build_and_fit.toml` and runs `build-cf`; use `--stage all` to run
`build-cf` followed by the slower `fit` stage.
Direct shell invocations automatically re-enter `O2Physics/latest-master-o2`
with `alienv` when `ROOTSYS`, `ROOT_DYN_PATH`, or `root-config` are missing.

Supported operator controls:

- `--config <file.toml>` selects the TOML configuration file
- `--stage all|build-cf|fit` selects the workflow stage; default is `build-cf`
- `--input-cf-root <path>` forwards the fit-side CF input override
- `--binary <path>` points to a non-default `exp_femto_1d` executable

The helper expects a built executable. Override the runtime module with
`EXP_FEMTO_1D_O2_MODULE` only when validating against another compatible
O2Physics stack.

## Output Contract

Build output:

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw1d`
- `slices/<slice_id>/ME_raw1d`
- `slices/<slice_id>/CF1D`
- `cent_slices/<cent_id>/<region_name>/CFByMtCanvas`

Fit output:

- `meta/FitCatalog`
- `fits/<slice_id>/DataCF`
- `fits/<slice_id>/FitFunction`
- `fits/<slice_id>/Baseline`
- `fits/<slice_id>/PureFemto`
- `fits/<slice_id>/FitCanvas`
- `summary/by_region/<group_id>/SourceSize|Norm|P1|P2`
- `fit_summary.tsv`

## v1 Limits

- like-sign `pi-pi` only
- fixed `baseline x CATS` model
- no unlike-sign support
- no source-family switching
- no legacy flat object-name compatibility layer

## Validation Expectations

- config parse test should always run
- ROOT-backed tests should run under a fully entered O2Physics environment
- real-data equivalence against the legacy macro remains a required follow-up
- real-data fit stability on at least one `MinBias` and one EP-differential slice remains required
