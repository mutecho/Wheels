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

## Output Contract

Build output:

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw1d`
- `slices/<slice_id>/ME_raw1d`
- `slices/<slice_id>/CF1D`

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
