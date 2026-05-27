# Decisions

## DEC-001 Adopt `project-state/`

- status: accepted
- reason: the task changes algorithm flow, config contracts, output schema, test
  strategy, and user run commands

## DEC-002 v1 supports only like-sign `pi-pi`

- status: accepted
- reason: locked by `PLAN.md`

## DEC-003 v1 model is `baseline polynomial x CATS`

- status: accepted
- reason: locked by `PLAN.md`

## DEC-004 Depend only on installed public CATS headers and `libCATS.dylib`

- status: accepted
- reason: the task explicitly forbids private helper dependencies from the CATS source tree

## DEC-005 Do not preserve legacy flat output object names

- status: accepted
- reason: structured `meta/ + slices/ + fits/ + summary/` output is the chosen public schema

## DEC-006 Keep event-plane slicing fixed to `MinBias / InPlane / OutOfPlane`

- status: accepted
- reason: locked by `PLAN.md`

## DEC-007 Add centrality-level CF overlay canvases to build output

- status: accepted
- reason: users need a direct per-centrality, per-event-plane visual comparison
  of all mT-bin `CF1D` histograms without changing the per-slice histogram
  contract

## DEC-008 Hide CF-by-mT overlay markers by default

- status: accepted
- reason: dense markers and error bars obscure the CF trend; `CFByMtCanvas`
  should default to line-only trend overlays while preserving an explicit
  `build.cf_by_mt_show_markers` opt-in

## DEC-009 Keep integrated ME as default and make EP-following ME opt-in

- status: accepted
- reason: existing CF production reused one MinBias mixed-event denominator per
  centrality/mT group; preserving that as `build.split_mixed_event_by_phi =
  false` keeps backward compatibility while allowing an explicit
  event-plane-following denominator mode for differential checks

## DEC-010 Rebin CF operands without changing raw SE/ME outputs

- status: accepted
- reason: `build.cf_rebin_factor` lets `CF1D` be built from coarser same-event
  and mixed-event operands for statistical stability, while `SE_raw1d` and
  `ME_raw1d` remain raw-window histograms for legacy comparisons and diagnosis
