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
