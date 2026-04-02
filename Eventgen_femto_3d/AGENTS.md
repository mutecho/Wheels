# Project Engineering Workflow

This repository uses a structured engineering workflow. Do not jump directly to implementation unless the task is already clear and bounded.

## Workflow skills available

Use these workflow skills when appropriate:
- analyze-requirement
- implement-feature
- generate-test-matrix
- run-eval-and-debug
- review-change

## Default workflow

Not every task needs every stage, but use this as the default model:

- unclear or unstable request -> `analyze-requirement`
- bounded implementation task -> `implement-feature`
- coverage planning needed -> `generate-test-matrix`
- existing failing validation or repro -> `run-eval-and-debug`
- concrete diff ready for assessment -> `review-change`

## ROOT

Use ROOT when needed, ROOT can be accessed by related skills, including:

- o2physics-root
- root-file-inspector
- root-histogram-ops


## Required handoff artifacts

When a stage is used, preserve its structured output for the next stage:

- analyze-requirement -> structured task card
- implement-feature -> implementation handoff
- generate-test-matrix -> structured test matrix
- run-eval-and-debug -> structured debug report
- review-change -> structured review report

Do not collapse these into a loose summary.

## Required workflow rules

- If scope, acceptance criteria, or constraints are unstable, start with `analyze-requirement`.
- Do not begin bounded implementation without a clear task definition.
- For medium- or high-risk changes, produce a test matrix before claiming readiness.
- If there is already a failing test, eval script, or repro path, use `run-eval-and-debug` rather than speculative debugging.
- After implementation, produce a reviewable handoff before concluding the task.
- Do not claim success without validation evidence.
- If validation could not run, mark the result as unverified.

## Build, test, and verification

Before considering a change complete, run the relevant commands for this repo:
- build: <fill>
- unit tests: <fill>
- integration tests: <fill>
- lint/typecheck: <fill>

## Done means

A task is only done when:
- the requested behavior is implemented or the limitation is stated clearly
- relevant validation has run, or the result is explicitly marked unverified
- risks and assumptions are visible
- the next reviewer can inspect the change without reconstructing intent from scratch

## Purpose and refrence

Purpose and requirements can be found in the Code_structure directory