# ROOT Runtime Diagnosis For Agents

## Purpose

This note captures a confirmed diagnostic result for `Exp_femto_3d`:

- the refactored project code is not the root cause of the observed `THnSparse`
  runtime failures
- the failure mode is caused by entering the ROOT/O2 environment incorrectly
  from a sandboxed tool process
- agents should treat this as an environment-entry problem first, not as a
  project-code regression

## Problem Statement

ROOT-dependent tests such as `root_runtime_probe`, `slice_catalog_roundtrip_test`,
 and `workflow_smoke_test` previously failed with errors like:

- `/opt/homebrew/bin/alienv: line 315: /dev/fd/62: Operation not permitted`
- `root: command not found` or `root-config: command not found`
- `fatal error: module file '.../BUILD/.../ROOT/lib/Vc.pcm' not found`
- `ROOT PCM .../lib/libMathCore_rdict.pcm file does not exist`
- `no interpreter information for class THnSparseT<TArrayF>`
- segmentation faults during `THnSparseF::Projection(...)`

These errors appeared when tests were started from a sandboxed tool invocation.

## Confirmed Facts

### 1. The project links against the intended ROOT installation

Build metadata resolves ROOT libraries and commands under:

- `/Users/allenzhou/ALICE/sw/osx_arm64/ROOT/v6-36-10-alice1-local1`

The built test executables also contain the expected ROOT `rpath`, so normal
shared-library loading is not the first failure point.

### 2. The crash is triggered at the ROOT interpreter / PCM layer

The failing path is not plain C++ symbol resolution. The breakage starts when
ROOT needs Cling, PCM files, rootmaps, and dictionaries during operations such
as:

- reading persisted `THnSparseF`
- `THnSparseF::Projection(...)`
- ROOT class metadata resolution for histogram and sparse types

### 3. The local ROOT installation works when entered correctly

Running outside the sandbox with a successful:

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '...'
```

confirmed all of the following:

- `which root` resolves correctly
- `which root-config` resolves correctly
- `root-config --version` returns `6.36.10`
- a ROOT snippet constructing `THnSparseF` runs successfully
- the built executable `./workflow_smoke_test` passes successfully

This is the decisive evidence that the host ROOT installation is usable.

## Root Cause

### Primary cause

The sandboxed tool environment prevented `alienv` from completing its normal
environment bootstrap, specifically around its temporary init-script handling:

- `/dev/fd/...: Operation not permitted`

As a result, the spawned shell inherited an incomplete ROOT/O2 runtime
environment.

### What this breaks in practice

An incomplete `alienv` entry can leave the process in a misleading half-working
state:

- CMake can still configure because cached ROOT paths already exist
- executables can still start because `rpath` points at ROOT shared libraries
- but ROOT interpreter/module features fail because the full O2/ROOT runtime
  environment is not actually active

This explains why the observed failures only appear once the code touches
`THnSparse` projection and dictionary-heavy ROOT behavior.

### Why the PCM messages are misleading

Errors mentioning:

- `Vc.pcm`
- `libMathCore_rdict.pcm`
- `libHist_rdict.pcm`
- `libTree_rdict.pcm`

should not be interpreted as proof that the project code requested the wrong
files. In this incident they are downstream symptoms of an environment that was
not entered cleanly enough for ROOT's runtime C++ module system.

## Differential Diagnosis

Agents should use the following decision tree.

### Case A: sandboxed `alienv` shows `/dev/fd/... Operation not permitted`

Interpretation:

- environment entry failed before ROOT was fully usable

Action:

- do not diagnose the project code yet
- rerun the same ROOT/O2 command outside the sandbox or with escalation

### Case B: `which root` or `root-config` fails after `alienv setenv ...`

Interpretation:

- environment entry is incomplete

Action:

- classify as environment-entry failure
- do not attribute later PCM/dictionary errors to the project

### Case C: outside the sandbox the same ROOT snippet and the built test pass

Interpretation:

- the project code is not the root cause
- sandbox execution mode is the differentiator

Action:

- report the issue as a tooling/runtime-entry problem
- keep code changes minimal or avoid unnecessary code-side workarounds

## Recommended Resolution

### For agents running ROOT commands

Prefer this sequence:

1. Run ROOT/O2 commands with the standard O2Physics wrapper.
2. If `alienv` emits `/dev/fd/... Operation not permitted`, treat that as an
   environment-entry failure.
3. Rerun the same command outside the sandbox or with escalated execution.
4. Only investigate project code after confirming the command still fails in a
   clean, non-sandboxed ROOT/O2 environment.

### Canonical command pattern

```bash
/bin/zsh -lc "alienv setenv O2Physics/latest-master-o2 -c sh -lc 'cd /abs/path && <command>'"
```

### Verification command pattern

Use a minimal ROOT check before blaming the project:

```bash
/bin/zsh -lc "alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  which root
  which root-config
  root-config --version
  root -l -b -q -e \"THnSparseF s(\\\"s\\\",\\\"s\\\",1,(int[]){4},(double[]){0.0},(double[]){1.0}); cout << \\\"ok\\\" << endl; gSystem->Exit(0);\"
'"
```

If this succeeds, ROOT itself is operational enough for `Exp_femto_3d`.

## Guidance For The ROOT-Calling Skill

When the ROOT-calling skill or any future agent wrapper evaluates failures, it
should follow these rules:

- classify `/dev/fd/... Operation not permitted` from `alienv` as an
  environment-entry failure
- do not infer that ROOT is broken solely from PCM or dictionary errors that
  occur after a failed sandboxed `alienv` entry
- do not infer that the project code is broken solely from `THnSparse`
  projection crashes observed in that state
- rerun important ROOT-dependent commands outside the sandbox before escalating
  the diagnosis to code-level debugging
- distinguish:
  - shared library loading problems
  - ROOT interpreter / PCM / dictionary problems
  - actual project logic regressions

## Project-Specific Implication

For `Exp_femto_3d`, failures in:

- `root_runtime_probe`
- `slice_catalog_roundtrip_test`
- `workflow_smoke_test`

must be interpreted in context:

- if they fail only in sandboxed execution, the result is not evidence of a
  regression in `RunBuildCf` or `RunFit`
- if they also fail in a clean non-sandboxed O2Physics environment, then code
  debugging becomes justified

## Current Best Practice

When validating this project, prefer:

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build &&
  ctest --output-on-failure
'
```

or run from an already-entered interactive O2Physics shell.

## Status

- diagnosis confidence: high
- local ROOT installation status: confirmed usable in a clean non-sandboxed
  O2Physics environment
- project-code culpability for the observed sandbox failure: not supported by
  current evidence
