# ROOT Runtime Agent Note

`Exp_femto_1d` depends on the full ROOT runtime from the ALICE/O2Physics
environment and on an installed CATS runtime. When the shell has not entered
that environment, CMake cannot resolve `ROOTConfig.cmake` and ROOT-backed tests
are not authoritative.

## Expected Validation Environment

Recommended authoritative flow:

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cmake -S /Users/allenzhou/Research_software/Code_base/Exp_femto_1d \
        -B /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/build &&
  cmake --build /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/build &&
  cd /Users/allenzhou/Research_software/Code_base/Exp_femto_1d/build &&
  ctest --output-on-failure
'
```

## Sandbox Caveat

If an agent or tool runs `alienv` inside a sandboxed context, environment
bootstrap may fail before ROOT is fully usable. Typical symptoms include:

- missing `ROOTConfig.cmake`
- dictionary or PCM load failures
- projection calls failing after a partial ROOT startup

Those failures should be treated as environment-entry problems first, not as
proof of a code regression.

## Practical Guidance

- trust non-sandboxed O2/ROOT validation over sandbox noise
- keep `tests/run_root_guard.sh` in place so ROOT-backed tests skip instead of
  reporting false failures when the runtime is absent
- do not mark real-data equivalence or fit-quality closure as complete until a
  non-sandboxed O2Physics run confirms them
