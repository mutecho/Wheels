# Legacy Reference

This directory is reserved for legacy macro helpers that should remain out of the
main CMake build and executable entry points.

The original ROOT macro entry point is preserved here as
`legacy/3d_cf_from_exp.cpp` so the refactored project can keep a stable
reference implementation without coupling the new CLI workflow to macro-style
execution.
