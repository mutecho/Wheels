#!/usr/bin/env bash

set -euo pipefail

# Resolve paths from the scripts directory back to the project root.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
binary_path="${project_root}/bin/eventgen_femto_3d"

# Stop early when the project has not been built for the current ROOT runtime.
if [[ ! -x "${binary_path}" ]]; then
  echo "Binary not found: ${binary_path}" >&2
  echo "Build the project first with CMake." >&2
  exit 1
fi

# Run from the project root so relative output paths stay predictable.
cd "${project_root}"

# Quote all forwarded arguments before handing them to alienv's shell entry.
quoted_binary="$(printf '%q' "${binary_path}")"
quoted_args=()
for arg in "$@"; do
  quoted_args+=("$(printf '%q' "${arg}")")
done
command_line="${quoted_binary}"
if [[ ${#quoted_args[@]} -gt 0 ]]; then
  command_line+=" ${quoted_args[*]}"
fi

# Enter the O2Physics ROOT runtime and execute the real binary.
exec alienv setenv O2Physics/latest-master-o2 -c sh -lc "${command_line}"
