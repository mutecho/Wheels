#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
binary_path="${script_dir}/build/eventgen_femto_3d"

if [[ ! -x "${binary_path}" ]]; then
  echo "Binary not found: ${binary_path}" >&2
  echo "Build the project first with CMake." >&2
  exit 1
fi

cd "${script_dir}"

quoted_binary="$(printf '%q' "${binary_path}")"
quoted_args=()
for arg in "$@"; do
  quoted_args+=("$(printf '%q' "${arg}")")
done
command_line="${quoted_binary}"
if [[ ${#quoted_args[@]} -gt 0 ]]; then
  command_line+=" ${quoted_args[*]}"
fi

exec alienv setenv O2Physics/latest-master-o2 -c sh -lc "${command_line}"
