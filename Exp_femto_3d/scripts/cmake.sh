#!/usr/bin/env bash
set -euo pipefail

# Resolve paths from this script so the helper works from any current directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
build_dir="${EXP_FEMTO_3D_BUILD_DIR:-${project_root}/build}"
jobs="${EXP_FEMTO_3D_BUILD_JOBS:-8}"
clean_first="${EXP_FEMTO_3D_CLEAN_FIRST:-1}"

# Reconfigure the selected build tree before compiling.
cmake -S "${project_root}" -B "${build_dir}"

# The project writes executables into source-tree bin/, so clean by default to
# avoid accepting a newer binary produced by another build tree.
build_args=(--parallel "${jobs}")
if [[ "${clean_first}" != "0" ]]; then
  build_args=(--clean-first "${build_args[@]}")
fi
cmake --build "${build_dir}" "${build_args[@]}"

# If this build was configured with CATS, verify the operator binary reflects it.
binary_path="${project_root}/bin/exp_femto_3d"
link_script="${build_dir}/CMakeFiles/exp_femto_3d.dir/link.txt"
if [[ -f "${link_script}" ]] && grep -q 'libCATS' "${link_script}"; then
  if ! otool -L "${binary_path}" | grep -q 'libCATS'; then
    echo "[error] ${binary_path} is not linked against CATS after the build." >&2
    echo "        Rerun with EXP_FEMTO_3D_CLEAN_FIRST=1 or inspect ${link_script}." >&2
    exit 1
  fi
fi
