#!/usr/bin/env bash
set -euo pipefail

# Resolve paths from this script so the helper works from any current directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
build_dir="${EVENTGEN_FEMTO_3D_BUILD_DIR:-${project_root}/build}"
jobs="${EVENTGEN_FEMTO_3D_BUILD_JOBS:-8}"
clean_first="${EVENTGEN_FEMTO_3D_CLEAN_FIRST:-0}"
cache_file="${build_dir}/CMakeCache.txt"
refresh_root_cache=0
current_root_cmake_dir=""
configured_root_dir=""

if command -v root-config >/dev/null 2>&1; then
  current_root_prefix="$(root-config --prefix)"
  candidate_root_cmake_dir="${current_root_prefix}/cmake"
  if [[ -f "${candidate_root_cmake_dir}/ROOTConfig.cmake" ]]; then
    current_root_cmake_dir="${candidate_root_cmake_dir}"
  fi
fi

if [[ -n "${current_root_cmake_dir}" && -f "${cache_file}" ]]; then
  configured_root_dir="$(awk -F= '/^ROOT_DIR:[^=]*=/{print $2; exit}' "${cache_file}")"
  if [[ -n "${configured_root_dir}" && "${configured_root_dir}" != "${current_root_cmake_dir}" ]]; then
    refresh_root_cache=1
    echo "[info] ROOT CMake cache changed:" >&2
    echo "       cached:  ${configured_root_dir}" >&2
    echo "       current: ${current_root_cmake_dir}" >&2
  fi
fi

# Reconfigure the selected build tree before compiling.
cmake_args=(-S "${project_root}" -B "${build_dir}")
if [[ -n "${current_root_cmake_dir}" ]]; then
  cmake_args=(-DROOT_DIR:PATH="${current_root_cmake_dir}" "${cmake_args[@]}")
fi
if [[ "${refresh_root_cache}" == "1" ]]; then
  cmake_args=(-U "ROOT_*" -U "ROOT_DIR" "${cmake_args[@]}")
fi
cmake "${cmake_args[@]}"

# If the active ROOT runtime changed, refresh link script timestamps so CMake
# relinks generated binaries without recompiling unchanged translation units.
if [[ "${refresh_root_cache}" == "1" && -d "${build_dir}/CMakeFiles" ]]; then
  find "${build_dir}/CMakeFiles" -path '*/link.txt' -type f -exec cmake -E touch {} +
fi

# Default to CMake's incremental dependency checks. Set
# EVENTGEN_FEMTO_3D_CLEAN_FIRST=1 only for a deliberately full clean rebuild.
build_args=(--parallel "${jobs}")
if [[ "${clean_first}" == "1" ]]; then
  build_args=(--clean-first "${build_args[@]}")
fi
cmake --build "${build_dir}" "${build_args[@]}"
