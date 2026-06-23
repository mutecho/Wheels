#!/usr/bin/env bash
set -euo pipefail

# Resolve paths from this script so the helper works from any current directory.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
build_dir="${EXP_FEMTO_3D_BUILD_DIR:-${project_root}/build}"
jobs="${EXP_FEMTO_3D_BUILD_JOBS:-8}"
clean_first="${EXP_FEMTO_3D_CLEAN_FIRST:-0}"
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
if [[ "${refresh_root_cache}" == "1" ]]; then
  cmake_args=(-U "ROOT_*" -U "ROOT_DIR" -DROOT_DIR:PATH="${current_root_cmake_dir}" "${cmake_args[@]}")
fi
cmake "${cmake_args[@]}"

# If the active ROOT runtime changed, refresh link script timestamps so CMake
# relinks generated binaries without recompiling unchanged translation units.
if [[ "${refresh_root_cache}" == "1" && -d "${build_dir}/CMakeFiles" ]]; then
  find "${build_dir}/CMakeFiles" -path '*/link.txt' -type f -exec cmake -E touch {} +
fi

binary_path="${project_root}/bin/exp_femto_3d"
link_script="${build_dir}/CMakeFiles/exp_femto_3d.dir/link.txt"
if [[ -f "${link_script}" && -f "${binary_path}" ]] && grep -q 'libCATS' "${link_script}"; then
  if ! otool -L "${binary_path}" | grep -q 'libCATS'; then
    echo "[info] CATS link rule changed; refreshing CATS-linked targets." >&2
    while IFS= read -r cats_link_script; do
      cmake -E touch "${cats_link_script}"
    done < <(find "${build_dir}/CMakeFiles" -path '*/link.txt' -type f -exec grep -l 'libCATS' {} +)
  fi
fi

# Default to CMake's incremental dependency checks. Set
# EXP_FEMTO_3D_CLEAN_FIRST=1 only for a deliberately full clean rebuild.
build_args=(--parallel "${jobs}")
if [[ "${clean_first}" == "1" ]]; then
  build_args=(--clean-first "${build_args[@]}")
fi
cmake --build "${build_dir}" "${build_args[@]}"

# If this build was configured with CATS, verify the operator binary reflects it.
if [[ -f "${link_script}" ]] && grep -q 'libCATS' "${link_script}"; then
  if ! otool -L "${binary_path}" | grep -q 'libCATS'; then
    echo "[error] ${binary_path} is not linked against CATS after the build." >&2
    echo "        Rerun with EXP_FEMTO_3D_CLEAN_FIRST=1 or inspect ${link_script}." >&2
    exit 1
  fi
fi
