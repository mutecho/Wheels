#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_exp_femto_3d_PROFILE.sh [--tier scout|focused-1d|focused-2d|strict]
                              [--stage all|build-cf|fit]
                              [--profile-estimate-only]
                              [--model full|diag] [--input-cf-root <cf.root>]

Options:
  --tier <tier>              Profile tier. Default: scout.
  --stage <stage>            Workflow stage. Default: fit.
  --profile-estimate-only    Validate the selected profile contract and print
                             its estimated work without writing outputs.
  --model <full|diag>        Override the configured fit model for fit stages.
  --input-cf-root <path>     Override the configured CF ROOT input for fit stages.
  --binary <path>            exp_femto_3d executable path.
  -h, --help                 Show this help.

The strict tier is the preserved LIKELYHOODTEST configuration. The other tiers
are profile_only runs with separate output and checkpoint names.
USAGE
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_path="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
project_root="$(cd "${script_dir}/.." && pwd)"
binary_path="${project_root}/bin/exp_femto_3d"
tier="scout"
stage="fit"
model_override=""
input_cf_root=""
profile_estimate_only=false
o2_module="${EXP_FEMTO_3D_O2_MODULE:-O2Physics/latest-master-o2}"
original_arg_count=$#
original_args=("$@")

root_runtime_ready() {
  [[ -n "${ROOTSYS:-}" ]] || return 1
  [[ -n "${ROOT_DYN_PATH:-}" ]] || return 1
  command -v root-config >/dev/null 2>&1 || return 1
  root-config --version >/dev/null 2>&1
}

enter_o2_runtime_if_needed() {
  if root_runtime_ready; then
    return
  fi

  if [[ "${EXP_FEMTO_3D_O2_RUNTIME_ENTERED:-0}" == "1" ]]; then
    echo "[error] O2Physics ROOT runtime is not available after alienv setup." >&2
    exit 1
  fi
  if ! command -v alienv >/dev/null 2>&1; then
    echo "[error] O2Physics ROOT runtime is not available and alienv was not found." >&2
    exit 1
  fi

  export EXP_FEMTO_3D_O2_RUNTIME_ENTERED=1
  if (( original_arg_count > 0 )); then
    exec alienv setenv "${o2_module}" -c bash -lc \
      'exec "$0" "$@"' "${script_path}" "${original_args[@]}"
  fi
  exec alienv setenv "${o2_module}" -c bash -lc \
    'exec "$0" "$@"' "${script_path}"
}

absolute_existing_path() {
  local path="$1"
  if [[ "${path}" == /* ]]; then
    printf '%s\n' "${path}"
    return
  fi
  printf '%s/%s\n' "$(cd "$(dirname "${path}")" && pwd)" "$(basename "${path}")"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tier)
      [[ $# -ge 2 ]] || { echo "[error] Missing value after --tier." >&2; exit 1; }
      tier="$2"
      shift 2
      ;;
    --stage)
      [[ $# -ge 2 ]] || { echo "[error] Missing value after --stage." >&2; exit 1; }
      stage="$2"
      shift 2
      ;;
    --profile-estimate-only)
      profile_estimate_only=true
      shift
      ;;
    --model)
      [[ $# -ge 2 ]] || { echo "[error] Missing value after --model." >&2; exit 1; }
      model_override="$2"
      shift 2
      ;;
    --input-cf-root)
      [[ $# -ge 2 ]] || { echo "[error] Missing value after --input-cf-root." >&2; exit 1; }
      input_cf_root="$2"
      shift 2
      ;;
    --binary)
      [[ $# -ge 2 ]] || { echo "[error] Missing value after --binary." >&2; exit 1; }
      binary_path="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[error] Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

case "${tier}" in
  scout)
    config_path="${project_root}/config/oo_build_and_fit_6bins_profile_scout.toml"
    ;;
  focused-1d)
    config_path="${project_root}/config/oo_build_and_fit_6bins_profile_focused_1d.toml"
    ;;
  focused-2d)
    config_path="${project_root}/config/oo_build_and_fit_6bins_profile_focused_2d.toml"
    ;;
  strict)
    config_path="${project_root}/config/oo_build_and_fit_6bins_checklikelihood.toml"
    ;;
  *)
    echo "[error] Unsupported tier: ${tier}" >&2
    echo "        Expected scout, focused-1d, focused-2d, or strict." >&2
    exit 1
    ;;
esac

case "${stage}" in
  all|build-cf|fit) ;;
  *)
    echo "[error] Unsupported stage: ${stage}" >&2
    exit 1
    ;;
esac
case "${model_override}" in
  ""|full|diag) ;;
  *)
    echo "[error] Unsupported fit model: ${model_override}" >&2
    exit 1
    ;;
esac
if [[ "${stage}" == "build-cf" && ( -n "${model_override}" || -n "${input_cf_root}" || "${profile_estimate_only}" == true ) ]]; then
  echo "[error] Fit-only options require --stage all or --stage fit." >&2
  exit 1
fi
if [[ ! -f "${config_path}" ]]; then
  echo "[error] Config file does not exist: ${config_path}" >&2
  exit 1
fi
if [[ ! -x "${binary_path}" ]]; then
  echo "[error] Executable is missing or not executable: ${binary_path}" >&2
  exit 1
fi
if [[ -n "${input_cf_root}" && ! -f "${input_cf_root}" ]]; then
  echo "[error] Input CF ROOT override does not exist: ${input_cf_root}" >&2
  exit 1
fi

config_path="$(absolute_existing_path "${config_path}")"
binary_path="$(absolute_existing_path "${binary_path}")"
if [[ -n "${input_cf_root}" ]]; then
  input_cf_root="$(absolute_existing_path "${input_cf_root}")"
fi

enter_o2_runtime_if_needed
cd "${project_root}"

run_build_cf() {
  echo "[run] tier=${tier} stage=build-cf" >&2
  "${binary_path}" build-cf --config "${config_path}"
}

run_fit() {
  local command=("${binary_path}" fit --config "${config_path}")
  if [[ -n "${model_override}" ]]; then
    command+=(--model "${model_override}")
  fi
  if [[ -n "${input_cf_root}" ]]; then
    command+=(--input-cf-root "${input_cf_root}")
  fi
  if [[ "${profile_estimate_only}" == true ]]; then
    command+=(--profile-estimate-only)
  fi
  echo "[run] tier=${tier} stage=fit" >&2
  "${command[@]}"
}

case "${stage}" in
  all)
    run_build_cf
    run_fit
    ;;
  build-cf)
    run_build_cf
    ;;
  fit)
    run_fit
    ;;
esac
