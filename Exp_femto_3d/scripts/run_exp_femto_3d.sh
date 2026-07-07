#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_exp_femto_3d.sh [--config <file.toml>] [--stage all|build-cf|fit]
  run_exp_femto_3d.sh --stage fit --config <file.toml> [--model full|diag] [--input-cf-root <cf.root>]

Options:
  --config <file.toml>      TOML configuration file.
                            Default: <project>/config/oo_build_and_fit.toml
  --stage <stage>           Workflow stage to run: all, build-cf, or fit.
                            Default: all
  --model <full|diag>       Override the configured fit model for fit stages.
  --input-cf-root <path>    Override the configured CF ROOT input for fit stages.
  --binary <path>           exp_femto_3d executable path.
                            Default: <project>/bin/exp_femto_3d
  -h, --help                Show this help.

When needed, this helper re-enters the O2Physics ROOT runtime with alienv before
running the ROOT-backed executable.
USAGE
}

# Resolve project-local defaults from the script path so callers can run from any cwd.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_path="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
project_root="$(cd "${script_dir}/.." && pwd)"
# config_path="${project_root}/config/oo_build_and_fit.toml"
config_path="${project_root}/config/NeNe_build_and_fit.toml"
binary_path="${project_root}/bin/exp_femto_3d"
stage="all"
# stage="build"
# stage="fit"
model_override=""
input_cf_root=""
o2_module="${EXP_FEMTO_3D_O2_MODULE:-O2Physics/latest-master-o2}"
original_arg_count=$#
original_args=("$@")

# Treat the runtime as ready only when ROOT command-line tools and module paths are active.
root_runtime_ready() {
  [[ -n "${ROOTSYS:-}" ]] || return 1
  [[ -n "${ROOT_DYN_PATH:-}" ]] || return 1
  command -v root-config >/dev/null 2>&1 || return 1
  root-config --version >/dev/null 2>&1
}

# Re-exec through alienv for direct shell invocations outside the O2Physics runtime.
enter_o2_runtime_if_needed() {
  if root_runtime_ready; then
    return
  fi

  if [[ "${EXP_FEMTO_3D_O2_RUNTIME_ENTERED:-0}" == "1" ]]; then
    echo "[error] O2Physics ROOT runtime is not available after alienv setup." >&2
    echo "        Expected ROOTSYS, ROOT_DYN_PATH, and a working root-config." >&2
    exit 1
  fi

  if ! command -v alienv >/dev/null 2>&1; then
    echo "[error] O2Physics ROOT runtime is not available and alienv was not found." >&2
    echo "        Enter O2Physics manually, or install alienv on PATH." >&2
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

# Convert existing files to absolute paths before the script changes directory.
absolute_existing_path() {
  local path="$1"
  local dir
  local base

  if [[ "${path}" == /* ]]; then
    printf '%s\n' "${path}"
    return
  fi

  dir="$(dirname "${path}")"
  base="$(basename "${path}")"
  printf '%s/%s\n' "$(cd "${dir}" && pwd)" "${base}"
}

# Parse operator overrides before validating paths or starting ROOT-backed work.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 ]]; then
        echo "[error] Missing value after --config." >&2
        exit 1
      fi
      config_path="$2"
      shift 2
      ;;
    --stage)
      if [[ $# -lt 2 ]]; then
        echo "[error] Missing value after --stage." >&2
        exit 1
      fi
      stage="$2"
      shift 2
      ;;
    --model)
      if [[ $# -lt 2 ]]; then
        echo "[error] Missing value after --model." >&2
        exit 1
      fi
      model_override="$2"
      shift 2
      ;;
    --input-cf-root)
      if [[ $# -lt 2 ]]; then
        echo "[error] Missing value after --input-cf-root." >&2
        exit 1
      fi
      input_cf_root="$2"
      shift 2
      ;;
    --binary)
      if [[ $# -lt 2 ]]; then
        echo "[error] Missing value after --binary." >&2
        exit 1
      fi
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

# Validate the selected run contract before entering the expensive workflow.
case "${stage}" in
  all|build-cf|fit)
    ;;
  *)
    echo "[error] Unsupported stage: ${stage}" >&2
    echo "        Expected one of: all, build-cf, fit." >&2
    exit 1
    ;;
esac

case "${model_override}" in
  ""|full|diag)
    ;;
  *)
    echo "[error] Unsupported fit model: ${model_override}" >&2
    echo "        Expected one of: full, diag." >&2
    exit 1
    ;;
esac

if [[ "${stage}" == "build-cf" && -n "${model_override}${input_cf_root}" ]]; then
  echo "[error] --model and --input-cf-root only apply to fit stages." >&2
  exit 1
fi

if [[ ! -f "${config_path}" ]]; then
  echo "[error] Config file does not exist: ${config_path}" >&2
  exit 1
fi

if [[ ! -x "${binary_path}" ]]; then
  echo "[error] Executable is missing or not executable: ${binary_path}" >&2
  echo "        Build the project first under the O2Physics ROOT runtime." >&2
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

# Run from the project root so relative outputs in configs stay predictable.
cd "${project_root}"

# Keep stage boundaries visible so long ROOT work is easy to follow in logs.
announce_stage() {
  echo "[run] starting ${1}" >&2
}

# Keep the build-cf invocation identical to the public CLI contract.
run_build_cf() {
  announce_stage "build-cf"
  "${binary_path}" build-cf --config "${config_path}"
}

# Keep fit-specific overrides local to the fit stage.
run_fit() {
  local command=("${binary_path}" fit --config "${config_path}")

  announce_stage "fit"
  if [[ -n "${model_override}" ]]; then
    command+=(--model "${model_override}")
  fi
  if [[ -n "${input_cf_root}" ]]; then
    command+=(--input-cf-root "${input_cf_root}")
  fi

  "${command[@]}"
}

# Execute the requested stage sequence and let the CLI print compact summaries.
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
