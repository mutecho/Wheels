#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_exp_femto_1d.sh [--config <file.toml>] [--stage all|build-cf|fit]
  run_exp_femto_1d.sh --stage fit --config <file.toml> [--input-cf-root <cf.root>]

Options:
  --config <file.toml>      TOML configuration file.
                            Default: <project>/config/pbpb_build_and_fit.toml
  --stage <stage>           Workflow stage to run: all, build-cf, or fit.
                            Default: build-cf
  --input-cf-root <path>    Override the configured CF ROOT input for fit.
  --binary <path>           exp_femto_1d executable path.
                            Default: <project>/bin/exp_femto_1d
  -h, --help                Show this help.

When needed, this helper re-enters the O2Physics ROOT runtime with alienv before
running the ROOT-backed executable.
USAGE
}

# Resolve project-local defaults from the script path so callers can run from any cwd.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
config_path="${project_root}/config/pbpb_build_and_fit.toml"
binary_path="${project_root}/bin/exp_femto_1d"
stage="build-cf"
input_cf_root=""
o2_module="${EXP_FEMTO_1D_O2_MODULE:-O2Physics/latest-master-o2}"
original_arg_count=$#
original_args=("$@")

# Treat the runtime as ready only when the ROOT command-line tools and dynamic
# module path come from a fully entered O2Physics environment.
root_runtime_ready() {
  [[ -n "${ROOTSYS:-}" ]] || return 1
  [[ -n "${ROOT_DYN_PATH:-}" ]] || return 1
  command -v root-config >/dev/null 2>&1 || return 1
  root-config --version >/dev/null 2>&1
}

# Re-exec through alienv for direct shell invocations that have not entered the
# O2Physics runtime yet; a marker prevents recursive retries after alienv fails.
enter_o2_runtime_if_needed() {
  if root_runtime_ready; then
    return
  fi

  if [[ "${EXP_FEMTO_1D_O2_RUNTIME_ENTERED:-0}" == "1" ]]; then
    echo "[error] O2Physics ROOT runtime is not available after alienv setup." >&2
    echo "        Expected ROOTSYS, ROOT_DYN_PATH, and a working root-config." >&2
    exit 1
  fi

  if ! command -v alienv >/dev/null 2>&1; then
    echo "[error] O2Physics ROOT runtime is not available and alienv was not found." >&2
    echo "        Enter O2Physics manually, or install alienv on PATH." >&2
    exit 1
  fi

  export EXP_FEMTO_1D_O2_RUNTIME_ENTERED=1
  exec alienv setenv "${o2_module}" -c bash -lc \
    'exec "$0" "$@"' "${BASH_SOURCE[0]}" "$@"
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

if [[ ! -f "${config_path}" ]]; then
  echo "[error] Config file does not exist: ${config_path}" >&2
  exit 1
fi

if [[ ! -x "${binary_path}" ]]; then
  echo "[error] Executable is missing or not executable: ${binary_path}" >&2
  echo "        Build the project first under the O2Physics ROOT runtime." >&2
  exit 1
fi

if (( original_arg_count > 0 )); then
  enter_o2_runtime_if_needed "${original_args[@]}"
else
  enter_o2_runtime_if_needed
fi

# Keep stage boundaries visible so a long fit stage is not mistaken for a hung shell.
announce_stage() {
  echo "[run] starting ${1}" >&2
}

# Keep the build-cf invocation identical to the public CLI contract.
run_build_cf() {
  announce_stage "build-cf"
  "${binary_path}" build-cf --config "${config_path}"
}

# Keep fit override handling local to the fit stage.
run_fit() {
  announce_stage "fit"
  if [[ -n "${input_cf_root}" ]]; then
    "${binary_path}" fit --config "${config_path}" --input-cf-root "${input_cf_root}"
    return
  fi

  "${binary_path}" fit --config "${config_path}"
}

# Execute the requested stage sequence and let the CLI print its compact summaries.
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
