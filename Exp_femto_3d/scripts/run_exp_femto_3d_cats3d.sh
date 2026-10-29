#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Edit this path to select the TOML configuration.
# config_path="${script_dir}/../config/oo_build_and_fit_6bins_cats3d.toml"
config_path="${script_dir}/../config/oo_cats_6bins_shared_lambda_cuttest_scan.toml"
log_directory="${script_dir}/../logs/cats3d"

if [[ "${1:-}" == '--help' || "${1:-}" == '-h' ]]; then
  cat <<'USAGE'
Usage: run_exp_femto_3d_cats3d.sh [--dry-run] [workflow arguments...]

Run the selected CATS3D configuration with live output and slice-level progress.
Edit config_path and log_directory near the top of this script.
--dry-run must be first; it previews the launch without starting ROOT or fitting.
Overrides such as --config FILE, --model diag, --input-cf-root FILE and --binary FILE
are forwarded to run_exp_femto_3d.sh. The selected stage defaults to fit.
USAGE
  printf '\nSelected configuration: %s\nLog directory: %s\n' "${config_path}" "${log_directory}"
  exit 0
fi

live_options=(--log-directory "${log_directory}")
if [[ "${1:-}" == '--dry-run' ]]; then
  live_options+=(--dry-run)
  shift
fi

exec bash "${script_dir}/run_cats3d_with_progress.sh" "${live_options[@]}" -- \
  --config "${config_path}" --stage fit "$@"
