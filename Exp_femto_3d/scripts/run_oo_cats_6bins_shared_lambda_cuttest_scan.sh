#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
config_path="${project_root}/config/oo_cats_6bins_shared_lambda_cuttest_scan.toml"
input_cf_root="/Users/allenzhou/ALICE/alidata/femtoep_res/OO/EP_dependence_CF_3d_25ae_pass2_loosecut05thres_6phi.root"
log_directory="${project_root}/logs/cats3d"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'USAGE'
Usage: run_oo_cats_6bins_shared_lambda_cuttest_scan.sh [--dry-run] [fit arguments...]

Scan the final-lambda coordinates in the TOML with CATS3D, using existing cuttest CF.
The dedicated scan_only configuration preserves formal fit ROOT/TSV/report files.
Edit coordinates and fit_selection in config/oo_cats_6bins_shared_lambda_cuttest_scan.toml.
Native kernel input and memory budgets are configured in its [fit.cats3d] section.
Output and slice-level progress are displayed live; logs are saved in logs/cats3d/.

--dry-run must be first; it prints the command without starting ROOT or fitting.
Additional fit arguments, such as --model diag or --input-cf-root PATH, are forwarded.
USAGE
  exit 0
fi

dry_run=false
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=true
  shift
fi

live_options=(--log-directory "${log_directory}")
if [[ "${dry_run}" == true ]]; then
  live_options+=(--dry-run)
fi

exec bash "${script_dir}/run_cats3d_with_progress.sh" "${live_options[@]}" -- \
  --config "${config_path}" --stage fit --input-cf-root "${input_cf_root}" "$@"
