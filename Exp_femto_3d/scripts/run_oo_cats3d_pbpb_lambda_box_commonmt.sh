#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
run_tag="oo_cats3d_pbpb_lambda_box_commonmt_20261023"
config_path="${project_root}/config/oo_cats3d_6phi_pbpb_lambda_box_commonmt.toml"
log_directory="${project_root}/logs/cats3d/${run_tag}"
stage="all"

if [[ "${1:-}" == '--help' || "${1:-}" == '-h' ]]; then
  cat <<'USAGE'
Usage: run_oo_cats3d_pbpb_lambda_box_commonmt.sh [--dry-run] [workflow arguments...]

Build a new box-normalized OO CF, then fit CATS3D with the prepared PbPb purity tables.
Uses the three common mT intervals and isolated CF/data/report paths in the TOML.
Output is relayed live through the O2Physics guarded executor; logs use a dedicated directory.
The prepared purity ROOT is consumed as-is; this launcher does not regenerate it.

--dry-run must be first; it previews the command without creating directories or starting ROOT.
Workflow arguments such as --stage fit, --model diag and --config FILE override the defaults.
USAGE
  printf '\nConfiguration: %s\nDefault stage: %s\nLog directory: %s\n' \
    "${config_path}" "${stage}" "${log_directory}"
  exit 0
fi

live_options=(--log-directory "${log_directory}")
if [[ "${1:-}" == '--dry-run' ]]; then
  live_options+=(--dry-run)
  shift
else
  # cf_root_name is absolute and outside output_directory. ROOT needs its parent
  # directory to exist before build-cf; these paths match this launcher's TOML.
  mkdir -p "${project_root}/cf/${run_tag}" \
    "${project_root}/data/${run_tag}" "${project_root}/res/${run_tag}"
fi

exec bash "${script_dir}/run_cats3d_with_progress.sh" "${live_options[@]}" -- \
  --config "${config_path}" --stage "${stage}" "$@"
