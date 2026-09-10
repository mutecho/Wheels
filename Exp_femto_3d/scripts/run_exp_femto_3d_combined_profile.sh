#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"

# The standard runner owns ROOT runtime entry and CLI validation. This wrapper
# only selects the public combined-profile production configuration.
exec "${script_dir}/run_exp_femto_3d.sh" \
  --stage fit \
  --config "${project_root}/config/oo_build_and_fit_6bins_combined_profile.toml" \
  "$@"
