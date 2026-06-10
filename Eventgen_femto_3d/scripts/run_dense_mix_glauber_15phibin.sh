#!/usr/bin/env bash

set -euo pipefail

# Reject external arguments so this remains a true one-click analysis entry.
if [[ "$#" -ne 0 ]]; then
  echo "Usage: $0" >&2
  echo "This script has all input config paths embedded internally." >&2
  exit 2
fi

# Resolve repository paths from this script location.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
runner="${script_dir}/run_eventgen_femto_3d.sh"
config_path="${project_root}/config/bw_dense_mix_glauber_15phibin copy.toml"
output_dir="${project_root}/res"

# Keep the dense-mix Glauber inputs and outputs configured here.
b1_input="/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b1_dense_mix_glauber.root"
b3_input="/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b3_dense_mix_glauber.root"
b6_input="/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b6_dense_mix_glauber.root"
b8_input="/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b8_dense_mix_glauber.root"
b1_output="${output_dir}/PbPb_b1_dense_mix_glauber_15phibin.root"
b3_output="${output_dir}/PbPb_b3_dense_mix_glauber_15phibin.root"
b6_output="${output_dir}/PbPb_b6_dense_mix_glauber_15phibin.root"
b8_output="${output_dir}/PbPb_b8_dense_mix_glauber_15phibin.root"

# Ensure the embedded runner and TOML config are available before launching ROOT.
if [[ ! -x "${runner}" ]]; then
  echo "Runner not executable: ${runner}" >&2
  exit 1
fi
if [[ ! -f "${config_path}" ]]; then
  echo "Config file not found: ${config_path}" >&2
  exit 1
fi

# Run the b1 analysis using the internally configured input and output paths.
echo "Running dense-mix Glauber b1 15-phi-bin analysis..."
"${runner}" \
  --config "${config_path}" \
  --input-root "${b1_input}" \
  --output-root "${b1_output}"

# Run the b3 analysis by overriding the internally configured input and output.
echo "Running dense-mix Glauber b3 15-phi-bin analysis..."
"${runner}" \
  --config "${config_path}" \
  --input-root "${b3_input}" \
  --output-root "${b3_output}"

# Run the b6 analysis by overriding the internally configured input and output.
echo "Running dense-mix Glauber b6 15-phi-bin analysis..."
"${runner}" \
  --config "${config_path}" \
  --input-root "${b6_input}" \
  --output-root "${b6_output}"

# Run the b8 analysis by overriding the internally configured input and output.
echo "Running dense-mix Glauber b8 15-phi-bin analysis..."
"${runner}" \
  --config "${config_path}" \
  --input-root "${b8_input}" \
  --output-root "${b8_output}"

# Print the generated artifacts for quick confirmation in terminal logs.
echo "Dense-mix Glauber 15-phi-bin analysis completed."
echo "b1 output: ${b1_output}"
echo "b3 output: ${b3_output}"
echo "b6 output: ${b6_output}"
echo "b8 output: ${b8_output}"
