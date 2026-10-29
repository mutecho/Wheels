#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
# Defaults for the OO cuttest scan. Both paths can be overridden below.
scan_root="/Users/allenzhou/ALICE/alidata/femtoep_res/OO/lambda_scan_oo_loosecut05thres_6phi_cats3d.root"
output_directory="${project_root}/res/cats3d_scan_report_$(date +%Y%m%d_%H%M%S)"
dry_run=false
summary_only=false

usage() {
  cat <<'USAGE'
Usage: run_cats3d_scan_report.sh [--scan-root FILE] [--output-directory NEW_DIRECTORY] [--summary-only] [--dry-run]

Generate PDF/ROOT reports, per-slice projection images, CSV and an offline HTML
index from a completed CATS3D scan. Uses saved parameters and the recorded native
kernel; does not rerun minimization or modify formal fit/scan outputs.

The scan must contain input provenance, and its CF and compatible kernel must be
available. Output directories must be new. ROOT work runs through the guarded
O2Physics executor. On a 16 GiB machine keep this serial (one kernel in memory).

--summary-only writes the familiar fit_report ROOT layout (FitCatalog,
source_parameters, R2_vs_phi, eps_vs_mt) from saved parameters. It needs no CF or
kernel file, and skips source folding, PNG/PDF/HTML output and slice fitting.
USAGE
}
while (($#)); do
  case "$1" in
    --scan-root|--output-directory)
      (($# >= 2)) || { echo "Missing value after $1" >&2; exit 2; }
      if [[ "$1" == --scan-root ]]; then scan_root="$2"; else output_directory="$2"; fi
      shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --summary-only) summary_only=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -f "${scan_root}" ]] || { echo "Scan ROOT not found: ${scan_root}" >&2; exit 2; }
[[ ! -e "${output_directory}" ]] || { echo "Use a new output directory: ${output_directory}" >&2; exit 2; }
scan_root="$(cd "$(dirname "${scan_root}")" && pwd)/$(basename "${scan_root}")"
if [[ "${output_directory}" != /* ]]; then output_directory="${PWD}/${output_directory}"; fi

printf -v configure_command '%q ' cmake -S "${project_root}" -B "${project_root}/build"
printf -v build_command '%q ' cmake --build "${project_root}/build" --target cats3d_scan_report --parallel 2
printf -v report_command '%q ' "${project_root}/bin/cats3d_scan_report" --scan-root "${scan_root}" --output-directory "${output_directory}"
if [[ "${summary_only}" == true ]]; then report_command+=" --summary-only"; fi
printf -v html_command '%q ' env -u PYTHONHOME -u PYTHONPATH python3 "${script_dir}/render_cats3d_scan_report.py" "${output_directory}"
command="set -e
${configure_command}
${build_command}
${report_command}"
if [[ "${dry_run}" == true ]]; then
  printf '%s\n' "${command}"
  if [[ "${summary_only}" != true ]]; then printf '%s\n' "${html_command}"; fi
  exit 0
fi

printf 'Scan input: %s\nReport output: %s\n' "${scan_root}" "${output_directory}"
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd "${project_root}" --command "${command}"
# O2's Python module variables can conflict with the ordinary Python executable.
# Rendering HTML needs only the standard library and runs in the parent runtime.
if [[ "${summary_only}" != true ]]; then
  env -u PYTHONHOME -u PYTHONPATH python3 "${script_dir}/render_cats3d_scan_report.py" "${output_directory}"
fi
