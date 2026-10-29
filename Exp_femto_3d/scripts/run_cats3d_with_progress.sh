#!/usr/bin/env bash
set -euo pipefail

# Keep ROOT environment validation in the shared executor. A macOS PTY makes
# the existing automatic progress renderer active; a FIFO forwards its output
# while the executor keeps its own environment/status record separately.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"
executor="${CODEX_HOME:-${HOME}/.codex}/skills/cern_root/o2physics-root/scripts/run_root_command.sh"
log_directory="${project_root}/logs/cats3d"
dry_run=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --log-directory)
      [[ $# -ge 2 && -n "$2" ]] || { echo '[error] --log-directory requires a path.' >&2; exit 2; }
      log_directory="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=true
      shift
      ;;
    --)
      shift
      break
      ;;
    *)
      echo '[error] Expected --log-directory PATH, --dry-run, or -- before workflow arguments.' >&2
      exit 2
      ;;
  esac
done
workflow_arguments=("$@")

# Help must remain lightweight and must not start a ROOT environment or create logs.
for argument in "${workflow_arguments[@]}"; do
  if [[ "${argument}" == '--help' || "${argument}" == '-h' ]]; then
    exec bash "${script_dir}/run_exp_femto_3d.sh" "${workflow_arguments[@]}"
  fi
done

[[ -f "${executor}" ]] || { echo "[error] O2Physics guarded executor not found: ${executor}" >&2; exit 1; }
[[ "$(uname -s)" == 'Darwin' && -x /usr/bin/script ]] || {
  echo '[error] This live-output helper requires the macOS /usr/bin/script PTY utility.' >&2
  exit 1
}
if [[ "${log_directory}" != /* ]]; then
  log_directory="${PWD}/${log_directory}"
fi
run_name="cats3d_$(date +%Y%m%d_%H%M%S)_$$"
log_file="${log_directory}/${run_name}.log"
runtime_log="${log_directory}/${run_name}.runtime.log"
pty_arguments=(/usr/bin/script -q -F -a "${log_file}"
  bash "${script_dir}/run_exp_femto_3d.sh" "${workflow_arguments[@]}")
printf -v pty_command '%q ' "${pty_arguments[@]}"

if [[ "${dry_run}" == true ]]; then
  printf '[dry-run] Live run log: %s\n' "${log_file}"
  printf '[dry-run] Runtime/status log: %s\n' "${runtime_log}"
  printf '[dry-run] Guarded command (live relay is connected at launch):\n'
  printf '%q ' bash "${executor}" --cwd "${project_root}" --command "${pty_command}"
  printf '\n'
  exit 0
fi

mkdir -p "${log_directory}"
relay_directory="$(mktemp -d "${TMPDIR:-/tmp}/cats3d-live.XXXXXX")"
relay_fifo="${relay_directory}/output.pipe"
relay_pid=''
cleanup() {
  exec 3>&-
  if [[ -n "${relay_pid}" ]]; then
    kill "${relay_pid}" 2>/dev/null || true
    wait "${relay_pid}" 2>/dev/null || true
  fi
  rm -f "${relay_fifo}"
  rmdir "${relay_directory}" 2>/dev/null || true
}
trap cleanup EXIT
mkfifo "${relay_fifo}"

printf '[run] Command: ' | tee -a "${log_file}"
printf '%q ' bash "${script_dir}/run_exp_femto_3d.sh" "${workflow_arguments[@]}" | tee -a "${log_file}"
printf '\n[run] Live run log: %s\n[run] Runtime/status log: %s\n' "${log_file}" "${runtime_log}"

# Hold the FIFO open until the executor returns, including environment-entry
# failures before its command starts. Close this extra descriptor in both
# children so the reader receives EOF and drains all output at completion.
exec 3<>"${relay_fifo}"
cat "${relay_fifo}" 3>&- &
relay_pid=$!
printf -v redirect ' > %q 2>&1' "${relay_fifo}"
exit_code=0
if bash "${executor}" --cwd "${project_root}" --command "${pty_command}${redirect}" \
    3>&- >"${runtime_log}" 2>&1; then
  exit_code=0
else
  exit_code=$?
fi
exec 3>&-
wait "${relay_pid}" || true
relay_pid=''

runtime_status="$(sed -n 's/^STATUS: //p' "${runtime_log}" | tail -n 1)"
printf '\n[run] Finished: exit=%s, runtime=%s\n' "${exit_code}" "${runtime_status:-unavailable}" | tee -a "${log_file}"
if [[ "${exit_code}" -ne 0 ]]; then
  tail -n 50 "${runtime_log}" >&2
fi
exit "${exit_code}"
