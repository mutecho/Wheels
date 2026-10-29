#!/usr/bin/env bash
set -euo pipefail
task_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
printf -v compare_command '%q ' bash ./build_and_run.sh "$@"
exec bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd "$task_dir" --command "$compare_command"
