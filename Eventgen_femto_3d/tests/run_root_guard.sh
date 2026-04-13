#!/bin/sh

probe="$1"
shift

"$probe" >/dev/null 2>&1
probe_status=$?
if [ "$probe_status" -ne 0 ]; then
  echo "Skipping ROOT-dependent test because the ROOT runtime probe failed with exit code ${probe_status}."
  exit 77
fi

exec "$@"
