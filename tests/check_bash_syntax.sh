#!/usr/bin/env bash
set -u

status=0
for script in scripts/*.sh; do
  name=$(basename "$script")
  if [[ "$name" == "readme.2.1.sh" ]]; then
    echo "SKIP $script - change-log note, not runnable shell"
    continue
  fi
  if bash -n "$script"; then
    echo "OK   $script"
  else
    echo "FAIL $script"
    status=1
  fi
done
exit "$status"
