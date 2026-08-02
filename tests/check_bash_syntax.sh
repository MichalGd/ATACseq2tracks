#!/usr/bin/env bash
set -u

status=0
for script in atacseq2tracks.sh fastq2tracks.sh scripts/*.sh; do
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

if ! cmp -s atacseq2tracks.sh fastq2tracks.sh; then
  echo "FAIL entry points differ: atacseq2tracks.sh and fastq2tracks.sh"
  status=1
else
  echo "OK   entry points are identical"
fi

python3 - <<'PY' || status=1
import ast
from pathlib import Path

for path in sorted(Path("scripts").glob("*.py")):
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    print(f"OK   {path}")
PY
exit "$status"
