#!/usr/bin/env bash
set -u

status=0
PYTHON_BIN="${PYTHON_BIN:-python3}"
for script in atacseq2tracks.sh fastq2tracks.sh scripts/*.sh utilities/*.sh; do
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

if bash tests/test_v3_2_0_hotfix.sh; then
  echo "OK   v3.2.0 hotfix regression checks"
else
  echo "FAIL v3.2.0 hotfix regression checks"
  status=1
fi

if bash tests/test_windows_input_sanitizer.sh; then
  echo "OK   Windows input sanitizer regression checks"
else
  echo "FAIL Windows input sanitizer regression checks"
  status=1
fi

if bash tests/test_track_normalization.sh; then
  echo "OK   track normalization regression checks"
else
  echo "FAIL track normalization regression checks"
  status=1
fi

if bash tests/test_differential_accessibility.sh; then
  echo "OK   differential accessibility regression checks"
else
  echo "FAIL differential accessibility regression checks"
  status=1
fi

if bash tests/test_hg38_ccre_reference.sh; then
  echo "OK   ENCODE4 hg38 cCRE reference utility checks"
else
  echo "FAIL ENCODE4 hg38 cCRE reference utility checks"
  status=1
fi

"$PYTHON_BIN" - <<'PY' || status=1
import ast
from pathlib import Path

for path in sorted(Path("scripts").glob("*.py")):
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    print(f"OK   {path}")
PY

if command -v Rscript >/dev/null 2>&1; then
  for script in scripts/*.R; do
    if Rscript -e 'parse(file=commandArgs(trailingOnly=TRUE)[1])' "$script" >/dev/null; then
      echo "OK   $script"
    else
      echo "FAIL $script"
      status=1
    fi
  done
else
  echo "WARN Rscript unavailable; R syntax checks skipped"
fi
exit "$status"
