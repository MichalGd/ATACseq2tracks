#!/usr/bin/env bash
# ATACseq2tracks v3.1.x — DiffBind differential accessibility analysis wrapper
# Usage: bash scripts/diffbind_analysis.sh <diffbind_dir> <out_dir>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$F2T_CONFIG"
fi
DIFFBIND_DIR="${1:-}"
OUT_DIR="${2:-}"
if [[ -z "$DIFFBIND_DIR" || -z "$OUT_DIR" ]]; then
    echo "Usage: bash scripts/diffbind_analysis.sh <diffbind_dir> <out_dir>" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"
R_CMD="${R_BIN:-Rscript}"
DIFFBIND_SUMMITS="${DIFFBIND_SUMMITS:-100}"
[[ "$DIFFBIND_SUMMITS" =~ ^[0-9]+$ ]] \
    || { echo "ERROR: DIFFBIND_SUMMITS must be a non-negative integer" >&2; exit 1; }
command -v "$R_CMD" >/dev/null 2>&1 || { echo "ERROR: R command not found: $R_CMD" >&2; exit 1; }
shopt -s nullglob
found=0
runnable=0
for SS in "$DIFFBIND_DIR"/diffbind_samplesheet_*_*.csv; do
    found=$((found + 1))
    if (( $(wc -l < "$SS") <= 1 )); then
        echo "WARNING: skipping empty DiffBind sample sheet: $SS" >&2
        continue
    fi
    runnable=$((runnable + 1))
    BASE=$(basename "$SS" .csv)
    SAMPLE_OUT="$OUT_DIR/$BASE"
    mkdir -p "$SAMPLE_OUT"
    echo "=== Running DiffBind differential analysis for: $BASE ==="
    "$R_CMD" "$SCRIPT_DIR/diffbind_analysis.R" "$SS" "$SAMPLE_OUT" "$DIFFBIND_SUMMITS"
done
shopt -u nullglob
(( found > 0 )) || { echo "ERROR: no DiffBind sample sheets found in $DIFFBIND_DIR" >&2; exit 1; }
(( runnable > 0 )) || { echo "ERROR: all DiffBind sample sheets are empty" >&2; exit 1; }
