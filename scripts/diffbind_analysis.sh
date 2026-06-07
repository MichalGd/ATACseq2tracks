#!/usr/bin/env bash
# ATACseq2tracks v3.1.x — DiffBind differential accessibility analysis wrapper
# Usage: bash scripts/diffbind_analysis.sh <diffbind_dir> <out_dir>
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIFFBIND_DIR="${1:-}"
OUT_DIR="${2:-}"
if [[ -z "$DIFFBIND_DIR" || -z "$OUT_DIR" ]]; then
    echo "Usage: bash scripts/diffbind_analysis.sh <diffbind_dir> <out_dir>" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"
shopt -s nullglob
for SS in "$DIFFBIND_DIR"/diffbind_samplesheet_*_*.csv; do
    BASE=$(basename "$SS" .csv)
    SAMPLE_OUT="$OUT_DIR/$BASE"
    mkdir -p "$SAMPLE_OUT"
    echo "=== Running DiffBind differential analysis for: $BASE ==="
    Rscript "$SCRIPT_DIR/diffbind_analysis.R" "$SS" "$SAMPLE_OUT"
done
shopt -u nullglob
