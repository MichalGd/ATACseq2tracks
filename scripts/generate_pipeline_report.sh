#!/bin/bash
# ATACseq2tracks v4.0.0 — Unified MultiQC pipeline report wrapper
# Usage: bash scripts/generate_pipeline_report.sh <outDir> [reportDir] [format]
set -euo pipefail
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to atacseq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

OUT_DIR="$1"; REPORT_DIR="${2:-${OUT_DIR}/reports}"; FORMAT="${3:-html}"
mkdir -p "$REPORT_DIR"
SUMMARY_DIR="${OUT_DIR}/reports"
mkdir -p "$SUMMARY_DIR"

# MultiQC does not natively parse these result formats. Build a stable pair-level
# TSV plus a compact standalone HTML table before invoking MultiQC.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 "${SCRIPT_DIR}/summarize_differential_accessibility.py" "$OUT_DIR" "$SUMMARY_DIR"
DA_SUMMARY="${SUMMARY_DIR}/differential_accessibility_summary.tsv"
DA_HTML="${SUMMARY_DIR}/differential_accessibility_summary.html"

multiqc "$OUT_DIR" -o "$REPORT_DIR" -n "fastq2tracks_unified_$(date +%Y%m%d)" \
    --data-format tsv --export -f
echo "Unified MultiQC report in: $REPORT_DIR"
echo "Differential-accessibility summary: $DA_SUMMARY"
echo "Differential-accessibility HTML: $DA_HTML"
