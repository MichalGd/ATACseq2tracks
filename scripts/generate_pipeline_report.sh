#!/bin/bash
# ATACseq2tracks v3.2.0 — Unified MultiQC pipeline report wrapper
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

# Keep all independent differential-accessibility summaries beside the
# MultiQC report. MultiQC does not natively parse these result formats.
DA_SUMMARY="${REPORT_DIR}/differential_accessibility_summary.tsv"
printf 'module\tanalysis\tall_tested_sites\tsignificant_sites\tsummary_file\n' > "$DA_SUMMARY"
while IFS= read -r -d '' summary; do
    analysis="$(basename "$(dirname "$summary")")"
    all_sites="$(awk -F': ' '$1=="All tested sites"{print $2; exit}' "$summary")"
    significant="$(awk -F': ' '$1=="FDR <= 0.05 sites"{print $2; exit}' "$summary")"
    printf 'DiffBind\t%s\t%s\t%s\t%s\n' \
        "$analysis" "${all_sites:-NA}" "${significant:-NA}" "$summary" >> "$DA_SUMMARY"
done < <(find "${OUT_DIR}/diffbind_results" -type f -name 'diffbind_summary.txt' -print0 2>/dev/null)

for peak_type in broad narrow; do
    DESEQ2ATAC_SUMMARY="${OUT_DIR}/deseq2atac/${peak_type}/deseq2atac_summary.txt"
    if [[ -s "$DESEQ2ATAC_SUMMARY" ]]; then
        all_sites="$(awk -F': ' '$1=="Nonzero tested regions"{print $2; exit}' "$DESEQ2ATAC_SUMMARY")"
        significant="$(awk -F': ' '$1=="Significant regions"{print $2; exit}' "$DESEQ2ATAC_SUMMARY")"
        printf 'DESeq2ATAC\t%s_consensus\t%s\t%s\t%s\n' "$peak_type" \
            "${all_sites:-NA}" "${significant:-NA}" "$DESEQ2ATAC_SUMMARY" >> "$DA_SUMMARY"
    fi
done

multiqc "$OUT_DIR" -o "$REPORT_DIR" -n "fastq2tracks_unified_$(date +%Y%m%d)" \
    --data-format tsv --export -f
echo "Unified MultiQC report in: $REPORT_DIR"
echo "Differential-accessibility summary: $DA_SUMMARY"
