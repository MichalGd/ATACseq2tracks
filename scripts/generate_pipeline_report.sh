#!/bin/bash
# ATACseq2tracks v4.2.0 — Unified MultiQC pipeline report wrapper
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
SPIKEIN_TABLE="${OUT_DIR}/spikein/tables/spikein_normalization.tsv"
SPIKEIN_WARNINGS="${OUT_DIR}/spikein/tables/spikein_warnings.tsv"
if [[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
    [[ -s "$SPIKEIN_TABLE" && -s "$SPIKEIN_WARNINGS" ]] || {
        echo "ERROR: enabled dm6 spike-in calibration tables are missing" >&2
        exit 1
    }
fi
SPIKEIN_TABLE="${OUT_DIR}/spikein/tables/spikein_normalization.tsv"
SPIKEIN_WARNINGS="${OUT_DIR}/spikein/tables/spikein_warnings.tsv"
if [[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
    [[ -s "$SPIKEIN_TABLE" && -s "$SPIKEIN_WARNINGS" ]] || {
        echo "ERROR: enabled dm6 spike-in calibration tables are missing" >&2
        exit 1
    }
fi

# MultiQC 1.35 can misparse deepTools plotPCA tables (null point names), merge
# repeated plotPCA/plotCorrelation sample IDs and reject deepTools RGB triplets.
# Keep the authoritative deepTools tables and plots untouched, exclude only its
# native MultiQC parser, and embed the already-rendered static QC plots as
# self-contained MultiQC custom-content images.
CUSTOM_DIR="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-multiqc.XXXXXX")"
cleanup_custom_dir() { rm -rf -- "$CUSTOM_DIR"; }
trap cleanup_custom_dir EXIT

QC_PLOTS="${OUT_DIR}/qc_post_alignment/plots"
CUSTOM_IMAGE_COUNT=0
stage_qc_image() {
    local source="$1" destination="$2"
    [[ -s "$source" ]] || return 0
    cp -- "$source" "${CUSTOM_DIR}/${destination}"
    CUSTOM_IMAGE_COUNT=$((CUSTOM_IMAGE_COUNT + 1))
}

stage_qc_image "${QC_PLOTS}/pca_bins.png" \
    "ATACseq2tracks_deepTools_PCA_genome_wide_mqc.png"
stage_qc_image "${QC_PLOTS}/pca_peaks.png" \
    "ATACseq2tracks_deepTools_PCA_consensus_peaks_mqc.png"
stage_qc_image "${QC_PLOTS}/correlation_heatmap_pearson.png" \
    "ATACseq2tracks_deepTools_Pearson_genome_wide_mqc.png"
stage_qc_image "${QC_PLOTS}/correlation_heatmap_spearman.png" \
    "ATACseq2tracks_deepTools_Spearman_genome_wide_mqc.png"
stage_qc_image "${QC_PLOTS}/correlation_heatmap_pearson_peaks.png" \
    "ATACseq2tracks_deepTools_Pearson_consensus_peaks_mqc.png"
stage_qc_image "${QC_PLOTS}/fingerprint.png" \
    "ATACseq2tracks_deepTools_fingerprint_mqc.png"
stage_qc_image "${QC_PLOTS}/heatmap_signal_over_peaks.png" \
    "ATACseq2tracks_deepTools_signal_heatmap_mqc.png"
stage_qc_image "${QC_PLOTS}/profile_signal_over_peaks.png" \
    "ATACseq2tracks_deepTools_signal_profile_mqc.png"

if (( CUSTOM_IMAGE_COUNT == 0 )); then
    echo "WARNING: no standalone deepTools QC plots were available for MultiQC custom content" >&2
else
    echo "MultiQC: staging ${CUSTOM_IMAGE_COUNT} standalone deepTools QC plots as custom content"
fi

REPORT_NAME="fastq2tracks_unified_$(date +%Y%m%d)"
MULTIQC_LOG="${REPORT_DIR}/${REPORT_NAME}.multiqc.log"
if ! multiqc "$OUT_DIR" "$CUSTOM_DIR" -o "$REPORT_DIR" -n "$REPORT_NAME" \
    --exclude deeptools --ignore "*/reports/*" \
    --cl-config 'ignore_images: false' \
    --data-format tsv --export -f 2>&1 | tee "$MULTIQC_LOG"; then
    echo "ERROR: MultiQC command failed; inspect $MULTIQC_LOG" >&2
    exit 1
fi

if grep -Eiq "Oops! The .* MultiQC module broke|ValidationError" "$MULTIQC_LOG"; then
    echo "ERROR: MultiQC reported a module or validation failure; inspect $MULTIQC_LOG" >&2
    exit 1
fi

# MultiQC 1.35 can log non-fatal mqc_colour errors while exporting plots when
# built-in palettes contain comma-separated RGB triplets. MultiQC still writes
# the HTML report and exported images. Accept only that precise known pattern;
# continue to reject any other colour-conversion error.
COLOUR_ERRORS="$(grep -Ei "Error converting colou?r" "$MULTIQC_LOG" || true)"
if [[ -n "$COLOUR_ERRORS" ]]; then
    UNKNOWN_COLOUR_ERRORS="$(printf '%s\n' "$COLOUR_ERRORS" | \
        grep -Eiv "mqc_colour.*Error converting colou?r ['\"]?[0-9]{1,3},[0-9]{1,3},[0-9]{1,3}['\"]? to RGB" || true)"
    if [[ -n "$UNKNOWN_COLOUR_ERRORS" ]]; then
        echo "ERROR: MultiQC reported an unexpected colour-conversion failure; inspect $MULTIQC_LOG" >&2
        exit 1
    fi
    echo "WARNING: MultiQC reported known non-fatal RGB-triplet export messages; exported plots will be validated" >&2
fi

[[ -s "${REPORT_DIR}/${REPORT_NAME}.html" ]] || {
    echo "ERROR: MultiQC did not create a non-empty report: ${REPORT_DIR}/${REPORT_NAME}.html" >&2
    exit 1
}
EXPORTED_PLOTS_DIR="${REPORT_DIR}/${REPORT_NAME}_plots"
if ! find "$EXPORTED_PLOTS_DIR" -type f \
        \( -name '*.png' -o -name '*.svg' -o -name '*.pdf' \) \
        -size +0c -print -quit 2>/dev/null | grep -q .; then
    echo "ERROR: MultiQC did not create non-empty exported report images: $EXPORTED_PLOTS_DIR" >&2
    exit 1
fi
[[ -s "$DA_SUMMARY" && -s "$DA_HTML" ]] || {
    echo "ERROR: differential-accessibility summary output is missing or empty" >&2
    exit 1
}

echo "Unified MultiQC report in: $REPORT_DIR"
echo "Differential-accessibility summary: $DA_SUMMARY"
echo "Differential-accessibility HTML: $DA_HTML"
if [[ -s "$SPIKEIN_TABLE" ]]; then
    echo "Drosophila spike-in normalization: $SPIKEIN_TABLE"
    echo "Drosophila spike-in warnings: $SPIKEIN_WARNINGS"
fi
if [[ -s "$SPIKEIN_TABLE" ]]; then
    echo "Drosophila spike-in normalization: $SPIKEIN_TABLE"
    echo "Drosophila spike-in warnings: $SPIKEIN_WARNINGS"
fi
