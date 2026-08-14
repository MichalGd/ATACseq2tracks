#!/usr/bin/env bash
# Independent broad- and narrow-consensus DESeq2 differential accessibility wrapper.
# Usage: deseq2atac_analysis.sh <samplesheet.csv> <bam_dir> <peaks_dir> <out_dir>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$F2T_CONFIG"
fi

SAMPLESHEET="${1:?samplesheet.csv required}"
BAM_DIR="${2:?filtered BAM directory required}"
PEAKS_DIR="${3:?peaks directory required}"
OUT_DIR="${4:?output directory required}"
R_CMD="${R_BIN:-Rscript}"

command -v "$R_CMD" >/dev/null 2>&1 \
    || { echo "ERROR: R command not found: $R_CMD" >&2; exit 1; }
[[ -s "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet missing or empty: $SAMPLESHEET" >&2; exit 1; }

GENOME="$(python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    genomes = sorted({row["genome"].strip().lower() for row in csv.DictReader(handle)
                      if row["genome"].strip()})
if len(genomes) != 1:
    raise SystemExit("DESeq2ATAC requires exactly one genome per run")
print(genomes[0])
PY
)"

case "$GENOME" in
    hg38) BLACKLIST="${BLACKLIST_HG38:-}" ;;
    mm39) BLACKLIST="${BLACKLIST_MM39:-}" ;;
    *) echo "ERROR: unsupported genome for DESeq2ATAC: $GENOME" >&2; exit 1 ;;
esac
[[ -s "$BLACKLIST" ]] || { echo "ERROR: blacklist missing for $GENOME: $BLACKLIST" >&2; exit 1; }

MIN_SUPPORT="${DESEQ2ATAC_MIN_SAMPLES:-2}"
ALPHA="${DESEQ2ATAC_ALPHA:-0.05}"
BLOCK_COLUMN="${DESEQ2ATAC_BLOCK_COLUMN:-}"
REFERENCE="${DESEQ2ATAC_REFERENCE_CONDITION:-}"

required_outputs=(
    deseq2atac_consensus_peaks.bed
    deseq2atac_consensus_peaks_with_support.tsv.gz
    deseq2atac_raw_counts.tsv.gz
    deseq2atac_normalized_counts.tsv.gz
    deseq2atac_results_all.tsv.gz
    deseq2atac_results_significant.tsv.gz
    deseq2atac_sample_metadata.tsv
    deseq2atac_library_summary.tsv
    deseq2atac_size_factors.tsv
    deseq2atac_analysis_object.rds
    deseq2atac_session_info.txt
)

validate_analysis() {
    local peak_type="$1" type_out="$2" output stem extension significant_count
    [[ -s "$type_out/deseq2atac_summary.txt" ]] \
        || { echo "ERROR: DESeq2ATAC $peak_type did not produce its completion summary" >&2; return 1; }
    grep -q '^Status: SUCCESS$' "$type_out/deseq2atac_summary.txt" \
        || { echo "ERROR: DESeq2ATAC $peak_type summary does not report success" >&2; return 1; }
    grep -q "^Peak type: ${peak_type}$" "$type_out/deseq2atac_summary.txt" \
        || { echo "ERROR: DESeq2ATAC $peak_type summary has the wrong peak type" >&2; return 1; }

    for output in "${required_outputs[@]}"; do
        [[ -s "$type_out/$output" ]] \
            || { echo "ERROR: DESeq2ATAC $peak_type output missing or empty: $type_out/$output" >&2; return 1; }
    done
    for output in "$type_out"/*.tsv.gz; do
        gzip -t "$output" \
            || { echo "ERROR: unreadable compressed DESeq2ATAC $peak_type output: $output" >&2; return 1; }
    done
    for stem in library_sizes_and_size_factors sample_correlation sample_distance \
        pca dispersion_estimates ma volcano; do
        for extension in png pdf; do
            [[ -s "$type_out/plots/${stem}.${extension}" ]] \
                || { echo "ERROR: DESeq2ATAC $peak_type figure missing or empty: ${stem}.${extension}" >&2; return 1; }
        done
    done

    significant_count="$(awk -F': ' '$1=="Significant regions"{print $2; exit}' \
        "$type_out/deseq2atac_summary.txt")"
    if [[ "${significant_count:-0}" -gt 0 ]]; then
        [[ -s "$type_out/plots/significant_site_overview.png" && \
           -s "$type_out/plots/significant_site_overview.pdf" ]] \
            || { echo "ERROR: DESeq2ATAC $peak_type significant-site overview is missing" >&2; return 1; }
    fi
}

mkdir -p "$OUT_DIR"
if [[ -s "$OUT_DIR/deseq2atac_summary.txt" ]]; then
    echo "WARNING: legacy single-analysis DESeq2ATAC files remain in $OUT_DIR;" \
         "new broad/narrow results use subdirectories and do not read those files." >&2
fi
COMPARISON_SUMMARY="$OUT_DIR/deseq2atac_peak_type_summary.tsv"
printf 'peak_type\tconsensus_regions\tnonzero_tested_regions\tsignificant_regions\tsummary_file\n' \
    > "$COMPARISON_SUMMARY"

for PEAK_TYPE in broad narrow; do
    TYPE_OUT="$OUT_DIR/$PEAK_TYPE"
    mkdir -p "$TYPE_OUT"
    echo "=== Running DESeq2ATAC with ${PEAK_TYPE} peaks ==="
    "$R_CMD" "$SCRIPT_DIR/deseq2atac_analysis.R" \
        "$SAMPLESHEET" "$BAM_DIR" "$PEAKS_DIR" "$TYPE_OUT" \
        "$GENOME" "$BLACKLIST" "$MIN_SUPPORT" "$ALPHA" \
        "$BLOCK_COLUMN" "$REFERENCE" "${MIN_MAPQ:-30}" "$PEAK_TYPE"
    validate_analysis "$PEAK_TYPE" "$TYPE_OUT"

    consensus_count="$(awk -F': ' '$1=="Consensus regions"{print $2; exit}' \
        "$TYPE_OUT/deseq2atac_summary.txt")"
    tested_count="$(awk -F': ' '$1=="Nonzero tested regions"{print $2; exit}' \
        "$TYPE_OUT/deseq2atac_summary.txt")"
    significant_count="$(awk -F': ' '$1=="Significant regions"{print $2; exit}' \
        "$TYPE_OUT/deseq2atac_summary.txt")"
    printf '%s\t%s\t%s\t%s\t%s\n' "$PEAK_TYPE" "${consensus_count:-NA}" \
        "${tested_count:-NA}" "${significant_count:-NA}" \
        "$TYPE_OUT/deseq2atac_summary.txt" >> "$COMPARISON_SUMMARY"
done

[[ -s "$COMPARISON_SUMMARY" ]] \
    || { echo "ERROR: DESeq2ATAC peak-type summary was not created" >&2; exit 1; }
