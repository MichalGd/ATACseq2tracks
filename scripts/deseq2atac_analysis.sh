#!/usr/bin/env bash
# Independent broad/narrow all-pair DESeq2ATAC analysis wrapper.
# Usage: deseq2atac_analysis.sh <samplesheet.csv> <bam_dir> <peaks_dir> <out_dir>
set -uo pipefail

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
    hg38)
        BLACKLIST="${BLACKLIST_HG38:-}"
        GTF="${GTF_HUMAN:-}"
        CCRE="${CCRE_BED_HG38:-}"
        CCRE_SOURCE="${CCRE_SOURCE_HG38:-ENCODE4_GRCh38}"
        ;;
    mm39)
        BLACKLIST="${BLACKLIST_MM39:-}"
        GTF="${GTF_MOUSE:-}"
        CCRE="${CCRE_BED_MM39:-}"
        CCRE_SOURCE="${CCRE_SOURCE_MM39:-ENCODE3_mm10_liftOver_mm39}"
        ;;
    *) echo "ERROR: unsupported genome for DESeq2ATAC: $GENOME" >&2; exit 1 ;;
esac
[[ -s "$BLACKLIST" ]] || { echo "ERROR: blacklist missing for $GENOME: $BLACKLIST" >&2; exit 1; }

MIN_SUPPORT="${DESEQ2ATAC_MIN_SAMPLES:-2}"
ALPHA="${DESEQ2ATAC_ALPHA:-0.05}"
BLOCK_COLUMN="${DESEQ2ATAC_BLOCK_COLUMN:-}"
REFERENCE="${DESEQ2ATAC_REFERENCE_CONDITION:-}"
CONDITION_ORDER="${DIFFERENTIAL_CONDITION_ORDER:-}"
if [[ -n "$CONDITION_ORDER" ]]; then
    REFERENCE=""
fi
MIN_ABS_LOG2FC="${DIFFERENTIAL_MIN_ABS_LOG2FC:-0}"
ANNOTATE="${RUN_SIMPLE_PEAK_ANNOTATION:-true}"
CCRE_ANNOTATE="${RUN_CCRE_ANNOTATION:-true}"
case "$CCRE_ANNOTATE" in
    true) ;;
    false) CCRE=""; CCRE_SOURCE="" ;;
    *) echo "ERROR: RUN_CCRE_ANNOTATION must be true or false" >&2; exit 1 ;;
esac
if [[ "$ANNOTATE" == "true" && "$CCRE_ANNOTATE" == "true" && ! -s "$CCRE" ]]; then
    echo "ERROR: cCRE annotation is enabled but the $GENOME cCRE BED is missing: $CCRE" >&2
    echo "Set RUN_CCRE_ANNOTATION=false for GTF-only annotation" >&2
    exit 1
fi
PROMOTER_UPSTREAM="${PEAK_ANNOTATION_PROMOTER_UPSTREAM:-2000}"
PROMOTER_DOWNSTREAM="${PEAK_ANNOTATION_PROMOTER_DOWNSTREAM:-500}"

validate_analysis() {
    local peak_type="$1" type_out="$2" status output
    [[ -s "$type_out/deseq2atac_summary.txt" ]] \
        || { echo "ERROR: DESeq2ATAC $peak_type did not produce a summary" >&2; return 1; }
    status="$(awk -F': ' '$1=="Status"{print $2; exit}' "$type_out/deseq2atac_summary.txt")"
    [[ "$status" == "SUCCESS" || "$status" == "SKIPPED" ]] \
        || { echo "ERROR: DESeq2ATAC $peak_type status is ${status:-missing}" >&2; return 1; }

    for output in deseq2atac_consensus_peaks.bed \
        deseq2atac_consensus_peaks_with_support.tsv.gz \
        deseq2atac_raw_counts.tsv.gz \
        deseq2atac_all_sample_metadata.tsv \
        differential_accessibility_condition_eligibility.tsv \
        differential_accessibility_comparisons.tsv \
        deseq2atac_session_info.txt; do
        [[ -s "$type_out/$output" ]] \
            || { echo "ERROR: DESeq2ATAC $peak_type output missing: $type_out/$output" >&2; return 1; }
    done
    for output in "$type_out"/*.tsv.gz "$type_out"/comparisons/*/*.tsv.gz; do
        [[ -e "$output" ]] || continue
        gzip -t "$output" \
            || { echo "ERROR: unreadable compressed DESeq2ATAC output: $output" >&2; return 1; }
    done

    if [[ "$status" == "SUCCESS" ]]; then
        for output in deseq2atac_normalized_counts.tsv.gz deseq2atac_sample_metadata.tsv \
            deseq2atac_all_sample_metadata.tsv deseq2atac_library_summary.tsv \
            deseq2atac_size_factors.tsv deseq2atac_analysis_object.rds; do
            [[ -s "$type_out/$output" ]] \
                || { echo "ERROR: successful DESeq2ATAC $peak_type output missing: $output" >&2; return 1; }
        done
        for stem in library_sizes_and_size_factors sample_correlation sample_distance \
            pca dispersion_estimates; do
            for extension in png pdf; do
                [[ -s "$type_out/plots/${stem}.${extension}" ]] \
                    || { echo "ERROR: DESeq2ATAC $peak_type figure missing: ${stem}.${extension}" >&2; return 1; }
            done
        done
    fi
}

write_failure_summary() {
    local output="$1" peak_type="$2" message="$3"
    local table="$output/differential_accessibility_comparisons.tsv"
    [[ -s "$table" ]] && return
    printf 'module\tpeak_type\tcomparison_id\tnumerator\treference\tnumerator_replicates\treference_replicates\tconsensus_regions\ttested_sites\tsignificant_sites\thigher_in_numerator\thigher_in_reference\talpha\tmin_abs_log2fc\tstatus\tresults_all\tresults_significant\tsummary_file\tmessage\n' > "$table"
    printf 'DESeq2ATAC\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\tFAILED\tNA\tNA\t%s\t%s\n' \
        "$peak_type" "$ALPHA" "$MIN_ABS_LOG2FC" \
        "$output/deseq2atac_summary.txt" "$message" >> "$table"
}

mkdir -p "$OUT_DIR"
PEAK_TYPE_SUMMARY="$OUT_DIR/deseq2atac_peak_type_summary.tsv"
printf 'peak_type\tstatus\tconsensus_regions\tplanned_comparisons\tsuccessful_comparisons\tfailed_comparisons\tsummary_file\n' \
    > "$PEAK_TYPE_SUMMARY"

failures=0
for PEAK_TYPE in broad narrow; do
    TYPE_OUT="$OUT_DIR/$PEAK_TYPE"
    mkdir -p "$TYPE_OUT"
    echo "=== Running DESeq2ATAC with ${PEAK_TYPE} peaks (all eligible condition pairs) ==="
    if ! "$R_CMD" "$SCRIPT_DIR/deseq2atac_analysis.R" \
        "$SAMPLESHEET" "$BAM_DIR" "$PEAKS_DIR" "$TYPE_OUT" \
        "$GENOME" "$BLACKLIST" "$MIN_SUPPORT" "$ALPHA" \
        "$BLOCK_COLUMN" "$REFERENCE" "${MIN_MAPQ:-30}" "$PEAK_TYPE" \
        "$CONDITION_ORDER" "$MIN_ABS_LOG2FC" "$ANNOTATE" "$GTF" "$CCRE" \
        "$CCRE_SOURCE" "$PROMOTER_UPSTREAM" "$PROMOTER_DOWNSTREAM"; then
        echo "ERROR: DESeq2ATAC $PEAK_TYPE failed; the other peak type remains eligible to run" >&2
        write_failure_summary "$TYPE_OUT" "$PEAK_TYPE" \
            "DESeq2ATAC model or export failed; inspect the R error output"
        failures=$((failures + 1))
        printf '%s\tFAILED\tNA\tNA\tNA\tNA\t%s\n' \
            "$PEAK_TYPE" "$TYPE_OUT/deseq2atac_summary.txt" >> "$PEAK_TYPE_SUMMARY"
        continue
    fi
    if ! validate_analysis "$PEAK_TYPE" "$TYPE_OUT"; then
        write_failure_summary "$TYPE_OUT" "$PEAK_TYPE" "DESeq2ATAC output validation failed"
        failures=$((failures + 1))
        printf '%s\tFAILED_VALIDATION\tNA\tNA\tNA\tNA\t%s\n' \
            "$PEAK_TYPE" "$TYPE_OUT/deseq2atac_summary.txt" >> "$PEAK_TYPE_SUMMARY"
        continue
    fi

    status="$(awk -F': ' '$1=="Status"{print $2; exit}' "$TYPE_OUT/deseq2atac_summary.txt")"
    consensus="$(awk -F': ' '$1=="Consensus regions"{print $2; exit}' "$TYPE_OUT/deseq2atac_summary.txt")"
    planned="$(awk -F': ' '$1=="Planned pairwise comparisons"{print $2; exit}' "$TYPE_OUT/deseq2atac_summary.txt")"
    successful="$(awk -F': ' '$1=="Successful comparisons"{print $2; exit}' "$TYPE_OUT/deseq2atac_summary.txt")"
    failed="$(awk -F': ' '$1=="Failed comparisons"{print $2; exit}' "$TYPE_OUT/deseq2atac_summary.txt")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$PEAK_TYPE" "$status" \
        "${consensus:-NA}" "${planned:-0}" "${successful:-0}" "${failed:-0}" \
        "$TYPE_OUT/deseq2atac_summary.txt" >> "$PEAK_TYPE_SUMMARY"
done

(( failures == 0 )) || exit 1
