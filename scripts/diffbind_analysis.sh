#!/usr/bin/env bash
# Run independent broad/narrow DiffBind models with all eligible condition pairs.
# Usage: bash scripts/diffbind_analysis.sh <diffbind_dir> <out_dir>
set -uo pipefail

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

R_CMD="${R_BIN:-Rscript}"
DIFFBIND_SUMMITS="${DIFFBIND_SUMMITS:-100}"
DIFFBIND_ALPHA="${DIFFBIND_ALPHA:-0.05}"
MIN_ABS_LOG2FC="${DIFFERENTIAL_MIN_ABS_LOG2FC:-0}"
CONDITION_ORDER="${DIFFERENTIAL_CONDITION_ORDER:-${DESEQ2ATAC_REFERENCE_CONDITION:-}}"
ANNOTATE="${RUN_SIMPLE_PEAK_ANNOTATION:-true}"
CCRE_ANNOTATE="${RUN_CCRE_ANNOTATION:-true}"
PROMOTER_UPSTREAM="${PEAK_ANNOTATION_PROMOTER_UPSTREAM:-2000}"
PROMOTER_DOWNSTREAM="${PEAK_ANNOTATION_PROMOTER_DOWNSTREAM:-500}"

[[ "$DIFFBIND_SUMMITS" =~ ^[0-9]+$ ]] \
    || { echo "ERROR: DIFFBIND_SUMMITS must be a non-negative integer" >&2; exit 1; }
command -v "$R_CMD" >/dev/null 2>&1 \
    || { echo "ERROR: R command not found: $R_CMD" >&2; exit 1; }

mkdir -p "$OUT_DIR"
write_failure_summary() {
    local output="$1" peak_type="$2" message="$3"
    local table="$output/differential_accessibility_comparisons.tsv"
    [[ -s "$table" ]] && return
    printf 'module\tpeak_type\tcomparison_id\tnumerator\treference\tnumerator_replicates\treference_replicates\tconsensus_regions\ttested_sites\tsignificant_sites\thigher_in_numerator\thigher_in_reference\talpha\tmin_abs_log2fc\tstatus\tresults_all\tresults_significant\tsummary_file\tmessage\n' > "$table"
    printf 'DiffBind\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\tFAILED\tNA\tNA\t%s\t%s\n' \
        "$peak_type" "$DIFFBIND_ALPHA" "$MIN_ABS_LOG2FC" \
        "$output/diffbind_summary.txt" "$message" >> "$table"
}
found=0
runnable=0
failures=0
shopt -s nullglob
for SS in "$DIFFBIND_DIR"/diffbind_samplesheet_*_*.csv; do
    found=$((found + 1))
    if (( $(wc -l < "$SS") <= 1 )); then
        echo "WARNING: skipping empty DiffBind sample sheet: $SS" >&2
        continue
    fi
    runnable=$((runnable + 1))
    BASE="$(basename "$SS" .csv)"
    GENOME="$(sed -E 's/^diffbind_samplesheet_(hg38|mm39)_.*/\1/' <<< "$BASE")"
    PEAK_TYPE="${BASE##*_}"
    case "$GENOME" in
        hg38)
            GTF="${GTF_HUMAN:-}"
            BLACKLIST="${BLACKLIST_HG38:-}"
            CCRE="${CCRE_BED_HG38:-}"
            CCRE_SOURCE="${CCRE_SOURCE_HG38:-ENCODE4_GRCh38}"
            ;;
        mm39)
            GTF="${GTF_MOUSE:-}"
            BLACKLIST="${BLACKLIST_MM39:-}"
            CCRE="${CCRE_BED_MM39:-}"
            CCRE_SOURCE="${CCRE_SOURCE_MM39:-ENCODE3_mm10_liftOver_mm39}"
            ;;
        *)
            echo "ERROR: cannot infer genome from DiffBind sheet: $SS" >&2
            failures=$((failures + 1))
            continue
            ;;
    esac
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
    SAMPLE_OUT="$OUT_DIR/$BASE"
    mkdir -p "$SAMPLE_OUT"
    echo "=== Running DiffBind for $BASE (all eligible condition pairs) ==="
    if ! "$R_CMD" "$SCRIPT_DIR/diffbind_analysis.R" \
        "$SS" "$SAMPLE_OUT" "$DIFFBIND_SUMMITS" "$DIFFBIND_ALPHA" \
        "$MIN_ABS_LOG2FC" "$CONDITION_ORDER" "$GENOME" "$PEAK_TYPE" \
        "$ANNOTATE" "$GTF" "$CCRE" "$CCRE_SOURCE" \
        "$PROMOTER_UPSTREAM" "$PROMOTER_DOWNSTREAM" "$BLACKLIST"; then
        echo "ERROR: DiffBind failed for $BASE; other peak types remain eligible to run" >&2
        write_failure_summary "$SAMPLE_OUT" "$PEAK_TYPE" "DiffBind model or export failed; inspect diffbind_log.txt"
        failures=$((failures + 1))
    fi
done
shopt -u nullglob

(( found > 0 )) || { echo "ERROR: no DiffBind sample sheets found in $DIFFBIND_DIR" >&2; exit 1; }
(( runnable > 0 )) || { echo "ERROR: all DiffBind sample sheets are empty" >&2; exit 1; }
(( failures == 0 )) || exit 1
