#!/usr/bin/env bash
# ATAC-specific post-alignment QC: TSS enrichment and fragment periodicity.
set -euo pipefail

SAMPLESHEET="${1:?Usage: ataqv_qc_batch.sh <samplesheet> <bam_dir> <peaks_dir> <out_dir>}"
BAM_DIR="${2:?filtered BAM directory required}"
PEAKS_DIR="${3:?peaks directory required}"
OUT_DIR="${4:?output directory required}"

if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
    # shellcheck disable=SC1090
    source "${F2T_CONFIG}"
else
    SCRIPT_DIR_PRE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    DEFAULT_CONFIG="${SCRIPT_DIR_PRE}/../config/config.conf"
    [[ -f "$DEFAULT_CONFIG" ]] || { echo "ERROR: config.conf not found" >&2; exit 1; }
    # shellcheck disable=SC1090
    source "$DEFAULT_CONFIG"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARALLEL_HELPERS="${SCRIPT_DIR}/parallel_job_helpers.sh"
[[ -f "$PARALLEL_HELPERS" ]] || { echo "ERROR: parallel helper not found: $PARALLEL_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$PARALLEL_HELPERS"
ATAQV_JOBS="${ATAQV_PARALLEL_JOBS:-4}"
ATAQV_THREADS="${THREADS_ATAQV:-8}"
parallel_require_positive_integer ATAQV_PARALLEL_JOBS "$ATAQV_JOBS"
parallel_require_positive_integer THREADS_ATAQV "$ATAQV_THREADS"
TABLE_DIR="${OUT_DIR}/tables"
PLOT_DIR="${OUT_DIR}/plots"
METRICS_DIR="${OUT_DIR}/ataqv_metrics"
LOG_DIR="${OUT_DIR}/logs"
VIEWER_DIR="${OUT_DIR}/ataqv_viewer"
REFERENCE_DIR="${OUT_DIR}/reference"
mkdir -p "$TABLE_DIR" "$PLOT_DIR" "$METRICS_DIR" "$LOG_DIR" "$REFERENCE_DIR"

command -v ataqv >/dev/null 2>&1 || { echo "ERROR: ataqv is required" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is required" >&2; exit 1; }
python3 -c 'import matplotlib' >/dev/null 2>&1 || { echo "ERROR: Python matplotlib is required" >&2; exit 1; }

GENOME="$(tail -n +2 "$SAMPLESHEET" | awk -F',' 'NR==1 {gsub(/"/,"",$5); print $5}')"
case "${GENOME,,}" in
    hg38)
        ORGANISM="human"
        GTF="${GTF_HUMAN:-}"
        CONFIGURED_TSS="${TSS_BED_HG38:-}"
        BLACKLIST="${BLACKLIST_HG38:-}"
        AUTOSOME_MAX=22
        ;;
    mm39)
        ORGANISM="mouse"
        GTF="${GTF_MOUSE:-}"
        CONFIGURED_TSS="${TSS_BED_MM39:-}"
        BLACKLIST="${BLACKLIST_MM39:-}"
        AUTOSOME_MAX=19
        ;;
    *) echo "ERROR: unsupported genome for ATAC QC: $GENOME" >&2; exit 1 ;;
esac

TSS_BED="$CONFIGURED_TSS"
if [[ -z "$TSS_BED" ]]; then
    [[ -f "$GTF" ]] || { echo "ERROR: GTF is required to generate the TSS BED: $GTF" >&2; exit 1; }
    TSS_BED="${REFERENCE_DIR}/${GENOME}.tss.bed"
    python3 "${SCRIPT_DIR}/prepare_tss_bed.py" "$GTF" "$TSS_BED"
fi
[[ -s "$TSS_BED" ]] || { echo "ERROR: TSS BED is empty or missing: $TSS_BED" >&2; exit 1; }

AUTOSOMES="${REFERENCE_DIR}/${GENOME}.autosomal_references.txt"
FIRST_BAM="$(find "$BAM_DIR" -maxdepth 1 -name '*_dedup_blFilt.bam' -print -quit)"
[[ -n "$FIRST_BAM" ]] || { echo "ERROR: no filtered BAM found in $BAM_DIR" >&2; exit 1; }
[[ -f "${FIRST_BAM}.bai" ]] || samtools index -@ "${THREADS_ATAQV:-8}" "$FIRST_BAM"
samtools idxstats "$FIRST_BAM" | awk -v maximum="$AUTOSOME_MAX" '
    $1 != "*" { original=$1; name=$1; sub(/^chr/, "", name); if (name ~ /^[0-9]+$/ && name+0 >= 1 && name+0 <= maximum) print original }
' > "$AUTOSOMES"
[[ -s "$AUTOSOMES" ]] || { echo "ERROR: could not identify autosomes from $FIRST_BAM" >&2; exit 1; }

SELECTED_METRICS="${TABLE_DIR}/ataqv_selected_metrics.tsv"
PERIODICITY_METRICS="${TABLE_DIR}/nucleosome_periodicity_metrics.tsv"
TIMING_TSV="${TABLE_DIR}/ataqv_job_timing.tsv"
printf 'sample_id\tmetric_path\tvalue\n' > "$SELECTED_METRICS"
printf 'sample_id\tmetric\tvalue\n' > "$PERIODICITY_METRICS"
printf 'scope\tlabel\tstart_epoch\tend_epoch\telapsed_seconds\tparallel_jobs\tthreads_per_job\tstatus\n' > "$TIMING_TSV"

declare -A SEEN
declare -a ATAQV_KEYS ATAQV_LAYOUTS ATAQV_MODES ATAQV_BAMS JSON_FILES
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment \
    cell_type rep tech_rep is_ctrl ctrl_id macs2_mode blacklist rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; rep="${rep//\"/}"; layout="${layout//\"/}"
    assay="${assay//\"/}"; is_ctrl="${is_ctrl//\"/}"; macs2_mode="${macs2_mode//\"/}"
    [[ "${is_ctrl,,}" == "true" || "$is_ctrl" == "1" || "$is_ctrl" == "yes" ]] && continue
    [[ "${assay,,}" =~ ^atac ]] || continue
    KEY="${sid}_bioR${rep}"
    [[ -n "${SEEN[$KEY]+x}" ]] && continue
    SEEN["$KEY"]=1

    BAM="${BAM_DIR}/${KEY}_dedup_blFilt.bam"
    [[ -f "$BAM" ]] || { echo "ERROR: filtered BAM not found: $BAM" >&2; exit 1; }
    ATAQV_KEYS+=("$KEY")
    ATAQV_LAYOUTS+=("${layout^^}")
    ATAQV_MODES+=("${macs2_mode,,}")
    ATAQV_BAMS+=("$BAM")
done < <(tail -n +2 "$SAMPLESHEET")

if (( ${#ATAQV_KEYS[@]} == 0 )); then
    echo "[ATAC-QC] No non-control ATAC-seq samples found; nothing to do"
    exit 0
fi

WORK_DIR="$(mktemp -d "${OUT_DIR}/.ataqv-workers.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

run_ataqv_sample_impl() {
    local i="$1" key bam layout mode peak json log sample_selected
    local selected_rows periodicity_rows
    local -a command
    key="${ATAQV_KEYS[$i]}"; bam="${ATAQV_BAMS[$i]}"
    layout="${ATAQV_LAYOUTS[$i]}"; mode="${ATAQV_MODES[$i]}"
    selected_rows="${WORK_DIR}/${i}.selected.rows"
    periodicity_rows="${WORK_DIR}/${i}.periodicity.rows"
    : > "$selected_rows"; : > "$periodicity_rows"

    [[ -f "${bam}.bai" ]] || samtools index -@ "$ATAQV_THREADS" "$bam" || return 1
    peak="${PEAKS_DIR}/per_replicate/${key}/narrow/${key}_peaks.narrowPeak"
    if [[ "$mode" == "broad" || ! -s "$peak" ]]; then
        peak="${PEAKS_DIR}/per_replicate/${key}/broad/${key}_peaks.broadPeak"
    fi

    json="${METRICS_DIR}/${key}.ataqv.json.gz"
    log="${LOG_DIR}/${key}.ataqv.log"
    command=(ataqv --threads "$ATAQV_THREADS" --name "$key"
        --metrics-file "$json" --tss-file "$TSS_BED"
        --tss-extension "${ATAQV_TSS_EXTENSION:-1000}"
        --autosomal-reference-file "$AUTOSOMES" --ignore-read-groups)
    [[ -s "$peak" ]] && command+=(--peak-file "$peak")
    [[ -s "$BLACKLIST" ]] && command+=(--excluded-region-file "$BLACKLIST")
    command+=("$ORGANISM" "$bam")

    "${command[@]}" > "$log" 2>&1 || return 1
    [[ -s "$json" ]] || { echo "ERROR: ataqv did not create metrics for $key" >&2; return 1; }
    sample_selected="${TABLE_DIR}/${key}.ataqv_selected_metrics.tsv"
    python3 "${SCRIPT_DIR}/extract_ataqv_metrics.py" "$json" "$key" "$sample_selected" || return 1
    tail -n +2 "$sample_selected" > "$selected_rows"

    if [[ "$layout" == "PE" ]]; then
        python3 "${SCRIPT_DIR}/plot_fragment_periodicity.py" \
            --bam "$bam" --sample "$key" --out-dir "$PLOT_DIR" \
            --max-fragment-length "${FRAGMENT_PLOT_MAX_BP:-1000}" || return 1
        tail -n +2 "${PLOT_DIR}/${key}.nucleosome_periodicity_metrics.tsv" > "$periodicity_rows"
    else
        printf '%s\tperiodicity_status\tNA_not_applicable_to_single_end\n' "$key" > "$periodicity_rows"
    fi
}

run_ataqv_sample() {
    local i="$1" key start end status=0
    key="${ATAQV_KEYS[$i]}"; start="$(date +%s)"
    run_ataqv_sample_impl "$i" > "${LOG_DIR}/${key}.worker.log" 2>&1 || status=$?
    end="$(date +%s)"
    parallel_write_timing_row "${WORK_DIR}/${i}.timing" ataqv_sample "$key" \
        "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" \
        "$ATAQV_JOBS" "$ATAQV_THREADS"
    return "$status"
}

echo "[ATAC-QC] Running ${#ATAQV_KEYS[@]} samples with up to $ATAQV_JOBS concurrent jobs ($ATAQV_THREADS threads/job)"
parallel_pool_init "$ATAQV_JOBS"
for i in "${!ATAQV_KEYS[@]}"; do
    echo "[ATAC-QC] queued: ${ATAQV_KEYS[$i]}"
    parallel_pool_submit "${ATAQV_KEYS[$i]}" run_ataqv_sample "$i"
done
if ! parallel_pool_wait_all; then
    echo "ERROR: ataqv sample job(s) failed: $(parallel_failed_labels_csv)" >&2
    exit 1
fi

for i in "${!ATAQV_KEYS[@]}"; do
    key="${ATAQV_KEYS[$i]}"
    JSON_FILES+=("${METRICS_DIR}/${key}.ataqv.json.gz")
    cat "${WORK_DIR}/${i}.selected.rows" >> "$SELECTED_METRICS"
    cat "${WORK_DIR}/${i}.periodicity.rows" >> "$PERIODICITY_METRICS"
    cat "${WORK_DIR}/${i}.timing" >> "$TIMING_TSV"
    echo "[ATAC-QC] complete: $key"
done

if [[ "${GENERATE_ATAQV_VIEWER:-true}" == "true" ]]; then
    if command -v mkarv >/dev/null 2>&1; then
        mkarv -f -r calculate -p calculate "$VIEWER_DIR" "${JSON_FILES[@]}" \
            > "${LOG_DIR}/mkarv.log" 2>&1 \
            || echo "WARNING: mkarv viewer generation failed; JSON and static QC outputs remain valid" >&2
    else
        echo "WARNING: mkarv is unavailable; skipping the optional interactive viewer" >&2
    fi
fi

echo "[ATAC-QC] Complete"
echo "  TSS/ataqv metrics : $SELECTED_METRICS"
echo "  Periodicity       : $PERIODICITY_METRICS"
echo "  Job timing        : $TIMING_TSV"
echo "  Static plots      : $PLOT_DIR"
echo "  Compressed JSON   : $METRICS_DIR"
