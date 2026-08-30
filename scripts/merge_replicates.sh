#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - merge filtered biological samples for group-level CPM tracks
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; IN_DIR="${2:?filtered BAM directory required}"; OUT_DIR="${3:?output directory required}"; GENOME="${4:?genome required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARALLEL_HELPERS="${SCRIPT_DIR}/parallel_job_helpers.sh"
[[ -f "$PARALLEL_HELPERS" ]] || { echo "ERROR: parallel helper not found: $PARALLEL_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$PARALLEL_HELPERS"
MERGE_JOBS="${MERGE_PARALLEL_JOBS:-2}"
parallel_require_positive_integer MERGE_PARALLEL_JOBS "$MERGE_JOBS"
LOG="${OUT_DIR}/merge_replicates.log"; TIMING_TSV="${OUT_DIR}/merge_job_timing.tsv"; mkdir -p "$OUT_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$LOG"; }
printf 'scope\tlabel\tstart_epoch\tend_epoch\telapsed_seconds\tparallel_jobs\tthreads_per_job\tstatus\n' > "$TIMING_TSV"
declare -A group_bams=() group_layout=() seen=()
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctrl rest; do
    [[ "$sid" == "sample_id" ]] && continue
    for name in sid layout genome factor condition treatment cell_type rep is_ctrl; do value="${!name}"; printf -v "$name" '%s' "${value//\"/}"; done
    [[ "$genome" != "$GENOME" || "${is_ctrl,,}" =~ ^(true|1|yes)$ ]] && continue
    key="${sid}_bioR${rep}"; [[ -n "${seen[$key]+x}" ]] && continue; seen["$key"]=1
    bam="${IN_DIR}/${key}_dedup_blFilt.bam"; [[ -s "$bam" ]] || { log "ERROR: filtered BAM missing: $bam"; exit 1; }
    group="${factor}__${condition}__${treatment}__${cell_type}"
    if [[ -n "${group_layout[$group]:-}" && "${group_layout[$group]}" != "${layout^^}" ]]; then
        log "ERROR: cannot merge PE and SE samples into one signal track: $group"
        exit 1
    fi
    group_layout["$group"]="${layout^^}"
    group_bams["$group"]+="$bam "
done < "$SAMPLESHEET"

run_merge_group() {
    local group="$1" merged_name merged_bam start end status=0 job_log timing_file
    local -a bams
    read -r -a bams <<< "${group_bams[$group]}"
    (( ${#bams[@]} >= 2 )) || return 0
    merged_name="${group//__/_}_merged"; merged_bam="${OUT_DIR}/${merged_name}.bam"
    job_log="${OUT_DIR}/.${merged_name}.merge_job.log"
    timing_file="${OUT_DIR}/.${merged_name}.merge_job_timing.tsv"
    start="$(date +%s)"
    samtools merge -@ "${THREADS_SAMTOOLS:-2}" -f "$merged_bam" "${bams[@]}" \
        > "$job_log" 2>&1 || status=$?
    if (( status == 0 )); then
        samtools index -@ "${THREADS_SAMTOOLS:-2}" "$merged_bam" >> "$job_log" 2>&1 || status=$?
    fi
    if (( status == 0 )); then
        bash "${SCRIPT_DIR}/genomecoverage_single.sh" "$merged_bam" "$GENOME" "$OUT_DIR" \
            "${group_layout[$group]}" >> "$job_log" 2>&1 || status=$?
    fi
    end="$(date +%s)"
    parallel_write_timing_row "$timing_file" merge_and_cpm "$merged_name" \
        "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" \
        "$MERGE_JOBS" "${THREADS_SAMTOOLS:-2}"
    return "$status"
}

mapfile -t GROUP_KEYS < <(printf '%s\n' "${!group_bams[@]}" | LC_ALL=C sort)
parallel_pool_init "$MERGE_JOBS"
for group in "${GROUP_KEYS[@]}"; do
    read -r -a bams <<< "${group_bams[$group]}"
    if (( ${#bams[@]} < 2 )); then
        log "SKIP group with one biological sample: $group"
        continue
    fi
    log "QUEUE replicate merge: $group"
    parallel_pool_submit "$group" run_merge_group "$group"
done
merge_status=0
parallel_pool_wait_all || merge_status=$?
for group in "${GROUP_KEYS[@]}"; do
    merged_name="${group//__/_}_merged"
    job_log="${OUT_DIR}/.${merged_name}.merge_job.log"
    timing_file="${OUT_DIR}/.${merged_name}.merge_job_timing.tsv"
    [[ -f "$job_log" ]] && cat "$job_log" >> "$LOG"
    [[ -f "$timing_file" ]] && cat "$timing_file" >> "$TIMING_TSV"
    rm -f "$job_log" "$timing_file"
done
(( merge_status == 0 )) || { log "FATAL: replicate merge job(s) failed: $(parallel_failed_labels_csv)"; exit 1; }
log "Merge timing: $TIMING_TSV"
log "Group merge and CPM track generation complete"
