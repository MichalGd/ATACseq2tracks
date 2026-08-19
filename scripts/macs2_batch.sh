#!/usr/bin/env bash
# ATACseq2tracks v4.0.0 - strict per-sample and pooled MACS3 batch
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; BAM_DIR="${2:?filtered BAM directory required}"; PEAKS_DIR="${3:?peaks directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARALLEL_HELPERS="${SCRIPT_DIR}/parallel_job_helpers.sh"
[[ -f "$PARALLEL_HELPERS" ]] || { echo "ERROR: parallel helper not found: $PARALLEL_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$PARALLEL_HELPERS"
MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"
POOLED_JOBS="${POOLED_MACS_PARALLEL_JOBS:-2}"
parallel_require_positive_integer THREADS_PARALLEL_JOBS "$MAX_JOBS"
parallel_require_positive_integer POOLED_MACS_PARALLEL_JOBS "$POOLED_JOBS"
LOG="${PEAKS_DIR}/macs3_batch.log"; TIMING_TSV="${PEAKS_DIR}/macs3_job_timing.tsv"; mkdir -p "$PEAKS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$LOG"; }
printf 'scope\tlabel\tstart_epoch\tend_epoch\telapsed_seconds\tparallel_jobs\tthreads_per_job\tstatus\n' > "$TIMING_TSV"
declare -A ip_bam=() ip_ctrl=() ip_mode=() ip_genome=() ip_layout=() seen=()
declare -A group_bams=() group_mode=() group_genome=() group_layout=()
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctrl ctrl_id peak_mode blacklist rest; do
    [[ "$sid" == "sample_id" ]] && continue
    for name in sid layout genome factor condition treatment cell_type rep is_ctrl ctrl_id peak_mode; do
        value="${!name}"
        printf -v "$name" '%s' "${value//\"/}"
    done
    [[ "${is_ctrl,,}" =~ ^(true|1|yes)$ ]] && continue
    key="${sid}_bioR${rep}"; [[ -n "${seen[$key]+x}" ]] && continue; seen["$key"]=1
    bam="${BAM_DIR}/${key}_dedup_blFilt.bam"; [[ -s "$bam" ]] || { log "ERROR: filtered BAM missing: $bam"; exit 1; }
    control="none"
    if [[ -n "$ctrl_id" ]]; then
        for candidate in "${BAM_DIR}/${ctrl_id}_bioR${rep}_dedup_blFilt.bam" "${BAM_DIR}/${ctrl_id}_bioR1_dedup_blFilt.bam"; do
            [[ -s "$candidate" ]] && { control="$candidate"; break; }
        done
        [[ "$control" == "none" ]] && { log "ERROR: declared control not found for $key: $ctrl_id"; exit 1; }
    fi
    ip_bam["$key"]="$bam"; ip_ctrl["$key"]="$control"; ip_mode["$key"]="$peak_mode"; ip_genome["$key"]="$genome"; ip_layout["$key"]="$layout"
    group="${factor}__${condition}__${treatment}__${cell_type}__${genome}__${layout}"
    group_bams["$group"]+="$bam "
    if [[ -n "${group_mode[$group]:-}" && "${group_mode[$group]}" != "$peak_mode" ]]; then group_mode["$group"]="both"; else group_mode["$group"]="$peak_mode"; fi
    group_genome["$group"]="$genome"; group_layout["$group"]="$layout"
done < "$SAMPLESHEET"

run_replicate_macs() {
    local key="$1" out start end status=0
    out="${PEAKS_DIR}/per_replicate/${key}"; mkdir -p "$out"
    start="$(date +%s)"
    bash "${SCRIPT_DIR}/macs2_peaks.sh" "${ip_bam[$key]}" "${ip_ctrl[$key]}" \
        "$out" "${ip_mode[$key]}" "${ip_genome[$key]}" "$key" "${ip_layout[$key]}" \
        > "${out}/.macs3_job.log" 2>&1 || status=$?
    end="$(date +%s)"
    parallel_write_timing_row "${out}/.macs3_job_timing.tsv" per_replicate "$key" \
        "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" \
        "$MAX_JOBS" 1
    return "$status"
}

mapfile -t IP_KEYS < <(printf '%s\n' "${!ip_bam[@]}" | LC_ALL=C sort)
parallel_pool_init "$MAX_JOBS"
for key in "${IP_KEYS[@]}"; do
    [[ "${ip_mode[$key],,}" == "none" ]] && continue
    log "QUEUE per-replicate MACS3: $key"
    parallel_pool_submit "$key" run_replicate_macs "$key"
done
replicate_status=0
parallel_pool_wait_all || replicate_status=$?
for key in "${IP_KEYS[@]}"; do
    out="${PEAKS_DIR}/per_replicate/${key}"
    [[ -f "${out}/.macs3_job.log" ]] && cat "${out}/.macs3_job.log" >> "$LOG"
    [[ -f "${out}/.macs3_job_timing.tsv" ]] && cat "${out}/.macs3_job_timing.tsv" >> "$TIMING_TSV"
    rm -f "${out}/.macs3_job.log" "${out}/.macs3_job_timing.tsv"
done
(( replicate_status == 0 )) || { log "FATAL: per-sample MACS3 job(s) failed: $(parallel_failed_labels_csv)"; exit 1; }

run_pooled_macs() {
    local group="$1" group_name out pooled_bam start end status=0
    local -a bams
    read -r -a bams <<< "${group_bams[$group]}"; (( ${#bams[@]} > 0 )) || return 0
    group_name="${group//__/_}_pooled"; out="${PEAKS_DIR}/pooled/${group}"; mkdir -p "$out"
    start="$(date +%s)"
    if (( ${#bams[@]} == 1 )); then
        pooled_bam="${bams[0]}"
    else
        pooled_bam="${out}/${group_name}.bam"
        samtools merge -@ "${THREADS_SAMTOOLS:-2}" -f "$pooled_bam" "${bams[@]}" \
            > "${out}/.pooled_macs3_job.log" 2>&1 || status=$?
        if (( status == 0 )); then
            samtools index -@ "${THREADS_SAMTOOLS:-2}" "$pooled_bam" \
                >> "${out}/.pooled_macs3_job.log" 2>&1 || status=$?
        fi
    fi
    # Standard ATAC-seq normally has no input control. A pooled control is deliberately not inferred here.
    if (( status == 0 )); then
        bash "${SCRIPT_DIR}/macs2_peaks.sh" "$pooled_bam" none "$out" \
            "${group_mode[$group]}" "${group_genome[$group]}" "$group_name" "${group_layout[$group]}" \
            >> "${out}/.pooled_macs3_job.log" 2>&1 || status=$?
    fi
    end="$(date +%s)"
    parallel_write_timing_row "${out}/.pooled_macs3_job_timing.tsv" pooled "$group_name" \
        "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" \
        "$POOLED_JOBS" "${THREADS_SAMTOOLS:-2}"
    return "$status"
}

mapfile -t GROUP_KEYS < <(printf '%s\n' "${!group_bams[@]}" | LC_ALL=C sort)
parallel_pool_init "$POOLED_JOBS"
for group in "${GROUP_KEYS[@]}"; do
    [[ "${group_mode[$group],,}" == "none" ]] && continue
    log "QUEUE pooled MACS3: ${group//__/_}_pooled"
    parallel_pool_submit "${group//__/_}_pooled" run_pooled_macs "$group"
done
pooled_status=0
parallel_pool_wait_all || pooled_status=$?
for group in "${GROUP_KEYS[@]}"; do
    out="${PEAKS_DIR}/pooled/${group}"
    [[ -f "${out}/.pooled_macs3_job.log" ]] && cat "${out}/.pooled_macs3_job.log" >> "$LOG"
    [[ -f "${out}/.pooled_macs3_job_timing.tsv" ]] && cat "${out}/.pooled_macs3_job_timing.tsv" >> "$TIMING_TSV"
    rm -f "${out}/.pooled_macs3_job.log" "${out}/.pooled_macs3_job_timing.tsv"
done
(( pooled_status == 0 )) || { log "FATAL: pooled MACS3 job(s) failed: $(parallel_failed_labels_csv)"; exit 1; }
log "MACS3 timing: $TIMING_TSV"
log "MACS3 peak calling complete"
