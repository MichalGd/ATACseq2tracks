#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - strict per-sample and pooled MACS3 batch
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; BAM_DIR="${2:?filtered BAM directory required}"; PEAKS_DIR="${3:?peaks directory required}"
MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"; LOG="${PEAKS_DIR}/macs3_batch.log"; mkdir -p "$PEAKS_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$LOG"; }
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

declare -a pids=() labels=(); failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: MACS3 failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}
for key in "${!ip_bam[@]}"; do
    [[ "${ip_mode[$key],,}" == "none" ]] && continue
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    out="${PEAKS_DIR}/per_replicate/${key}"; mkdir -p "$out"
    bash "$(dirname "$0")/macs2_peaks.sh" "${ip_bam[$key]}" "${ip_ctrl[$key]}" "$out" "${ip_mode[$key]}" "${ip_genome[$key]}" "$key" "${ip_layout[$key]}" >> "$LOG" 2>&1 &
    pids+=("$!"); labels+=("$key")
done
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} per-sample MACS3 job(s) failed"; exit 1; }

for group in "${!group_bams[@]}"; do
    [[ "${group_mode[$group],,}" == "none" ]] && continue
    read -r -a bams <<< "${group_bams[$group]}"; (( ${#bams[@]} > 0 )) || continue
    group_name="${group//__/_}_pooled"; out="${PEAKS_DIR}/pooled/${group}"; mkdir -p "$out"
    if (( ${#bams[@]} == 1 )); then pooled_bam="${bams[0]}"
    else pooled_bam="${out}/${group_name}.bam"; samtools merge -@ "${THREADS_SAMTOOLS:-2}" -f "$pooled_bam" "${bams[@]}"; samtools index -@ "${THREADS_SAMTOOLS:-2}" "$pooled_bam"; fi
    # Standard ATAC-seq normally has no input control. A pooled control is deliberately not inferred here.
    bash "$(dirname "$0")/macs2_peaks.sh" "$pooled_bam" none "$out" "${group_mode[$group]}" "${group_genome[$group]}" "$group_name" "${group_layout[$group]}" >> "$LOG" 2>&1
done
log "MACS3 peak calling complete"
