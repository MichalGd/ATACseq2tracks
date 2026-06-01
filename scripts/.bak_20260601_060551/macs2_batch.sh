#!/bin/bash
# fastq2tracks v3.0.4 — MACS2 batch (per-replicate + pooled; always both peak types)
# Usage: bash scripts/macs2_batch.sh <samplesheet.csv> <bamDir> <peaksDir>
set -euo pipefail

_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.conf not found. Export F2T_CONFIG or pass --config to fastq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

SAMPLESHEET="$1"; BAM_DIR="$2"; PEAKS_DIR="$3"
mkdir -p "$PEAKS_DIR"
LOG="${PEAKS_DIR}/macs2_batch_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }

# Pre-build lookup: sample_id -> replicate (for control BAM resolution)
declare -A sid_rep
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment \
    cell_type rep tech_rep is_ctr ctrl_id macs2_mode blacklist rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; rep="${rep//\"/}"
    sid_rep["$sid"]="$rep"
done < "$SAMPLESHEET"

declare -A ip_bams ip_ctrls ip_modes ip_genomes
declare -A grp_bams grp_ctrl grp_mode grp_genome

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment \
    cell_type rep tech_rep is_ctr ctrl_id macs2_mode blacklist rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; is_ctr="${is_ctr//\"/}"; ctrl_id="${ctrl_id//\"/}"
    rep="${rep//\"/}"; genome="${genome//\"/}"; macs2_mode="${macs2_mode//\"/}"
    factor="${factor//\"/}"; condition="${condition//\"/}"
    treatment="${treatment//\"/}"; cell_type="${cell_type//\"/}"

    [[ "${is_ctr,,}" == "true" || "$is_ctr" == "1" ]] && continue

    BAM="${BAM_DIR}/${sid}_bioR${rep}_dedup_blFilt.bam"
    [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${sid}_bioR${rep}_dedup.bam"
    [[ ! -f "$BAM" ]] && { log "WARN: BAM not found for $sid"; continue; }

    CTRL_BAM="none"
    if [[ -n "$ctrl_id" ]]; then
        ctrl_rep="${sid_rep[$ctrl_id]:-1}"
        CB="${BAM_DIR}/${ctrl_id}_bioR${ctrl_rep}_dedup_blFilt.bam"
        [[ ! -f "$CB" ]] && CB="${BAM_DIR}/${ctrl_id}_bioR${ctrl_rep}_dedup.bam"
        [[ -f "$CB" ]] && CTRL_BAM="$CB"
    fi

    ip_bams["$sid"]="$BAM"; ip_ctrls["$sid"]="$CTRL_BAM"
    ip_modes["$sid"]="$macs2_mode"; ip_genomes["$sid"]="$genome"
    GRP="${factor}__${condition}__${treatment}__${cell_type}__${genome}"
    grp_bams["$GRP"]+="$BAM "; grp_ctrl["$GRP"]="$CTRL_BAM"
    grp_mode["$GRP"]="$macs2_mode"; grp_genome["$GRP"]="$genome"
done < "$SAMPLESHEET"

log "=== Per-replicate MACS2 ==="
for sid in "${!ip_bams[@]}"; do
    [[ "${ip_modes[$sid],,}" == "none" ]] && { log "SKIP $sid (macs2_mode=none)"; continue; }
    OUT="${PEAKS_DIR}/per_replicate/${sid}"; mkdir -p "$OUT"
    (bash "$(dirname "$0")/macs2_peaks.sh" \
        "${ip_bams[$sid]}" "${ip_ctrls[$sid]}" "$OUT" \
        "both" "${ip_genomes[$sid]}" "$sid" >> "$LOG" 2>&1) &
done
wait

log "=== Pooled-replicate MACS2 ==="
for GRP in "${!grp_bams[@]}"; do
    BAMS=(${grp_bams[$GRP]})
    [[ ${#BAMS[@]} -eq 0 ]] && continue
    [[ "${grp_mode[$GRP],,}" == "none" ]] && continue
    OUT="${PEAKS_DIR}/pooled/${GRP}"; mkdir -p "$OUT"
    GRP_NAME="${GRP//__/_}_pooled"
    if [[ ${#BAMS[@]} -eq 1 ]]; then
        MERGED_BAM="${BAMS[0]}"
    else
        MERGED_BAM="${PEAKS_DIR}/pooled/${GRP_NAME}_tmp.bam"
        samtools merge -@ "${THREADS_SAMTOOLS}" -f "$MERGED_BAM" "${BAMS[@]}"
        samtools index -@ "${THREADS_SAMTOOLS}" "$MERGED_BAM"
    fi
    (bash "$(dirname "$0")/macs2_peaks.sh" \
        "$MERGED_BAM" "${grp_ctrl[$GRP]}" "$OUT" \
        "both" "${grp_genome[$GRP]}" "$GRP_NAME" >> "$LOG" 2>&1) &
done
wait
log "=== MACS2 batch complete ==="
