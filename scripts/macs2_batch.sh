#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — MACS2 batch (samplesheet-driven; replicates analyzed together)
# Usage: bash scripts/macs2_batch.sh <samplesheet.csv> <filteredBamDir> <peaksOutDir>
#
# IP samples sharing (factor, condition, treatment, cell_type) are pooled per
# condition across replicates for group-level peak calling.
# Additionally calls peaks on each individual replicate BAM.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

SAMPLESHEET="$1"; BAM_DIR="$2"; PEAKS_DIR="$3"
mkdir -p "$PEAKS_DIR"
LOG="${PEAKS_DIR}/macs2_batch_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }

declare -A ip_bams ip_ctrls ip_modes ip_genomes  # keyed by sample_id
declare -A grp_bams grp_ctrl grp_mode grp_genome  # keyed by group

# --- parse samplesheet ---
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctr ctrl_id macs2_mode blacklist rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; is_ctr="${is_ctr//\"/}"; ctrl_id="${ctrl_id//\"/}"
    genome="${genome//\"/}"; macs2_mode="${macs2_mode//\"/}"
    factor="${factor//\"/}"; condition="${condition//\"/}"
    treatment="${treatment//\"/}"; cell_type="${cell_type//\"/}"
    [[ "${is_ctr,,}" == "true" || "$is_ctr" == "1" ]] && continue

    BAM="${BAM_DIR}/${sid}_dedup_blFilt.bam"
    [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${sid}_dedup.bam"
    [[ ! -f "$BAM" ]] && { log "WARN: BAM not found for $sid"; continue; }

    CTRL_BAM="none"
    if [[ -n "$ctrl_id" ]]; then
        CB="${BAM_DIR}/${ctrl_id}_dedup_blFilt.bam"
        [[ ! -f "$CB" ]] && CB="${BAM_DIR}/${ctrl_id}_dedup.bam"
        [[ -f "$CB" ]] && CTRL_BAM="$CB"
    fi

    ip_bams["$sid"]="$BAM"
    ip_ctrls["$sid"]="$CTRL_BAM"
    ip_modes["$sid"]="$macs2_mode"
    ip_genomes["$sid"]="$genome"

    GRP="${factor}__${condition}__${treatment}__${cell_type}__${genome}"
    grp_bams["$GRP"]+="$BAM "
    grp_ctrl["$GRP"]="$CTRL_BAM"
    grp_mode["$GRP"]="$macs2_mode"
    grp_genome["$GRP"]="$genome"
done < "$SAMPLESHEET"

# --- per-replicate peak calling ---
log "=== Per-replicate MACS2 ==="
for sid in "${!ip_bams[@]}"; do
    OUT="${PEAKS_DIR}/per_replicate/${sid}"
    mkdir -p "$OUT"
    (bash "$(dirname "$0")/macs2_peaks.sh" \
        "${ip_bams[$sid]}" "${ip_ctrls[$sid]}" "$OUT" \
        "${ip_modes[$sid]}" "${ip_genomes[$sid]}" "$sid" >> "$LOG" 2>&1) &
done
wait

# --- grouped (pooled replicates) peak calling ---
log "=== Pooled-replicate MACS2 ==="
for GRP in "${!grp_bams[@]}"; do
    BAMS=(${grp_bams[$GRP]})
    [[ ${#BAMS[@]} -eq 0 ]] && continue
    OUT="${PEAKS_DIR}/pooled/${GRP}"
    mkdir -p "$OUT"
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
        "${grp_mode[$GRP]}" "${grp_genome[$GRP]}" "$GRP_NAME" >> "$LOG" 2>&1) &
done
wait
log "=== MACS2 batch complete ==="
