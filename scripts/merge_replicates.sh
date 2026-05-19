#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Merge biological/technical replicates
# Usage: bash scripts/merge_replicates.sh <samplesheet.csv> <filteredBamDir> <mergedBamDir> <genome>
#
# Merges BAMs for all samples sharing (factor, condition, treatment, cell_type) 
# grouped across replicates. Technical replicates (same bio replicate, multiple 
# tech_replicate values) are merged first.
# Only generates normalised bedGraph + bigWig for merged BAMs.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

SAMPLESHEET="$1"; IN_DIR="$2"; OUT_DIR="$3"; GENOME="${4:-hg38}"
mkdir -p "$OUT_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG="${OUT_DIR}/merge_replicates_${TIMESTAMP}.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }

declare -A group_bams   # group_key -> space-separated BAM paths

# --- PHASE 1: collect BAMs per group (excluding controls) ---
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctr ctrl_id rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; genome="${genome//\"/}"; factor="${factor//\"/}"
    condition="${condition//\"/}"; treatment="${treatment//\"/}"
    cell_type="${cell_type//\"/}"; is_ctr="${is_ctr//\"/}"
    [[ "${is_ctr,,}" == "true" || "${is_ctr}" == "1" ]] && continue
    [[ "$genome" != "$GENOME" ]] && continue

    BAM="${IN_DIR}/${sid}_dedup_blFilt.bam"
    [[ ! -f "$BAM" ]] && BAM="${IN_DIR}/${sid}_dedup.bam"
    [[ ! -f "$BAM" ]] && { log "WARN: BAM not found for $sid, skipping"; continue; }

    KEY="${factor}__${condition}__${treatment}__${cell_type}"
    group_bams["$KEY"]+="$BAM "
done < "$SAMPLESHEET"

# --- PHASE 2: merge each group ---
for KEY in "${!group_bams[@]}"; do
    BAMS=(${group_bams[$KEY]})
    [[ ${#BAMS[@]} -lt 2 ]] && { log "SKIP $KEY — only 1 sample, no merge needed"; continue; }
    MERGED_NAME="${KEY//__/_}_merged"
    MERGED_BAM="${OUT_DIR}/${MERGED_NAME}.bam"
    log "Merging ${#BAMS[@]} BAMs -> $MERGED_BAM"
    samtools merge -@ "${THREADS_SAMTOOLS}" -f "$MERGED_BAM" "${BAMS[@]}"
    samtools index -@ "${THREADS_SAMTOOLS}" "$MERGED_BAM"
    bash "$(dirname "$0")/genomecoverage_single.sh" "$MERGED_BAM" "$GENOME" "$OUT_DIR" >> "$LOG" 2>&1
    log "Merged coverage done: $MERGED_NAME"
done
log "=== Replicate merging complete ==="
