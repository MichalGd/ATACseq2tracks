#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Picard deduplication batch
# Usage: bash scripts/picard_dedup_batch.sh <bamDir> <outDir> [max_jobs]
# Processes all *H.bam files (after samtools addreplacerg) or *.sorted_stChr.bam
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

IN="$1"; OUT="$2"; MAX_JOBS="${3:-${THREADS_PARALLEL_JOBS}}"
mkdir -p "$OUT"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="${OUT}/picard_batch_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR/individual_jobs"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}

# First add RG tags if missing, then deduplicate
all_bams=("${IN}"/*.sorted_stChr.bam)
[[ ! -e "${all_bams[0]}" ]] && { log "No *.sorted_stChr.bam in $IN"; exit 1; }

log "Adding RG tags..."
rg_pids=()
for bam in "${all_bams[@]}"; do
    base=$(basename "$bam" .bam)
    out_rg="${IN}/${base}H.bam"
    [[ -f "$out_rg" ]] && continue
    samtools addreplacerg \
        -r "@RG\tID:RG1\tSM:${base}\tPL:Illumina\tLB:Library.fastq2tracks" \
        -o "$out_rg" "$bam" &
    rg_pids+=($!)
done
for p in "${rg_pids[@]+"${rg_pids[@]}"}"; do wait "$p"; done
log "RG tags done"

log "=== Picard deduplication ==="
for bam in "${IN}"/*.sorted_stChrH.bam; do
    [[ ! -f "$bam" ]] && continue
    base=$(basename "$bam")
    [[ -f "${OUT}/${base}_dedup.bam" ]] && { log "SKIP $base (dedup exists)"; continue; }
    wait_slot
    (bash "$(dirname "$0")/picard_dedup.sh" "$IN" "$OUT" "$base" \
         > "$LOG_DIR/individual_jobs/${base}.log" 2>&1) &
    pids+=($!); log "STARTED: $base"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Picard batch complete ==="
