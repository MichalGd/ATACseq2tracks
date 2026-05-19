#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Genome coverage batch (samplesheet-driven)
# Usage: bash scripts/genomecoverage_batch.sh <samplesheet.csv> <filteredBamDir> <outDir>
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

SAMPLESHEET="$1"; BAM_DIR="$2"; OUT_DIR="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"
mkdir -p "$OUT_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="${OUT_DIR}/genomecov_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}

log "=== Genome coverage batch ==="
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r sid fq1 fq2 layout genome rest; do
    sid="${sid//\"/}"; genome="${genome//\"/}"
    BAM="${BAM_DIR}/${sid}_dedup_blFilt.bam"
    [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${sid}_dedup.bam"
    [[ ! -f "$BAM" ]] && { log "SKIP $sid — BAM not found"; continue; }
    [[ -f "${OUT_DIR}/${sid}_dedup_blFilt_Snorm.bw" ]] && { log "SKIP $sid (bw exists)"; continue; }
    [[ -f "${OUT_DIR}/${sid}_dedup_Snorm.bw" ]] && { log "SKIP $sid (bw exists)"; continue; }
    wait_slot
    (bash "$(dirname "$0")/genomecoverage_single.sh" "$BAM" "$genome" "$OUT_DIR" \
         > "$LOG_DIR/${sid}.log" 2>&1) &
    pids+=($!); log "STARTED: $sid [$genome]"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Genome coverage batch complete ==="
