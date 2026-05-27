#!/bin/bash
# fastq2tracks v3.0.2 — Blacklist filtering batch
# Usage: bash scripts/blacklist_filter_batch.sh <samplesheet.csv> <inDir> <outDir>
set -euo pipefail
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to fastq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

SAMPLESHEET="$1"; IN_DIR="$2"; OUT_DIR="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"
mkdir -p "$OUT_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S"); LOG_DIR="${OUT_DIR}/blacklist_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }
declare -a pids=()
wait_for_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}
log "=== Blacklist filtering batch ==="
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r sample_id fastq_1 fastq_2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control control_id macs2_mode blacklist rest; do
    sample_id="${sample_id//"/}"; blacklist="${blacklist//"/}"
    [[ -z "$blacklist" ]] && { log "SKIP $sample_id — no blacklist defined"; continue; }
    IN_BAM="${IN_DIR}/${sample_id}_dedup.bam"
    [[ ! -f "$IN_BAM" ]] && { log "SKIP $sample_id — dedup BAM not found"; continue; }
    OUT_BAM="${OUT_DIR}/${sample_id}_dedup_blFilt.bam"
    [[ -f "$OUT_BAM" ]] && { log "SKIP $sample_id (filtered BAM exists)"; continue; }
    wait_for_slot
    (bash "$(dirname "$0")/blacklist_filter.sh" "$IN_BAM" "$blacklist" "$OUT_DIR"         > "$LOG_DIR/${sample_id}.log" 2>&1) &
    pids+=($!); log "STARTED: $sample_id"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Blacklist filtering complete ==="
