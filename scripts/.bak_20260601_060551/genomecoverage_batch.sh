#!/bin/bash
# fastq2tracks v3.0.2 — Genome coverage batch
# Usage: bash scripts/genomecoverage_batch.sh <samplesheet.csv> <bamDir> <outDir>
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

SAMPLESHEET="$1"; BAM_DIR="$2"; OUT_DIR="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"
mkdir -p "$OUT_DIR"
LOG="${OUT_DIR}/genomecoverage_batch_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }
declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}
log "=== Genome coverage batch ==="
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r sample_id fastq_1 fastq_2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control rest; do
    sample_id="${sample_id//\"/}"; genome="${genome//\"/}"
    BAM="${BAM_DIR}/${sample_id}_bioR${replicate}_dedup_blFilt.bam"
    [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${sample_id}_bioR${replicate}_dedup.bam"
    [[ ! -f "$BAM" ]] && { log "SKIP $sample_id — BAM not found"; continue; }
    [[ -f "${OUT_DIR}/${sample_id}_bioR${replicate}_dedup_blFilt_Snorm.bw" ]] && { log "SKIP $sample_id (bw exists)"; continue; }
    wait_slot
    (bash "$(dirname "$0")/genomecoverage_single.sh" "$BAM" "$genome" "$OUT_DIR" >> "$LOG" 2>&1) &
    pids+=($!); log "STARTED: $sample_id"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Coverage batch complete ==="
