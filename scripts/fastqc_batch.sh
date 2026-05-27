#!/bin/bash
# fastq2tracks v3.0.2 — FastQC batch
# Usage: bash scripts/fastqc_batch.sh <inputDir> <outputDir> [max_jobs]
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

IN="$1"; OUT="$2"; MAX_PARAM="${3:-${THREADS_PARALLEL_JOBS}}"
MAX_JOBS=$((2 * MAX_PARAM)); THREADS="${THREADS_FASTQC}"
mkdir -p "${OUT}/fastQC"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S"); LOG_DIR="${OUT}/fastqc_batch_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR/individual_jobs"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }
declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}
mapfile -t all_files < <(find "$IN" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort)
log "=== FastQC batch: ${#all_files[@]} files ==="
for fq in "${all_files[@]}"; do
    base=$(basename "$fq" .fq.gz); base=${base%.fastq.gz}
    [[ -d "${OUT}/fastQC/${base}_fastqc" ]] && { log "SKIP $base"; continue; }
    wait_slot
    (fastqc --outdir "${OUT}/fastQC/" --format fastq --threads "$THREADS" "$fq"         > "$LOG_DIR/individual_jobs/${base}.log" 2>&1) &
    pids+=($!); log "STARTED: $base"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== FastQC batch complete ==="
