#!/bin/bash
# fastq2tracks v3.0.4 — Picard deduplication batch
# Usage: bash scripts/picard_dedup_batch.sh <inDir> <outDir> [max_jobs]
set -euo pipefail

_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.conf not found." >&2; exit 1
        }
    fi
}
_load_config

IN="$1"; OUT="$2"; MAX_JOBS="${3:-${THREADS_PARALLEL_JOBS}}"
mkdir -p "$OUT"
LOG_DIR="${OUT}/../logs/picard"; mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=()
        for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}

# Collect only original BAMs — exclude *_rg.bam and *_dedup.bam
mapfile -t all_bams < <(find "$IN" -maxdepth 1 -name "*.bam" \
    ! -name "*_rg.bam" ! -name "*_dedup.bam" | sort)

[[ ${#all_bams[@]} -eq 0 ]] && { log "No source BAMs found in $IN"; exit 1; }
log "Found ${#all_bams[@]} BAMs to process"

# Step 1: Add RG tags -> <base>_rg.bam
log "=== Step 1: Adding RG tags ==="
rg_pids=()
for bam in "${all_bams[@]}"; do
    base=$(basename "$bam" .bam)
    out_rg="${IN}/${base}_rg.bam"
    [[ -f "$out_rg" ]] && { log "SKIP RG: $base"; continue; }
    samtools addreplacerg \
        -r "ID:RG1	SM:${base}	PL:Illumina	LB:Library.fastq2tracks" \
        -o "$out_rg" "$bam" &
    rg_pids+=($!)
    log "RG: $base"
done
for p in "${rg_pids[@]+"${rg_pids[@]}"}"; do wait "$p"; done
log "RG tags complete"

# Step 2: Picard dedup on *_rg.bam -> <OUT>/<base>_dedup.bam
log "=== Step 2: Picard deduplication ==="
for rg_bam in "${IN}"/*_rg.bam; do
    [[ ! -f "$rg_bam" ]] && continue
    base=$(basename "$rg_bam" _rg.bam)
    [[ -f "${OUT}/${base}_dedup.bam" ]] && { log "SKIP dedup: $base"; continue; }
    wait_slot
    (
        java -Xmx8g -jar "${PICARD_JAR}" MarkDuplicates \
            INPUT="$rg_bam" \
            OUTPUT="${OUT}/${base}_dedup.bam" \
            METRICS_FILE="${LOG_DIR}/${base}_dup_metrics.txt" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            > "$LOG_DIR/${base}_picard.log" 2>&1
        samtools index "${OUT}/${base}_dedup.bam" >> "$LOG_DIR/${base}_picard.log" 2>&1
    ) &
    pids+=($!)
    log "STARTED dedup: $base"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Picard batch complete ==="
