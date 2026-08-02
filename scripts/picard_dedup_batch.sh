#!/bin/bash
# ATACseq2tracks v3.2.0 - Picard deduplication batch
# PATCHED: samtools quickcheck guard before addreplacerg;
#          tolerant wait loop so one failed RG job doesn't kill the whole batch
# Usage: bash scripts/picard_dedup_batch.sh [max_jobs]
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
if [[ ${#pids[@]} -ge $MAX_JOBS ]]; then
    local first="${pids[0]}"
    wait "$first" || { log "FATAL: Picard job PID $first failed"; exit 1; }
    pids=("${pids[@]:1}")
fi
}

# Collect only original BAMs — exclude *_rg.bam and *_dedup.bam
mapfile -t all_bams < <(find "$IN" -maxdepth 1 -name "*.bam" \
! -name "*_rg.bam" ! -name "*_dedup.bam" | sort)

[[ ${#all_bams[@]} -eq 0 ]] && { log "No source BAMs found in $IN"; exit 1; }
log "Found ${#all_bams[@]} BAMs to process"

# Step 1: Add RG tags -> _rg.bam
log "=== Step 1: Adding RG tags ==="
rg_pids=()
failed_rg=()
for bam in "${all_bams[@]}"; do
base=$(basename "$bam" .bam)
out_rg="${IN}/${base}_rg.bam"
[[ -f "$out_rg" ]] && { log "SKIP RG: $base"; continue; }

# PATCH: validate BAM integrity before attempting to add RG tags
if ! samtools quickcheck "$bam" 2>/dev/null; then
    log "SKIP RG (corrupt BAM — quickcheck failed): $base"
    failed_rg+=("$base")
    continue
fi

(
samtools addreplacerg \
-r "ID:RG1 SM:${base} PL:Illumina LB:Library.fastq2tracks" \
-o "$out_rg" "$bam"
) &
rg_pids+=($!)
log "RG: $base"
done

# PATCH: tolerant wait — log failures instead of aborting entire batch
rg_fail=0
for p in "${rg_pids[@]+"${rg_pids[@]}"}"; do
wait "$p" || { log "WARNING: RG tag job PID $p failed"; rg_fail=$((rg_fail+1)); }
done
[[ $rg_fail -gt 0 ]] && { log "FATAL: $rg_fail read-group job(s) failed"; exit 1; }
log "RG tags complete"

# Step 2: Picard dedup on *_rg.bam -> /_dedup.bam
log "=== Step 2: Picard deduplication ==="
for rg_bam in "${IN}"/*_rg.bam; do
[[ ! -f "$rg_bam" ]] && continue
base=$(basename "$rg_bam" _rg.bam)
[[ -f "${OUT}/${base}_dedup.bam" ]] && { log "SKIP dedup: $base"; continue; }
wait_slot

(
java "${PICARD_XMX:--Xmx8g}" -jar "${PICARD_JAR}" MarkDuplicates \
INPUT="$rg_bam" \
OUTPUT="${OUT}/${base}_dedup.bam" \
METRICS_FILE="${LOG_DIR}/${base}_dup_metrics.txt" \
REMOVE_DUPLICATES=true \
OPTICAL_DUPLICATE_PIXEL_DISTANCE="${PICARD_OPTICAL_DISTANCE:-100}" \
TMP_DIR="${PICARD_TMP:-/tmp}" \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT \
> "$LOG_DIR/${base}_picard.log" 2>&1
samtools index "${OUT}/${base}_dedup.bam" >> "$LOG_DIR/${base}_picard.log" 2>&1
) &
pids+=($!)
log "STARTED dedup: $base"
done
dedup_fail=0
for p in "${pids[@]+"${pids[@]}"}"; do
    wait "$p" || { log "ERROR: Picard job PID $p failed"; dedup_fail=$((dedup_fail + 1)); }
done
[[ $dedup_fail -gt 0 ]] && { log "FATAL: $dedup_fail Picard job(s) failed"; exit 1; }
log "=== Picard batch complete ==="
