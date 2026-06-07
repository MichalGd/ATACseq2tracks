#!/bin/bash
# ATACseq2tracks v3.0.4 — Bowtie2 alignment batch (SE + PE, samplesheet-driven)
# PATCHED: added seen_keys guard to prevent duplicate alignment of tech-rep rows
set -euo pipefail

_load_config() {
local _c="${F2T_CONFIG:-}"
if [[ -n "$_c" && -f "$_c" ]]; then source "$_c"
else
local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
_c="${_d}/../config/config.conf"
[[ -f "$_c" ]] && source "$_c" || { echo "ERROR: config.conf not found" >&2; exit 1; }
fi
}

_load_config

SAMPLESHEET="$1"; TRIM_DIR="$2"; OUT_DIR="$3"
mkdir -p "$OUT_DIR"
LOG_DIR="${OUT_DIR}/bowtie2_logs"; mkdir -p "$LOG_DIR"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"; }

declare -a pids=()
wait_slot() {
while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
local new=()
for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
pids=("${new[@]+"${new[@]}"}"); sleep 2
done
}

log "=== Bowtie2 alignment batch ==="

# PATCH: track which KEYs have already been dispatched
declare -A seen_keys=()

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment celltype rep tech_rep is_ctr rest; do
[[ "$sid" == "sample_id" ]] && continue
sid="${sid//\"/}"; layout="${layout//\"/}"; genome="${genome//\"/}"
rep="${rep//\"/}"; tech_rep="${tech_rep//\"/}"
KEY="${sid}_bioR${rep}"

# PATCH: skip if this KEY was already dispatched (duplicate tech-rep row)
if [[ -n "${seen_keys[$KEY]+_}" ]]; then
    log "SKIP DUPLICATE KEY: $KEY (tech_rep=$tech_rep already handled)"
    continue
fi
seen_keys[$KEY]=1

BAM_OUT="${OUT_DIR}/${KEY}.bam"
[[ -f "$BAM_OUT" ]] && { log "SKIP $KEY (BAM exists)"; continue; }

case "$genome" in
hg38) INDEX="$INDEX_HG38" ;;
mm39) INDEX="$INDEX_MM39" ;;
*) log "WARN: unknown genome '$genome' for $KEY, skipping"; continue ;;
esac

wait_slot

JOB_LOG="${LOG_DIR}/${KEY}.log"
(
if [[ "$layout" == "PE" ]]; then
R1="${TRIM_DIR}/${KEY}_1_val_1.fq.gz"
R2="${TRIM_DIR}/${KEY}_2_val_2.fq.gz"
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
echo "ERROR: trimmed PE files not found for $KEY: $R1 / $R2" | tee -a "$JOB_LOG" >&2
exit 1
fi
bowtie2 -x "$INDEX" -1 "$R1" -2 "$R2" \
-p "${THREADS_ALIGN}" --no-mixed --no-discordant --dovetail \
2>>"$JOB_LOG" | samtools sort -@ "${THREADS_SAMTOOLS}" -o "$BAM_OUT"
else
R1="${TRIM_DIR}/${KEY}_trimmed.fq.gz"
if [[ ! -f "$R1" ]]; then
echo "ERROR: trimmed SE file not found for $KEY: $R1" | tee -a "$JOB_LOG" >&2
exit 1
fi
bowtie2 -x "$INDEX" -U "$R1" \
-p "${THREADS_ALIGN}" \
2>>"$JOB_LOG" | samtools sort -@ "${THREADS_SAMTOOLS}" -o "$BAM_OUT"
fi
samtools index "$BAM_OUT" >> "$JOB_LOG" 2>&1
) &
pids+=($!)
log "STARTED: $KEY [$layout] genome=$genome"
done < "$SAMPLESHEET"

for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== Bowtie2 alignment batch complete ==="
