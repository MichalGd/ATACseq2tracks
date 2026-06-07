#!/bin/bash
# ATACseq2tracks v3.0.4 — Genome coverage batch (tech-replicate aware)
# Usage: bash scripts/genomecoverage_batch.sh
set -euo pipefail
_load_config() {
  if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
    source "${F2T_CONFIG}"
  else
    local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    local _c="${_d}/../config/config.conf"
    [[ -f "$_c" ]] && source "$_c" || {
      echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to atacseq2tracks.sh." >&2
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
declare -A seen_keys=()
wait_slot() {
  while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
    local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
    pids=("${new[@]+"${new[@]}"}"); sleep 2
  done
}
log "=== Genome coverage batch (tech-rep aware) ==="
tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r sample_id fastq_1 fastq_2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control rest; do
  sample_id="${sample_id//\"/}"; genome="${genome//\"/}"
  replicate="${replicate//\"/}"
  KEY="${sample_id}_bioR${replicate}"
  [[ -n "${seen_keys[$KEY]+x}" ]] && continue
  seen_keys["$KEY"]=1
  BAM="${BAM_DIR}/${KEY}_dedup_blFilt.bam"
  [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${KEY}_dedup.bam"
  [[ ! -f "$BAM" ]] && { log "SKIP $KEY — BAM not found"; continue; }
  [[ -f "${OUT_DIR}/${KEY}_dedup_blFilt_Snorm.bw" ]] && { log "SKIP $KEY (bw exists)"; continue; }
  wait_slot
  (bash "$(dirname "$0")/genomecoverage_single.sh" "$BAM" "$genome" "$OUT_DIR" >> "$LOG" 2>&1) &
  pids+=($!); log "STARTED: $KEY"
done
FAIL=0
for p in "${pids[@]+"${pids[@]}"}"; do
    wait "$p" || { log "ERROR: coverage job failed (pid $p)"; FAIL=1; }
done
[[ $FAIL -eq 1 ]] && { log "FATAL: one or more coverage jobs failed"; exit 1; }
log "=== Coverage batch complete ==="
