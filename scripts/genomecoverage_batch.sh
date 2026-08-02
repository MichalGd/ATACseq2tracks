#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - strict RPM bigWig batch
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; BAM_DIR="${2:?filtered BAM directory required}"; OUT_DIR="${3:?output directory required}"
MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"; LOG="${OUT_DIR}/genomecoverage_batch.log"; mkdir -p "$OUT_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$LOG"; }
declare -a pids=() labels=(); declare -A seen_keys=(); failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: RPM bigWig failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}
while IFS=',' read -r sample_id fq1 fq2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control rest; do
    [[ "$sample_id" == "sample_id" ]] && continue
    sample_id="${sample_id//\"/}"; genome="${genome//\"/}"; replicate="${replicate//\"/}"
    key="${sample_id}_bioR${replicate}"; [[ -n "${seen_keys[$key]+x}" ]] && continue; seen_keys["$key"]=1
    bam="${BAM_DIR}/${key}_dedup_blFilt.bam"; [[ -s "$bam" ]] || { log "ERROR: filtered BAM not found: $bam"; exit 1; }
    [[ -s "${OUT_DIR}/${key}_dedup_blFilt_RPM.bw" ]] && { log "SKIP: $key"; continue; }
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    bash "$(dirname "$0")/genomecoverage_single.sh" "$bam" "$genome" "$OUT_DIR" >> "$LOG" 2>&1 &
    pids+=("$!"); labels+=("$key")
done < "$SAMPLESHEET"
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} RPM bigWig job(s) failed"; exit 1; }
log "RPM bigWig batch complete"
