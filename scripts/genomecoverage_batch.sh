#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - strict fragment/read CPM track batch
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
    if ! wait "$pid"; then log "ERROR: CPM track generation failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}
while IFS=',' read -r sample_id fq1 fq2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control rest; do
    [[ "$sample_id" == "sample_id" ]] && continue
    sample_id="${sample_id//\"/}"; genome="${genome//\"/}"; replicate="${replicate//\"/}"; layout="${layout//\"/}"
    key="${sample_id}_bioR${replicate}"; [[ -n "${seen_keys[$key]+x}" ]] && continue; seen_keys["$key"]=1
    bam="${BAM_DIR}/${key}_dedup_blFilt.bam"; [[ -s "$bam" ]] || { log "ERROR: filtered BAM not found: $bam"; exit 1; }
    outputs_complete=true
    [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" && ! -s "${OUT_DIR}/${key}_dedup_blFilt_CPM.bw" ]] && outputs_complete=false
    [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" && ! -s "${OUT_DIR}/${key}_dedup_blFilt_CPM.bedGraph" ]] && outputs_complete=false
    [[ "$outputs_complete" == "true" ]] && { log "SKIP: $key"; continue; }
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    bash "$(dirname "$0")/genomecoverage_single.sh" "$bam" "$genome" "$OUT_DIR" "$layout" >> "$LOG" 2>&1 &
    pids+=("$!"); labels+=("$key")
done < "$SAMPLESHEET"
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} CPM track job(s) failed"; exit 1; }
log "CPM track batch complete"
