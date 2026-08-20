#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - strict usable-read and blacklist-filter batch
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; IN_DIR="${2:?deduplicated BAM directory required}"; OUT_DIR="${3:?output directory required}"
MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"; LOG_DIR="${OUT_DIR}/filter_logs"; mkdir -p "$OUT_DIR" "$LOG_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "${LOG_DIR}/main.log"; }
declare -a pids=() labels=(); declare -A seen_keys=(); failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: filtering failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}
while IFS=',' read -r sample_id fq1 fq2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control control_id macs2_mode blacklist rest; do
    [[ "$sample_id" == "sample_id" ]] && continue
    sample_id="${sample_id//\"/}"; layout="${layout//\"/}"; genome="${genome//\"/}"; replicate="${replicate//\"/}"; blacklist="${blacklist//\"/}"
    key="${sample_id}_bioR${replicate}"; [[ -n "${seen_keys[$key]+x}" ]] && continue; seen_keys["$key"]=1
    input_bam="${IN_DIR}/${key}_dedup.bam"; [[ -s "$input_bam" ]] || { log "ERROR: deduplicated BAM not found: $input_bam"; exit 1; }
    [[ -s "$blacklist" ]] || { log "ERROR: blacklist not found for $key: $blacklist"; exit 1; }
    output_bam="${OUT_DIR}/${key}_dedup_blFilt.bam"
    if [[ -s "$output_bam" ]] && samtools quickcheck "$output_bam" 2>/dev/null; then log "SKIP: $key"; continue; fi
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    bash "$(dirname "$0")/blacklist_filter.sh" "$input_bam" "$blacklist" "$OUT_DIR" "$layout" "$genome" >"${LOG_DIR}/${key}.log" 2>&1 &
    pids+=("$!"); labels+=("$key")
done < "$SAMPLESHEET"
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} filtering job(s) failed"; exit 1; }
log "Usable-read and blacklist filtering complete"
