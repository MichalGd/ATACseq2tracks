#!/usr/bin/env bash
# ATACseq2tracks v4.0.0 - fail-fast FastQC batch
# Usage: fastqc_batch.sh <input_dir> <output_dir> [parallel_jobs] [directory|samplesheet]
set -euo pipefail

if [[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]]; then
    # shellcheck disable=SC1090
    source "$F2T_CONFIG"
fi

IN="${1:?input directory required}"
OUT="${2:?output directory required}"
MAX_JOBS="${3:-${THREADS_PARALLEL_JOBS:-2}}"
SOURCE_MODE="${4:-directory}"
THREADS="${THREADS_FASTQC:-2}"
QC_DIR="${OUT}/fastQC"
LOG_DIR="${OUT}/fastqc_logs"
mkdir -p "$QC_DIR" "$LOG_DIR"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "${LOG_DIR}/main.log"; }

declare -a all_files=()
if [[ "$SOURCE_MODE" == "samplesheet" ]]; then
    [[ -f "${SAMPLESHEET:-}" ]] || { log "ERROR: SAMPLESHEET is unavailable"; exit 1; }
    SHEET_DIR="$(cd "$(dirname "$SAMPLESHEET")" && pwd)"
    while IFS=',' read -r sid fq1 fq2 rest; do
        [[ "$sid" == "sample_id" ]] && continue
        for fq in "$fq1" "$fq2"; do
            fq="${fq//\"/}"
            [[ -z "$fq" ]] && continue
            if [[ -f "$fq" ]]; then all_files+=("$fq")
            elif [[ -f "${IN}/${fq}" ]]; then all_files+=("${IN}/${fq}")
            elif [[ -f "${SHEET_DIR}/${fq}" ]]; then all_files+=("${SHEET_DIR}/${fq}")
            else log "ERROR: FASTQ not found: $fq"; exit 1
            fi
        done
    done < "$SAMPLESHEET"
else
    mapfile -t all_files < <(find "$IN" -maxdepth 1 -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) | sort)
fi

mapfile -t all_files < <(printf '%s\n' "${all_files[@]}" | awk 'NF && !seen[$0]++')
(( ${#all_files[@]} > 0 )) || { log "ERROR: no FASTQ files selected"; exit 1; }

declare -a pids=() labels=()
failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: FastQC failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}")
    labels=("${labels[@]:1}")
}

for fq in "${all_files[@]}"; do
    base="$(basename "$fq")"; base="${base%.fastq.gz}"; base="${base%.fq.gz}"
    [[ -s "${QC_DIR}/${base}_fastqc.html" ]] && { log "SKIP: $base"; continue; }
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    fastqc --outdir "$QC_DIR" --format fastq --threads "$THREADS" "$fq" >"${LOG_DIR}/${base}.log" 2>&1 &
    pids+=("$!"); labels+=("$base")
done
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} FastQC job(s) failed"; exit 1; }
log "FastQC complete: ${#all_files[@]} files"
