#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Bowtie2 batch (samplesheet-driven, SE+PE, human+mouse)
# Usage: bash scripts/bowtie2_batch.sh <samplesheet.csv> <trimmedFastqDir> <outBamDir>
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

SAMPLESHEET="$1"
TRIM_DIR="$2"
OUT_DIR="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"

mkdir -p "$OUT_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="${OUT_DIR}/bowtie2_batch_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR/individual_jobs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a job_pids=()
declare -a job_names=()

wait_for_slot() {
    while [[ ${#job_pids[@]} -ge $MAX_JOBS ]]; do
        local new_pids=() new_names=()
        for i in "${!job_pids[@]}"; do
            if kill -0 "${job_pids[$i]}" 2>/dev/null; then
                new_pids+=("${job_pids[$i]}")
                new_names+=("${job_names[$i]}")
            else
                log "COMPLETED: ${job_names[$i]}"
            fi
        done
        job_pids=("${new_pids[@]+"${new_pids[@]}"}")
        job_names=("${new_names[@]+"${new_names[@]}"}")
        [[ ${#job_pids[@]} -ge $MAX_JOBS ]] && sleep 3
    done
}

log "=== Bowtie2 Batch === Samplesheet: $SAMPLESHEET"

tail -n +2 "$SAMPLESHEET" | while IFS=',' read -r sample_id fastq_1 fastq_2 layout genome assay factor condition treatment cell_type replicate tech_rep is_control control_id macs2_mode blacklist chipqc_anno output_prefix rest; do
    sample_id="${sample_id//\"/}"
    layout="${layout//\"/}"
    genome="${genome//\"/}"

    # resolve index
    if [[ "$genome" == "hg38" ]]; then INDEX="$INDEX_HG38"
    elif [[ "$genome" == "mm39" ]]; then INDEX="$INDEX_MM39"
    else log "SKIP $sample_id — unknown genome $genome"; continue; fi

    # resolve trimmed R1
    R1_file="${TRIM_DIR}/${sample_id}_1_val_1.fq.gz"
    [[ ! -f "$R1_file" ]] && R1_file="${TRIM_DIR}/${sample_id}_R1_001_val_1.fq.gz"
    [[ ! -f "$R1_file" ]] && R1_file="${TRIM_DIR}/${sample_id}_trimmed.fq.gz"
    [[ ! -f "$R1_file" ]] && { log "SKIP $sample_id — trimmed R1 not found"; continue; }

    OUT_BAM="${OUT_DIR}/${sample_id}.sorted_stChr.bam"
    [[ -f "$OUT_BAM" ]] && { log "SKIP $sample_id (BAM exists)"; continue; }

    wait_for_slot

    (
        JOB_LOG="$LOG_DIR/individual_jobs/${sample_id}.log"
        if [[ "$layout" == "PE" ]]; then
            R2_file="${TRIM_DIR}/${sample_id}_2_val_2.fq.gz"
            [[ ! -f "$R2_file" ]] && R2_file="${TRIM_DIR}/${sample_id}_R2_001_val_2.fq.gz"
            bash "$(dirname "$0")/bowtie2_align.sh" "$INDEX" "$R1_file" "$OUT_DIR" "PE" "$R2_file" >"$JOB_LOG" 2>&1
        else
            bash "$(dirname "$0")/bowtie2_align.sh" "$INDEX" "$R1_file" "$OUT_DIR" "SE" >"$JOB_LOG" 2>&1
        fi
    ) &
    job_pids+=($!)
    job_names+=("$sample_id")
    log "STARTED: $sample_id [${layout}] [${genome}] [PID:${job_pids[-1]}]"
done

# wait remaining
for pid in "${job_pids[@]+"${job_pids[@]}"}"; do wait "$pid" || true; done
log "=== Bowtie2 batch complete ==="
