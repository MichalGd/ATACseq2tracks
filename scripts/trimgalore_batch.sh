#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — TrimGalore batch (SE + PE, samplesheet-driven)
# Usage: bash scripts/trimgalore_batch.sh <samplesheet.csv> <rawFastqDir> <outDir>
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

SAMPLESHEET="$1"; IN="$2"; OUT="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"
mkdir -p "$OUT"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_DIR="${OUT}/trimgalore_batch_logs_${TIMESTAMP}"
mkdir -p "$LOG_DIR/individual_jobs"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}

# Collect per-sample tech replicates first → merge raw FASTQs if needed
declare -A tech_r1 tech_r2 tech_layout

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctr rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; fq1="${fq1//\"/}"; fq2="${fq2//\"/}"
    layout="${layout//\"/}"; tech_rep="${tech_rep//\"/}"

    # key without tech replicate to group them
    KEY="${sid%_tech*}"
    [[ -n "$tech_rep" ]] && KEY="${KEY}_bioR${rep//\"/}"

    tech_r1["$KEY"]+="${IN}/${fq1} "
    [[ -n "$fq2" && "$layout" == "PE" ]] && tech_r2["$KEY"]+="${IN}/${fq2} "
    tech_layout["$KEY"]="${layout}"
done < "$SAMPLESHEET"

log "=== TrimGalore batch ==="

for KEY in "${!tech_r1[@]}"; do
    R1s=(${tech_r1[$KEY]})
    LAYOUT="${tech_layout[$KEY]}"

    # check if trimmed output already exists
    if [[ "$LAYOUT" == "PE" ]]; then
        [[ -f "${OUT}/${KEY}_1_val_1.fq.gz" ]] && { log "SKIP $KEY (PE trimmed exists)"; continue; }
    else
        [[ -f "${OUT}/${KEY}_trimmed.fq.gz" ]] && { log "SKIP $KEY (SE trimmed exists)"; continue; }
    fi

    wait_slot
    (
        JOB_LOG="$LOG_DIR/individual_jobs/${KEY}.log"
        if [[ "$LAYOUT" == "PE" ]]; then
            R2s=(${tech_r2[$KEY]:-})
            # Merge tech replicates if multiple
            if [[ ${#R1s[@]} -gt 1 ]]; then
                cat "${R1s[@]}" > "/tmp/${KEY}_R1_merged.fq.gz"
                cat "${R2s[@]}" > "/tmp/${KEY}_R2_merged.fq.gz"
                R1_IN="/tmp/${KEY}_R1_merged.fq.gz"
                R2_IN="/tmp/${KEY}_R2_merged.fq.gz"
            else
                R1_IN="${R1s[0]}"; R2_IN="${R2s[0]:-}"
            fi
            trim_galore --paired --cores "${THREADS_TRIMGALORE}" \
                "$R1_IN" "$R2_IN" -o "$OUT" >> "$JOB_LOG" 2>&1
        else
            if [[ ${#R1s[@]} -gt 1 ]]; then
                cat "${R1s[@]}" > "/tmp/${KEY}_SE_merged.fq.gz"
                R1_IN="/tmp/${KEY}_SE_merged.fq.gz"
            else
                R1_IN="${R1s[0]}"
            fi
            trim_galore --cores "${THREADS_TRIMGALORE}" "$R1_IN" -o "$OUT" >> "$JOB_LOG" 2>&1
        fi
    ) &
    pids+=($!); log "STARTED: $KEY [$LAYOUT]"
done
for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== TrimGalore batch complete ==="
