#!/bin/bash
# ATACseq2tracks v3.0.4 — TrimGalore batch (SE + PE), trim_galore v2.2.0 Oxidized Edition
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

SAMPLESHEET_ARG="$1"; IN="$2"; OUT="$3"
MAX_JOBS="${THREADS_PARALLEL_JOBS}"
mkdir -p "$OUT"
LOG_DIR="${OUT}/../logs/trimgalore"; mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/main.log"; }

declare -a pids=()
wait_slot() {
    while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
        local new=()
        for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && new+=("$p"); done
        pids=("${new[@]+"${new[@]}"}"); sleep 2
    done
}

# stem_of: strip directory and both .fastq.gz / .fq.gz extensions
stem_of() { local b; b=$(basename "$1"); b="${b%.fastq.gz}"; echo "${b%.fq.gz}"; }

declare -A tech_r1 tech_r2 tech_layout

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment celltype rep tech_rep is_ctr rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; fq1="${fq1//\"/}"; fq2="${fq2//\"/}"
    layout="${layout//\"/}"; rep="${rep//\"/}"
    KEY="${sid}_bioR${rep}"
    tech_r1["$KEY"]+="${fq1} "
    [[ -n "$fq2" && "$layout" == "PE" ]] && tech_r2["$KEY"]+="${fq2} "
    tech_layout["$KEY"]="${layout}"
done < "$SAMPLESHEET_ARG"

log "=== TrimGalore batch: ${#tech_r1[@]} sample groups ==="

for KEY in "${!tech_r1[@]}"; do
    R1s=(${tech_r1[$KEY]})
    LAYOUT="${tech_layout[$KEY]}"

    if [[ "$LAYOUT" == "PE" ]]; then
        [[ -f "${OUT}/${KEY}_1_val_1.fq.gz" ]] && { log "SKIP $KEY (PE trimmed exists)"; continue; }
    else
        [[ -f "${OUT}/${KEY}_trimmed.fq.gz" ]] && { log "SKIP $KEY (SE trimmed exists)"; continue; }
    fi

    wait_slot
    (
        set -euo pipefail
        JOB_LOG="${LOG_DIR}/${KEY}.log"

        if [[ "$LAYOUT" == "PE" ]]; then
            R2s=(${tech_r2[$KEY]:-})

            # Merge tech reps or use single files
            if [[ ${#R1s[@]} -gt 1 ]]; then
                R1_IN="/tmp/${KEY}_R1.fq.gz"
                R2_IN="/tmp/${KEY}_R2.fq.gz"
                cat "${R1s[@]}" > "$R1_IN"
                cat "${R2s[@]}" > "$R2_IN"
                MERGED=1
            else
                R1_IN="${R1s[0]}"
                R2_IN="${R2s[0]}"
                MERGED=0
            fi

            R1_STEM=$(stem_of "$R1_IN")
            R2_STEM=$(stem_of "$R2_IN")

            trim_galore --paired --cores "${THREADS_TRIMGALORE}" \
                --output_dir "$OUT" "$R1_IN" "$R2_IN" >> "$JOB_LOG" 2>&1

            # Rename to canonical KEY names
            # v2.2.0 PE output: <stem>_val_1.fq.gz / <stem>_val_2.fq.gz
            if [[ -f "${OUT}/${R1_STEM}_val_1.fq.gz" ]]; then
                mv -f "${OUT}/${R1_STEM}_val_1.fq.gz" "${OUT}/${KEY}_1_val_1.fq.gz"
                mv -f "${OUT}/${R2_STEM}_val_2.fq.gz" "${OUT}/${KEY}_2_val_2.fq.gz"
            elif [[ -f "${OUT}/${R1_STEM}_trimmed.fq.gz" ]]; then
                # fallback: v2.2.0 may use _trimmed for both
                mv -f "${OUT}/${R1_STEM}_trimmed.fq.gz" "${OUT}/${KEY}_1_val_1.fq.gz"
                mv -f "${OUT}/${R2_STEM}_trimmed.fq.gz" "${OUT}/${KEY}_2_val_2.fq.gz"
            else
                echo "ERROR [$KEY]: cannot find PE trimmed output (looked for ${R1_STEM}_val_1 or _trimmed)" \
                    | tee -a "$JOB_LOG" >&2
                ls "$OUT"/*.fq.gz 2>/dev/null | tee -a "$JOB_LOG" >&2
                exit 1
            fi

            [[ $MERGED -eq 1 ]] && rm -f "$R1_IN" "$R2_IN"

        else
            if [[ ${#R1s[@]} -gt 1 ]]; then
                R1_IN="/tmp/${KEY}_SE.fq.gz"
                cat "${R1s[@]}" > "$R1_IN"
                MERGED=1
            else
                R1_IN="${R1s[0]}"
                MERGED=0
            fi

            R1_STEM=$(stem_of "$R1_IN")

            trim_galore --cores "${THREADS_TRIMGALORE}" \
                --output_dir "$OUT" "$R1_IN" >> "$JOB_LOG" 2>&1

            if [[ -f "${OUT}/${R1_STEM}_trimmed.fq.gz" ]]; then
                mv -f "${OUT}/${R1_STEM}_trimmed.fq.gz" "${OUT}/${KEY}_trimmed.fq.gz"
            else
                echo "ERROR [$KEY]: cannot find SE trimmed output (looked for ${R1_STEM}_trimmed)" \
                    | tee -a "$JOB_LOG" >&2
                ls "$OUT"/*.fq.gz 2>/dev/null | tee -a "$JOB_LOG" >&2
                exit 1
            fi

            [[ $MERGED -eq 1 ]] && rm -f "$R1_IN"
        fi

        log "DONE: $KEY [$LAYOUT]"
    ) &
    pids+=($!)
    log "STARTED: $KEY [$LAYOUT]"
done

for p in "${pids[@]+"${pids[@]}"}"; do wait "$p" || true; done
log "=== TrimGalore batch complete ==="
