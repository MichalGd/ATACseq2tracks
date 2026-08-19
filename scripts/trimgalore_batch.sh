#!/usr/bin/env bash
# ATACseq2tracks v4.0.0 - Trim Galore batch with strict input and job validation
set -euo pipefail

[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; IN="${2:?raw FASTQ directory required}"; OUT="${3:?output directory required}"
MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"
SHEET_DIR="$(cd "$(dirname "$SAMPLESHEET")" && pwd)"
LOG_DIR="$(dirname "$OUT")/logs/trimgalore"
mkdir -p "$OUT" "$LOG_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "${LOG_DIR}/main.log"; }
resolve_fastq() {
    local path="$1"
    if [[ -f "$path" ]]; then printf '%s\n' "$path"
    elif [[ -f "${IN}/${path}" ]]; then printf '%s\n' "${IN}/${path}"
    elif [[ -f "${SHEET_DIR}/${path}" ]]; then printf '%s\n' "${SHEET_DIR}/${path}"
    else return 1; fi
}
stem_of() { local name; name="$(basename "$1")"; name="${name%.fastq.gz}"; printf '%s\n' "${name%.fq.gz}"; }

declare -A tech_r1=() tech_r2=() tech_layout=()
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment celltype rep tech_rep is_ctr rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; fq1="${fq1//\"/}"; fq2="${fq2//\"/}"; layout="${layout//\"/}"; rep="${rep//\"/}"
    key="${sid}_bioR${rep}"
    r1_path="$(resolve_fastq "$fq1")" || { log "ERROR: FASTQ not found: $fq1"; exit 1; }
    tech_r1["$key"]+="${r1_path}"$'\n'
    if [[ "$layout" == "PE" ]]; then
        r2_path="$(resolve_fastq "$fq2")" || { log "ERROR: FASTQ not found: $fq2"; exit 1; }
        tech_r2["$key"]+="${r2_path}"$'\n'
    fi
    tech_layout["$key"]="$layout"
done < "$SAMPLESHEET"

declare -a pids=() labels=(); failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: trimming failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}

for key in "${!tech_r1[@]}"; do
    layout="${tech_layout[$key]}"
    [[ "$layout" == "PE" && -s "${OUT}/${key}_1_val_1.fq.gz" && -s "${OUT}/${key}_2_val_2.fq.gz" ]] && { log "SKIP: $key"; continue; }
    [[ "$layout" == "SE" && -s "${OUT}/${key}_trimmed.fq.gz" ]] && { log "SKIP: $key"; continue; }
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    (
        set -euo pipefail
        mapfile -t r1s < <(printf '%s' "${tech_r1[$key]}" | awk 'NF')
        tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks.${key}.XXXXXX")"
        trap 'rm -rf "$tmp_dir"' EXIT
        read -r -a extra_args <<< "${TRIMGALORE_EXTRA_ARGS:-}"
        if [[ "$layout" == "PE" ]]; then
            mapfile -t r2s < <(printf '%s' "${tech_r2[$key]}" | awk 'NF')
            (( ${#r1s[@]} == ${#r2s[@]} )) || { echo "ERROR: unequal R1/R2 count for $key" >&2; exit 1; }
            r1_in="${r1s[0]}"; r2_in="${r2s[0]}"
            if (( ${#r1s[@]} > 1 )); then
                r1_in="${tmp_dir}/${key}_R1.fastq.gz"; r2_in="${tmp_dir}/${key}_R2.fastq.gz"
                gzip -t "${r1s[@]}" "${r2s[@]}"
                command cat "${r1s[@]}" > "$r1_in"; command cat "${r2s[@]}" > "$r2_in"
            fi
            r1_stem="$(stem_of "$r1_in")"; r2_stem="$(stem_of "$r2_in")"
            trim_galore --paired --cores "${THREADS_TRIMGALORE:-2}" "${extra_args[@]}" --output_dir "$OUT" "$r1_in" "$r2_in" >"${LOG_DIR}/${key}.log" 2>&1
            mv "${OUT}/${r1_stem}_val_1.fq.gz" "${OUT}/${key}_1_val_1.fq.gz"
            mv "${OUT}/${r2_stem}_val_2.fq.gz" "${OUT}/${key}_2_val_2.fq.gz"
            gzip -t "${OUT}/${key}_1_val_1.fq.gz" "${OUT}/${key}_2_val_2.fq.gz"
        else
            r1_in="${r1s[0]}"
            if (( ${#r1s[@]} > 1 )); then
                r1_in="${tmp_dir}/${key}_SE.fastq.gz"; gzip -t "${r1s[@]}"; command cat "${r1s[@]}" > "$r1_in"
            fi
            r1_stem="$(stem_of "$r1_in")"
            trim_galore --cores "${THREADS_TRIMGALORE:-2}" "${extra_args[@]}" --output_dir "$OUT" "$r1_in" >"${LOG_DIR}/${key}.log" 2>&1
            mv "${OUT}/${r1_stem}_trimmed.fq.gz" "${OUT}/${key}_trimmed.fq.gz"; gzip -t "${OUT}/${key}_trimmed.fq.gz"
        fi
    ) &
    pids+=("$!"); labels+=("$key")
done
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} trimming job(s) failed"; exit 1; }
log "Trim Galore batch complete"
