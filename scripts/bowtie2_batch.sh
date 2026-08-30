#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - fail-fast Bowtie2 alignment batch
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; TRIM_DIR="${2:?trimmed FASTQ directory required}"; OUT_DIR="${3:?output directory required}"
LOG_DIR="${OUT_DIR}/bowtie2_logs"; MAX_JOBS="${THREADS_PARALLEL_JOBS:-2}"
mkdir -p "$OUT_DIR" "$LOG_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "${LOG_DIR}/main.log"; }
declare -a pids=() labels=(); declare -A seen_keys=(); failures=0
wait_one() {
    local pid="${pids[0]}" label="${labels[0]}"
    if ! wait "$pid"; then log "ERROR: alignment failed: $label"; failures=$((failures + 1)); fi
    pids=("${pids[@]:1}"); labels=("${labels[@]:1}")
}
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment celltype rep tech_rep is_ctr rest; do
    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; layout="${layout//\"/}"; genome="${genome//\"/}"; rep="${rep//\"/}"
    key="${sid}_bioR${rep}"; [[ -n "${seen_keys[$key]+x}" ]] && continue; seen_keys["$key"]=1
    case "$genome" in hg38) index="$INDEX_HG38" ;; mm39) index="$INDEX_MM39" ;; *) log "ERROR: unsupported genome: $genome"; exit 1 ;; esac
    bam_out="${OUT_DIR}/${key}.bam"
    if [[ -s "$bam_out" ]] && samtools quickcheck "$bam_out" 2>/dev/null; then log "SKIP: $key"; continue; fi
    while (( ${#pids[@]} >= MAX_JOBS )); do wait_one; done
    (
        set -euo pipefail
        staging="${OUT_DIR}/.${key}.staging.$$.bam"; trap 'rm -f "$staging" "${staging}.bai"' EXIT
        read -r -a extra_args <<< "${BOWTIE2_EXTRA_ARGS:---very-sensitive}"
        common=( -x "$index" -p "${THREADS_ALIGN:-2}" --rg-id "$key" --rg "SM:${key}" --rg 'PL:ILLUMINA' )
        if [[ "$layout" == "PE" ]]; then
            r1="${TRIM_DIR}/${key}_1_val_1.fq.gz"; r2="${TRIM_DIR}/${key}_2_val_2.fq.gz"
            [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: missing trimmed PE files for $key" >&2; exit 1; }
            bowtie2 "${common[@]}" "${extra_args[@]}" --no-mixed --no-discordant --dovetail -1 "$r1" -2 "$r2" 2>"${LOG_DIR}/${key}.log" | samtools sort -@ "${THREADS_SAMTOOLS:-2}" -o "$staging"
        else
            r1="${TRIM_DIR}/${key}_trimmed.fq.gz"; [[ -s "$r1" ]] || { echo "ERROR: missing trimmed SE file for $key" >&2; exit 1; }
            bowtie2 "${common[@]}" "${extra_args[@]}" -U "$r1" 2>"${LOG_DIR}/${key}.log" | samtools sort -@ "${THREADS_SAMTOOLS:-2}" -o "$staging"
        fi
        samtools quickcheck "$staging"; samtools index -@ "${THREADS_SAMTOOLS:-2}" "$staging"
        mv "$staging" "$bam_out"; mv "${staging}.bai" "${bam_out}.bai"
        samtools flagstat -@ "${THREADS_SAMTOOLS:-2}" "$bam_out" > "${LOG_DIR}/${key}.flagstat"
        samtools stats -@ "${THREADS_SAMTOOLS:-2}" "$bam_out" > "${LOG_DIR}/${key}.stats"
        trap - EXIT
    ) &
    pids+=("$!"); labels+=("$key")
done < "$SAMPLESHEET"
while (( ${#pids[@]} > 0 )); do wait_one; done
(( failures == 0 )) || { log "FATAL: ${failures} alignment job(s) failed"; exit 1; }
log "Bowtie2 alignment batch complete"
