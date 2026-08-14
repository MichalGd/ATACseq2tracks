#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - merge filtered biological samples for group-level CPM tracks
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
SAMPLESHEET="${1:?samplesheet required}"; IN_DIR="${2:?filtered BAM directory required}"; OUT_DIR="${3:?output directory required}"; GENOME="${4:?genome required}"
LOG="${OUT_DIR}/merge_replicates.log"; mkdir -p "$OUT_DIR"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$LOG"; }
declare -A group_bams=() group_layout=() seen=()
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctrl rest; do
    [[ "$sid" == "sample_id" ]] && continue
    for name in sid layout genome factor condition treatment cell_type rep is_ctrl; do value="${!name}"; printf -v "$name" '%s' "${value//\"/}"; done
    [[ "$genome" != "$GENOME" || "${is_ctrl,,}" =~ ^(true|1|yes)$ ]] && continue
    key="${sid}_bioR${rep}"; [[ -n "${seen[$key]+x}" ]] && continue; seen["$key"]=1
    bam="${IN_DIR}/${key}_dedup_blFilt.bam"; [[ -s "$bam" ]] || { log "ERROR: filtered BAM missing: $bam"; exit 1; }
    group="${factor}__${condition}__${treatment}__${cell_type}"
    if [[ -n "${group_layout[$group]:-}" && "${group_layout[$group]}" != "${layout^^}" ]]; then
        log "ERROR: cannot merge PE and SE samples into one signal track: $group"
        exit 1
    fi
    group_layout["$group"]="${layout^^}"
    group_bams["$group"]+="$bam "
done < "$SAMPLESHEET"
for group in "${!group_bams[@]}"; do
    read -r -a bams <<< "${group_bams[$group]}"
    (( ${#bams[@]} >= 2 )) || { log "SKIP group with one biological sample: $group"; continue; }
    merged_name="${group//__/_}_merged"; merged_bam="${OUT_DIR}/${merged_name}.bam"
    samtools merge -@ "${THREADS_SAMTOOLS:-2}" -f "$merged_bam" "${bams[@]}"
    samtools index -@ "${THREADS_SAMTOOLS:-2}" "$merged_bam"
    bash "$(dirname "$0")/genomecoverage_single.sh" "$merged_bam" "$GENOME" "$OUT_DIR" "${group_layout[$group]}" >> "$LOG" 2>&1
done
log "Group merge and CPM track generation complete"
