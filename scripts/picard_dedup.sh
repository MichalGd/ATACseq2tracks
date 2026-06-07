#!/bin/bash
# ATACseq2tracks v3.0.2 — Picard MarkDuplicates (single sample)
# Usage: bash scripts/picard_dedup.sh <inDir> <outDir> <bam_basename>
set -euo pipefail
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to atacseq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

IN="$1"; OUT="$2"; BAM="$3"
mkdir -p "$OUT"
java "${PICARD_XMX}" -jar "${PICARD_JAR}" MarkDuplicates     -INPUT  "${IN}/${BAM}"     -OUTPUT "${OUT}/${BAM}_dedup.bam"     -METRICS_FILE "${OUT}/${BAM}_dedup_rep.txt"     -OPTICAL_DUPLICATE_PIXEL_DISTANCE "${PICARD_OPTICAL_DISTANCE}"     -REMOVE_DUPLICATES true     -ASSUME_SORT_ORDER coordinate     -CREATE_INDEX true     -TMP_DIR "${PICARD_TMP}"
