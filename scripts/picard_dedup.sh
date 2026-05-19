#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Picard MarkDuplicates (single sample)
# Usage: bash scripts/picard_dedup.sh <inDir> <outDir> <bamFileName>
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

IN="$1"; OUT="$2"; BAM="$3"
mkdir -p "$OUT"

java "${PICARD_XMX}" -jar "${PICARD_JAR}" MarkDuplicates \
    -INPUT  "${IN}/${BAM}" \
    -OUTPUT "${OUT}/${BAM}_dedup.bam" \
    -METRICS_FILE "${OUT}/${BAM}_dedup_rep.txt" \
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE "${PICARD_OPTICAL_DISTANCE}" \
    -REMOVE_DUPLICATES true \
    -ASSUME_SORT_ORDER coordinate \
    -CREATE_INDEX true \
    -TMP_DIR "${PICARD_TMP}"
