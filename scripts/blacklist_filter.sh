#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Blacklist / exclusion region filtering
# Usage: bash scripts/blacklist_filter.sh <input.bam> <blacklist.bed> <output_dir>
#
# Removes reads overlapping blacklisted regions using bedtools intersect.
# Applied AFTER deduplication, BEFORE coverage/normalisation/MACS2.
# Output: <sample>_blFilt.bam + .bai index
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

INPUT_BAM="$1"
BLACKLIST_BED="$2"
OUT_DIR="$3"

SAMPLE=$(basename "${INPUT_BAM}" .bam)
FILTERED="${OUT_DIR}/${SAMPLE}_blFilt.bam"

mkdir -p "$OUT_DIR"

bedtools intersect -v -abam "$INPUT_BAM" -b "$BLACKLIST_BED" > "${FILTERED}.tmp.bam"
samtools sort -@ "${THREADS_SAMTOOLS}" "${FILTERED}.tmp.bam" -o "$FILTERED"
samtools index -@ "${THREADS_SAMTOOLS}" "$FILTERED"
rm -f "${FILTERED}.tmp.bam"

READS_BEFORE=$(samtools view -@ "${THREADS_SAMTOOLS}" -c "$INPUT_BAM")
READS_AFTER=$(samtools view  -@ "${THREADS_SAMTOOLS}" -c "$FILTERED")
echo "Blacklist filtered: $SAMPLE  before=${READS_BEFORE}  after=${READS_AFTER}  removed=$((READS_BEFORE-READS_AFTER))"
