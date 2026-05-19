#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Bowtie2 alignment (PE and SE)
# Usage: bash scripts/bowtie2_align.sh <index> <fastq_R1> <outFolder> <layout> [fastq_R2]
#   index      : bowtie2 index prefix
#   fastq_R1   : trimmed R1 (or only) FASTQ
#   outFolder  : output directory for BAM files
#   layout     : PE or SE
#   fastq_R2   : (required for PE) trimmed R2 FASTQ
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

INDEX="$1"
R1="$2"
OUT="$3"
LAYOUT="${4:-PE}"
R2="${5:-}"

SAMPLE=$(basename "$R1" | sed 's/_R1_001_val_1\.fq\.gz//; s/_1_val_1\.fq\.gz//; s/_trimmed\.fq\.gz//')
SAM="${OUT}/${SAMPLE}.sam"
BAM="${OUT}/${SAMPLE}.bam"
SORTED="${OUT}/${SAMPLE}.sorted_stChr.bam"
LOG="${OUT}/${SAMPLE}.bowtie2.txt"
METRICS="${OUT}/${SAMPLE}.bowtie2.metrics.txt"

mkdir -p "$OUT"

if [[ "$LAYOUT" == "PE" ]]; then
    [[ -z "$R2" ]] && { echo "ERROR: PE layout requires R2 FASTQ" >&2; exit 1; }
    bowtie2 --minins 20 --maxins 1200 --dovetail -q --phred33 \
        -p "${THREADS_ALIGN}" --no-unal \
        -x "$INDEX" -1 "$R1" -2 "$R2" \
        -S "$SAM" --met-file "$METRICS" 2>"$LOG"
else
    bowtie2 -q --phred33 \
        -p "${THREADS_ALIGN}" --no-unal \
        -x "$INDEX" -U "$R1" \
        -S "$SAM" --met-file "$METRICS" 2>"$LOG"
fi

samtools view -@ "${THREADS_SAMTOOLS}" -bS "$SAM" | \
    samtools sort -@ "${THREADS_SAMTOOLS}" -o "$SORTED"
samtools index -@ "${THREADS_SAMTOOLS}" "$SORTED"

mv "$LOG" "$METRICS" "$OUT/"
rm -f "$SAM" "$BAM"

echo "Alignment done: $SORTED"
