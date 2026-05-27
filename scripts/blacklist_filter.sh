#!/bin/bash
# fastq2tracks v3.0.2 — Blacklist filtering (single sample)
# Usage: bash scripts/blacklist_filter.sh <input.bam> <blacklist.bed> <outDir>
set -euo pipefail
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to fastq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

INPUT_BAM="$1"; BLACKLIST_BED="$2"; OUT_DIR="$3"
SAMPLE=$(basename "${INPUT_BAM}" .bam)
FILTERED="${OUT_DIR}/${SAMPLE}_blFilt.bam"
mkdir -p "$OUT_DIR"
bedtools intersect -v -abam "$INPUT_BAM" -b "$BLACKLIST_BED" > "${FILTERED}.tmp.bam"
samtools sort -@ "${THREADS_SAMTOOLS}" "${FILTERED}.tmp.bam" -o "$FILTERED"
samtools index -@ "${THREADS_SAMTOOLS}" "$FILTERED"
rm -f "${FILTERED}.tmp.bam"
READS_BEFORE=$(samtools view -@ "${THREADS_SAMTOOLS}" -c "$INPUT_BAM")
READS_AFTER=$(samtools  view -@ "${THREADS_SAMTOOLS}" -c "$FILTERED")
echo "Blacklist filtered: $SAMPLE before=${READS_BEFORE} after=${READS_AFTER} removed=$((READS_BEFORE-READS_AFTER))"
