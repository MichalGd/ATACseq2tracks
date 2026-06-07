#!/bin/bash
# ATACseq2tracks v3.0.2 — Blacklist filtering (single sample)
# Usage: bash scripts/blacklist_filter.sh <input.bam> <blacklist.bed> <outDir>
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

INPUT_BAM="$1"; BLACKLIST_BED="$2"; OUT_DIR="$3"
SAMPLE=$(basename "${INPUT_BAM}" .bam)
FILTERED="${OUT_DIR}/${SAMPLE}_blFilt.bam"
mkdir -p "$OUT_DIR"
# Atomic write: use a staging path; only rename to final on full success
STAGED="${FILTERED}.staging.bam"
SORTED_TMP="${FILTERED}.sort_tmp.bam"
rm -f "$STAGED" "$SORTED_TMP"

bedtools intersect -v -abam "$INPUT_BAM" -b "$BLACKLIST_BED" > "$STAGED" || {
    rm -f "$STAGED"; echo "ERROR: bedtools intersect failed for $SAMPLE" >&2; exit 1
}
samtools sort -@ "${THREADS_SAMTOOLS:-1}" "$STAGED" -o "$SORTED_TMP" || {
    rm -f "$STAGED" "$SORTED_TMP"; echo "ERROR: samtools sort failed for $SAMPLE" >&2; exit 1
}
samtools index -@ "${THREADS_SAMTOOLS:-1}" "$SORTED_TMP" || {
    rm -f "$STAGED" "$SORTED_TMP" "${SORTED_TMP}.bai"; echo "ERROR: samtools index failed for $SAMPLE" >&2; exit 1
}
# Atomic rename — only happens if all three steps above succeeded
mv "$SORTED_TMP" "$FILTERED"
mv "${SORTED_TMP}.bai" "${FILTERED}.bai" 2>/dev/null || samtools index -@ "${THREADS_SAMTOOLS:-1}" "$FILTERED"
rm -f "$STAGED"
READS_BEFORE=$(samtools view -@ "${THREADS_SAMTOOLS}" -c "$INPUT_BAM")
READS_AFTER=$(samtools  view -@ "${THREADS_SAMTOOLS}" -c "$FILTERED")
echo "Blacklist filtered: $SAMPLE before=${READS_BEFORE} after=${READS_AFTER} removed=$((READS_BEFORE-READS_AFTER))"
