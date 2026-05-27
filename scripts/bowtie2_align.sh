#!/bin/bash
# fastq2tracks v3.0.2 — Bowtie2 alignment (PE and SE)
# Usage: bash scripts/bowtie2_align.sh <index> <R1> <outDir> <PE|SE> [R2]
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

INDEX="$1"; R1="$2"; OUT="$3"; LAYOUT="${4:-PE}"; R2="${5:-}"
SAMPLE=$(basename "$R1" | sed 's/_R1_001_val_1\.fq\.gz//; s/_1_val_1\.fq\.gz//; s/_trimmed\.fq\.gz//')
SAM="${OUT}/${SAMPLE}.sam"; SORTED="${OUT}/${SAMPLE}.sorted_stChr.bam"
LOG="${OUT}/${SAMPLE}.bowtie2.txt"; METRICS="${OUT}/${SAMPLE}.bowtie2.metrics.txt"
mkdir -p "$OUT"
if [[ "$LAYOUT" == "PE" ]]; then
    [[ -z "$R2" ]] && { echo "ERROR: PE requires R2" >&2; exit 1; }
    bowtie2 --minins 20 --maxins 1200 --dovetail -q --phred33         -p "${THREADS_ALIGN}" --no-unal -x "$INDEX" -1 "$R1" -2 "$R2"         -S "$SAM" --met-file "$METRICS" 2>"$LOG"
else
    bowtie2 -q --phred33 -p "${THREADS_ALIGN}" --no-unal         -x "$INDEX" -U "$R1" -S "$SAM" --met-file "$METRICS" 2>"$LOG"
fi
samtools view -@ "${THREADS_SAMTOOLS}" -bS "$SAM" |     samtools sort -@ "${THREADS_SAMTOOLS}" -o "$SORTED"
samtools index -@ "${THREADS_SAMTOOLS}" "$SORTED"
rm -f "$SAM"
echo "Alignment done: $SORTED"
