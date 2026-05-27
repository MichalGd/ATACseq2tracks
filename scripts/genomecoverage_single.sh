#!/bin/bash
# fastq2tracks v3.0.2 — Genome coverage: bedGraph + RPM bigWig
# Usage: bash scripts/genomecoverage_single.sh <bam> <hg38|mm39> <outDir>
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

BAM="$1"; GENOME="$2"; OUT="$3"
SAMPLE=$(basename "$BAM" .bam)
mkdir -p "$OUT"
if   [[ "$GENOME" == "hg38" ]]; then CHROM_SIZES="$CHROM_SIZES_HUMAN"; STD_CHR="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
elif [[ "$GENOME" == "mm39" ]]; then CHROM_SIZES="$CHROM_SIZES_MOUSE"; STD_CHR="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
else echo "ERROR: Unknown genome $GENOME" >&2; exit 1; fi
samtools index -@ "${THREADS_BIGWIG}" "$BAM" 2>/dev/null || true
STD_BAM="${OUT}/${SAMPLE}_stChr.bam"
samtools view -@ "${THREADS_BIGWIG}" -b "$BAM" $STD_CHR -o "$STD_BAM"
samtools index -@ "${THREADS_BIGWIG}" "$STD_BAM"
bedtools genomecov -ibam "$BAM" -bga | gzip > "${OUT}/${SAMPLE}.bedGraph.gz"
READS=$(samtools view -@ "${THREADS_BIGWIG}" -c "$STD_BAM")
SCALE=$(awk "BEGIN {printf \"%.10f\", 1000000/$READS}")
bedtools genomecov -ibam "$STD_BAM" -bga -scale "$SCALE" > "${OUT}/${SAMPLE}_norm.bedGraph"
awk '$1~/^chr/' "${OUT}/${SAMPLE}_norm.bedGraph" | \
    LC_COLLATE=C sort -k1,1 -k2,2n > "${OUT}/${SAMPLE}_Snorm.bedGraph"
"${BEDGRAPH_TO_BIGWIG}" "${OUT}/${SAMPLE}_Snorm.bedGraph" "$CHROM_SIZES" "${OUT}/${SAMPLE}_Snorm.bw"
gzip "${OUT}/${SAMPLE}_Snorm.bedGraph"
rm -f "$STD_BAM" "${STD_BAM}.bai" "${OUT}/${SAMPLE}_norm.bedGraph"
echo "Coverage done: ${SAMPLE} reads=${READS} scale=${SCALE}"
