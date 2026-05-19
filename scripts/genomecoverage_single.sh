#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Genome coverage: raw bedGraph + RPM-normalised bedGraph + bigWig
# Usage: bash scripts/genomecoverage_single.sh <bam> <genome> <outDir>
#   genome: hg38 or mm39  (uses local chrom.sizes from config)
# Blacklist filtering must be applied BEFORE this script (use _blFilt.bam).
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

BAM="$1"; GENOME="$2"; OUT="$3"
SAMPLE=$(basename "$BAM" .bam)
mkdir -p "$OUT"

if [[ "$GENOME" == "hg38" ]]; then
    CHROM_SIZES="$CHROM_SIZES_HUMAN"
    STD_CHR="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
elif [[ "$GENOME" == "mm39" ]]; then
    CHROM_SIZES="$CHROM_SIZES_MOUSE"
    STD_CHR="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
else
    echo "ERROR: Unknown genome $GENOME" >&2; exit 1
fi

# re-index if needed
samtools index -@ "${THREADS_BIGWIG}" "$BAM" 2>/dev/null || true

# Standard-chromosome BAM (for normalisation)
STD_BAM="${OUT}/${SAMPLE}_stChr.bam"
samtools view -@ "${THREADS_BIGWIG}" -b "$BAM" $STD_CHR -o "$STD_BAM"
samtools index -@ "${THREADS_BIGWIG}" "$STD_BAM"

# Raw bedGraph (all chromosomes)
bedtools genomecov -ibam "$BAM" -bga | gzip > "${OUT}/${SAMPLE}.bedGraph.gz"

# RPM-normalised bigWig
READS=$(samtools view -@ "${THREADS_BIGWIG}" -c "$STD_BAM")
SCALE=$(awk "BEGIN {printf \"%.10f\", 1000000/$READS}")
bedtools genomecov -ibam "$STD_BAM" -bga -scale "$SCALE" > "${OUT}/${SAMPLE}_norm.bedGraph"
awk '$1~/^chr/' "${OUT}/${SAMPLE}_norm.bedGraph" | \
    LC_COLLATE=C sort -k1,1 -k2,2n > "${OUT}/${SAMPLE}_Snorm.bedGraph"
"${BEDGRAPH_TO_BIGWIG}" "${OUT}/${SAMPLE}_Snorm.bedGraph" "$CHROM_SIZES" "${OUT}/${SAMPLE}_Snorm.bw"
gzip "${OUT}/${SAMPLE}_Snorm.bedGraph"

rm -f "$STD_BAM" "${STD_BAM}.bai" "${OUT}/${SAMPLE}_norm.bedGraph"
echo "Coverage done: ${SAMPLE}  reads=${READS}  scale=${SCALE}"
