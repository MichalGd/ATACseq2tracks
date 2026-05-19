#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Global Configuration
# Run from the fastq2tracks/ root directory: bash fastq2tracks.sh config/samplesheet.csv
# =============================================================================

# ---------------------------------------------------------------------------
# Thread / worker settings  (70 physical / 140 logical cores, 500 GB RAM)
# ---------------------------------------------------------------------------
THREADS_ALIGN=16          # bowtie2 -p per sample (8 parallel = 128 total)
THREADS_SAMTOOLS=16       # samtools -@ per sample
THREADS_PARALLEL_JOBS=8   # max parallel samples per batch stage
THREADS_FASTQC=10         # fastqc --threads per job (2x jobs = 20 per slot)
THREADS_TRIMGALORE=8      # trim_galore --cores per job
THREADS_BIGWIG=16         # samtools / bedtools threads per job
THREADS_CHIPQC=20         # BiocParallel workers for ChIPQC

# ---------------------------------------------------------------------------
# Reference genomes — bowtie2 indices
# ---------------------------------------------------------------------------
INDEX_HG38="/home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38"
INDEX_MM39="/home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39"

# ---------------------------------------------------------------------------
# Annotation files
# ---------------------------------------------------------------------------
GTF_HUMAN="/home/micgdu/GenomicData/genomesDec2022/gencode.v42.primary_assembly.annotation.gtf"
CHROM_SIZES_HUMAN="/home/micgdu/GenomicData/genomesDec2022/hs38n.chrom.sizes"

GTF_MOUSE="/home/micgdu/GenomicData/genomesDec2022/gencode.vM31.primary_assembly.annotation.gtf"
CHROM_SIZES_MOUSE="/home/micgdu/GenomicData/genomesDec2022/mm39n.chrom.sizes"

# ---------------------------------------------------------------------------
# Blacklist / exclusion region BED files
# ---------------------------------------------------------------------------
BLACKLIST_HG38="/home/micgdu/data/shared/AnnotationHub_cache/blacklist_hg38_ENCFF356LFX.bed"
BLACKLIST_MM39="/home/micgdu/data/shared/AnnotationHub_cache/blacklist_mm39_Boyle.bed"

# ---------------------------------------------------------------------------
# ChIPQC annotation objects (pre-built RDS)
# ---------------------------------------------------------------------------
CHIPQC_ANNOTATION_HG38="/home/micgdu/data/shared/AnnotationHub_cache/anno_hg38_chipqc.rds"
CHIPQC_ANNOTATION_MM39="/home/micgdu/data/shared/AnnotationHub_cache/anno_mm39_chipqc.rds"

CHIPQC_BLACKLIST_HG38_RDS="/home/micgdu/data/shared/AnnotationHub_cache/blacklist_hg38.rds"
CHIPQC_BLACKLIST_MM39_RDS="/home/micgdu/data/shared/AnnotationHub_cache/blacklist_mm39.rds"

# ---------------------------------------------------------------------------
# Software paths
# ---------------------------------------------------------------------------
PICARD_JAR="/home/micgdu/software/picard.jar"
BEDGRAPH_TO_BIGWIG="/home/micgdu/kentutils/bedGraphToBigWig"
MACS2_BIN="macs2"
R_BIN="Rscript"
CONDA_ENV_ACTIVATE="/home/micgdu/myenv/bin/activate"

# ---------------------------------------------------------------------------
# Picard settings
# ---------------------------------------------------------------------------
PICARD_XMX="-Xmx32G"
PICARD_OPTICAL_DISTANCE=2500
PICARD_TMP="/tmp"

# ---------------------------------------------------------------------------
# MACS2 settings
# ---------------------------------------------------------------------------
MACS2_GENOME_HG38="hs"
MACS2_GENOME_MM39="mm"
MACS2_QVALUE=0.05
MACS2_BROAD_CUTOFF=0.1

# ---------------------------------------------------------------------------
# Output control
# ---------------------------------------------------------------------------
KEEP_INTERMEDIATE_BAMS=false   # set true to keep pre-dedup bams
KEEP_TRIMMED_FASTQ=false        # set true to keep trimmed fastqs
