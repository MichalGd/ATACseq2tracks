#!/bin/bash
# =============================================================================
# fastq2tracks v3.0
# ChIP-seq / NGS multi-genome track-generation and QC workflow
#
# Processes paired-end (PE) and single-end (SE) Illumina data.
# Produces: deduplicated BAMs, blacklist-filtered BAMs, RPM-normalised bigWig
# tracks, MACS2 narrow/broad peak calls, merged-replicate tracks, ChIPQC
# report, and DiffBind-compatible output.
#
# Supported genomes : hg38 (GRCh38), mm39 (GRCm39)
# Sample metadata   : defined in config/samplesheet.csv
# Configuration     : config/config.sh
#
# Usage:
#   bash fastq2tracks.sh config/samplesheet.csv <outputDir> [max_parallel_jobs]
#
# Must be run from the fastq2tracks/ root directory.
# =============================================================================
set -euo pipefail

SAMPLESHEET="${1:?ERROR: samplesheet required as first argument}"
OUTDIR="${2:?ERROR: output directory required as second argument}"
MAX_JOBS="${3:-8}"

CONFIG="config/config.sh"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
source "$CONFIG"
THREADS_PARALLEL_JOBS="$MAX_JOBS"

SCRIPT_DIR="$(cd "$(dirname "$0")/scripts" && pwd)"

# --- Step 0: Pre-flight smoke test ---
echo "=== [0] Pre-flight checks ==="
bash "${SCRIPT_DIR}/smoke_test.sh" "$SAMPLESHEET" "$CONFIG"

# --- Folder setup ---
mkdir -p \
    "${OUTDIR}/fastQC/fastQC_unTrimmed" \
    "${OUTDIR}/fastQC/fastQC_trimmed" \
    "${OUTDIR}/multiQC/multiQC_unTrimmed" \
    "${OUTDIR}/multiQC/multiQC_trimmed" \
    "${OUTDIR}/multiQC/multiQC_alignments" \
    "${OUTDIR}/multiQC/multiQC_deduplication" \
    "${OUTDIR}/trimmedFastq" \
    "${OUTDIR}/bams" \
    "${OUTDIR}/dedupBams" \
    "${OUTDIR}/filteredBams" \
    "${OUTDIR}/bedGraph" \
    "${OUTDIR}/NormBedGraph" \
    "${OUTDIR}/bigwig" \
    "${OUTDIR}/bigwig_merged" \
    "${OUTDIR}/peaks" \
    "${OUTDIR}/chipqc" \
    "${OUTDIR}/diffbind" \
    "${OUTDIR}/reports"

# --- Determine input FASTQ directory from samplesheet ---
# FASTQs are specified with absolute paths in the samplesheet; no fixed input folder needed.
RAW_FASTQ_DIR=$(dirname "$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f2 | tr -d '"')")

# --- Step 1: FastQC on raw ---
echo "=== [1] FastQC (raw) ==="
bash "${SCRIPT_DIR}/fastqc_batch.sh" "$RAW_FASTQ_DIR" "${OUTDIR}/fastQC/fastQC_unTrimmed" "$MAX_JOBS"
source "$CONDA_ENV_ACTIVATE"
multiqc "${OUTDIR}/fastQC/fastQC_unTrimmed" -n multiQC_unTrimmed \
    -o "${OUTDIR}/multiQC/multiQC_unTrimmed" --data-format tsv --export

# --- Step 2: Trimming ---
echo "=== [2] TrimGalore ==="
bash "${SCRIPT_DIR}/trimgalore_batch.sh" "$SAMPLESHEET" "$RAW_FASTQ_DIR" "${OUTDIR}/trimmedFastq"

# --- Step 3: FastQC on trimmed ---
echo "=== [3] FastQC (trimmed) ==="
bash "${SCRIPT_DIR}/fastqc_batch.sh" "${OUTDIR}/trimmedFastq" "${OUTDIR}/fastQC/fastQC_trimmed" "$MAX_JOBS"
multiqc "${OUTDIR}/fastQC/fastQC_trimmed" -n multiQC_trimmed \
    -o "${OUTDIR}/multiQC/multiQC_trimmed" --data-format tsv --export

# --- Step 4: Alignment ---
echo "=== [4] Bowtie2 alignment ==="
bash "${SCRIPT_DIR}/bowtie2_batch.sh" "$SAMPLESHEET" "${OUTDIR}/trimmedFastq" "${OUTDIR}/bams"
multiqc "${OUTDIR}/bams" -n multiQC_alignments \
    -o "${OUTDIR}/multiQC/multiQC_alignments" --data-format tsv --export

# --- Step 5: Deduplication ---
echo "=== [5] Picard deduplication ==="
bash "${SCRIPT_DIR}/picard_dedup_batch.sh" "${OUTDIR}/bams" "${OUTDIR}/dedupBams" "$MAX_JOBS"
multiqc "${OUTDIR}/dedupBams" -n multiQC_deduplication \
    -o "${OUTDIR}/multiQC/multiQC_deduplication" --data-format tsv --export

# --- Step 6: Blacklist filtering ---
echo "=== [6] Blacklist filtering ==="
bash "${SCRIPT_DIR}/blacklist_filter_batch.sh" "$SAMPLESHEET" "${OUTDIR}/dedupBams" "${OUTDIR}/filteredBams"

# --- Step 7: Genome coverage (individual samples) ---
echo "=== [7] Genome coverage (individual) ==="
bash "${SCRIPT_DIR}/genomecoverage_batch.sh" "$SAMPLESHEET" "${OUTDIR}/filteredBams" "${OUTDIR}/bigwig"
# Move bedGraphs
mv "${OUTDIR}/bigwig"/*.bedGraph.gz "${OUTDIR}/bedGraph/" 2>/dev/null || true
mv "${OUTDIR}/bigwig"/*Snorm*.bedGraph.gz "${OUTDIR}/NormBedGraph/" 2>/dev/null || true

# --- Step 8: Merge replicates + merged coverage ---
echo "=== [8] Replicate merging + merged tracks ==="
# Determine genomes present in samplesheet
for GENOME in hg38 mm39; do
    if grep -q "$GENOME" "$SAMPLESHEET"; then
        bash "${SCRIPT_DIR}/merge_replicates.sh" "$SAMPLESHEET" "${OUTDIR}/filteredBams" "${OUTDIR}/bigwig_merged" "$GENOME"
    fi
done

# --- Step 9: MACS2 peak calling ---
echo "=== [9] MACS2 peak calling ==="
bash "${SCRIPT_DIR}/macs2_batch.sh" "$SAMPLESHEET" "${OUTDIR}/filteredBams" "${OUTDIR}/peaks"

# --- Step 10: ChIPQC ---
echo "=== [10] ChIPQC ==="
for GENOME in hg38 mm39; do
    if grep -q "$GENOME" "$SAMPLESHEET"; then
        "${R_BIN}" "${SCRIPT_DIR}/run_chipqc.R" \
            "$SAMPLESHEET" "${OUTDIR}/filteredBams" "${OUTDIR}/peaks" \
            "${OUTDIR}/chipqc" "$GENOME" "${THREADS_CHIPQC}"
    fi
done

# --- Step 11: Prepare DiffBind samplesheet ---
echo "=== [11] DiffBind samplesheet preparation ==="
for GENOME in hg38 mm39; do
    if grep -q "$GENOME" "$SAMPLESHEET"; then
        "${R_BIN}" "${SCRIPT_DIR}/prepare_diffbind.R" \
            "$SAMPLESHEET" "${OUTDIR}/filteredBams" "${OUTDIR}/peaks" \
            "${OUTDIR}/diffbind" "$GENOME"
    fi
done

# --- Step 12: UCSC track lines ---
echo "=== [12] UCSC tracks ==="
bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" "${OUTDIR}/bigwig" "http://your-server.com/data"

# --- Step 13: MultiQC unified report ---
echo "=== [13] MultiQC unified report ==="
bash "${SCRIPT_DIR}/generate_pipeline_report.sh" "$OUTDIR" \
    "${OUTDIR}/reports/pipeline_report_$(date +%Y%m%d)" html

# --- Cleanup ---
if [[ "$KEEP_INTERMEDIATE_BAMS" == "false" ]]; then
    rm -f "${OUTDIR}/bams"/*.bam "${OUTDIR}/bams"/*.bai || true
fi
if [[ "$KEEP_TRIMMED_FASTQ" == "false" ]]; then
    rm -f "${OUTDIR}/trimmedFastq"/*.gz || true
fi

echo "=== fastq2tracks v3.0 complete ==="
echo "Outputs in: $OUTDIR"
