#!/bin/bash
# =============================================================================
# fastq2tracks v3.0.4
# ChIP-seq / NGS multi-genome track-generation and QC workflow
#
# Usage:
#   bash /path/to/fastq2tracks/fastq2tracks.sh --config /path/to/config.conf
#
# Checkpoint system: completed steps are skipped on re-run.
# To force a step to re-run:
#   rm /your/output/.checkpoints/stepN.done
# To force ALL steps to re-run:
#   rm -rf /your/output/.checkpoints/
# =============================================================================
set -euo pipefail

INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${INSTALL_DIR}/scripts"
F2T_CONFIG=""

usage() {
    echo "Usage: bash fastq2tracks.sh --config /absolute/path/to/config.conf"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)  F2T_CONFIG="$(realpath "$2")"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done

[[ -z "$F2T_CONFIG" ]] && { echo "ERROR: --config is required" >&2; usage; }
[[ -f "$F2T_CONFIG" ]] || { echo "ERROR: config file not found: $F2T_CONFIG" >&2; exit 1; }

export F2T_CONFIG
source "$F2T_CONFIG"

[[ -z "${SAMPLESHEET:-}" ]] && { echo "ERROR: SAMPLESHEET not set in $F2T_CONFIG" >&2; exit 1; }
[[ -z "${OUTPUT_DIR:-}"  ]] && { echo "ERROR: OUTPUT_DIR not set in $F2T_CONFIG"  >&2; exit 1; }
[[ -f "$SAMPLESHEET"     ]] || { echo "ERROR: SAMPLESHEET not found: $SAMPLESHEET" >&2; exit 1; }

echo "============================================================"
echo " fastq2tracks v3.0.4"
echo " Install dir : $INSTALL_DIR"
echo " Config      : $F2T_CONFIG"
echo " Samplesheet : $SAMPLESHEET"
echo " Output      : $OUTPUT_DIR"
echo " Max jobs    : ${THREADS_PARALLEL_JOBS}"
echo "============================================================"

# ── Checkpoint helpers ────────────────────────────────────────────────────────
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "$CHECKPOINT_DIR"

is_done()   { [[ -f "${CHECKPOINT_DIR}/step${1}.done" ]]; }
mark_done() {
    touch "${CHECKPOINT_DIR}/step${1}.done"
    echo "[CHECKPOINT] Step ${1} complete -- to re-run: rm ${CHECKPOINT_DIR}/step${1}.done"
}
skip_msg()  { echo "=== [${1}] SKIPPED (already complete -- checkpoint exists) ==="; }

# ── Create output folders ─────────────────────────────────────────────────────
mkdir -p \
    "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" \
    "${OUTPUT_DIR}/fastQC/fastQC_trimmed" \
    "${OUTPUT_DIR}/multiQC/multiQC_unTrimmed" \
    "${OUTPUT_DIR}/multiQC/multiQC_trimmed" \
    "${OUTPUT_DIR}/multiQC/multiQC_alignments" \
    "${OUTPUT_DIR}/multiQC/multiQC_deduplication" \
    "${OUTPUT_DIR}/trimmedFastq" \
    "${OUTPUT_DIR}/bams" \
    "${OUTPUT_DIR}/dedupBams" \
    "${OUTPUT_DIR}/filteredBams" \
    "${OUTPUT_DIR}/bedGraph" \
    "${OUTPUT_DIR}/NormBedGraph" \
    "${OUTPUT_DIR}/bigwig" \
    "${OUTPUT_DIR}/bigwig_merged" \
    "${OUTPUT_DIR}/peaks/per_replicate" \
    "${OUTPUT_DIR}/peaks/pooled" \
    "${OUTPUT_DIR}/chipqc" \
    "${OUTPUT_DIR}/diffbind" \
    "${OUTPUT_DIR}/reports"

[[ -n "${CONDA_ENV_ACTIVATE:-}" && -f "${CONDA_ENV_ACTIVATE}" ]] && source "${CONDA_ENV_ACTIVATE}"

RAW_FASTQ_DIR=$(dirname "$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f2 | tr -d '"')")

# ── Step 0: Pre-flight ────────────────────────────────────────────────────────
echo "=== [0] Pre-flight checks ==="
bash "${SCRIPT_DIR}/smoke_test.sh" "$SAMPLESHEET" "$F2T_CONFIG"

# ── Step 1: FastQC raw ────────────────────────────────────────────────────────
if is_done 1; then skip_msg 1; else
    echo "=== [1] FastQC (raw) ==="
    bash "${SCRIPT_DIR}/fastqc_batch.sh" \
        "$RAW_FASTQ_DIR" "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" "${THREADS_PARALLEL_JOBS}"
    multiqc "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" -n multiQC_unTrimmed \
        -o "${OUTPUT_DIR}/multiQC/multiQC_unTrimmed" --data-format tsv --export
    mark_done 1
fi

# ── Step 2: TrimGalore ────────────────────────────────────────────────────────
if is_done 2; then skip_msg 2; else
    echo "=== [2] TrimGalore ==="
    bash "${SCRIPT_DIR}/trimgalore_batch.sh" \
        "$SAMPLESHEET" "$RAW_FASTQ_DIR" "${OUTPUT_DIR}/trimmedFastq"
    mark_done 2
fi

# ── Step 3: FastQC trimmed ────────────────────────────────────────────────────
if is_done 3; then skip_msg 3; else
    echo "=== [3] FastQC (trimmed) ==="
    bash "${SCRIPT_DIR}/fastqc_batch.sh" \
        "${OUTPUT_DIR}/trimmedFastq" "${OUTPUT_DIR}/fastQC/fastQC_trimmed" "${THREADS_PARALLEL_JOBS}"
    multiqc "${OUTPUT_DIR}/fastQC/fastQC_trimmed" -n multiQC_trimmed \
        -o "${OUTPUT_DIR}/multiQC/multiQC_trimmed" --data-format tsv --export
    mark_done 3
fi

# ── Step 4: Alignment ─────────────────────────────────────────────────────────
if is_done 4; then skip_msg 4; else
    echo "=== [4] Bowtie2 alignment ==="
    bash "${SCRIPT_DIR}/bowtie2_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/trimmedFastq" "${OUTPUT_DIR}/bams"
    multiqc "${OUTPUT_DIR}/bams" -n multiQC_alignments \
        -o "${OUTPUT_DIR}/multiQC/multiQC_alignments" --data-format tsv --export
    mark_done 4
fi

# ── Step 5: Deduplication ─────────────────────────────────────────────────────
if is_done 5; then skip_msg 5; else
    echo "=== [5] Picard deduplication ==="
    bash "${SCRIPT_DIR}/picard_dedup_batch.sh" \
        "${OUTPUT_DIR}/bams" "${OUTPUT_DIR}/dedupBams" "${THREADS_PARALLEL_JOBS}"
    multiqc "${OUTPUT_DIR}/dedupBams" -n multiQC_deduplication \
        -o "${OUTPUT_DIR}/multiQC/multiQC_deduplication" --data-format tsv --export
    mark_done 5
fi

# ── Step 6: Blacklist filtering ───────────────────────────────────────────────
if is_done 6; then skip_msg 6; else
    echo "=== [6] Blacklist filtering ==="
    bash "${SCRIPT_DIR}/blacklist_filter_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/dedupBams" "${OUTPUT_DIR}/filteredBams"
    mark_done 6
fi

# ── Step 7: Genome coverage ───────────────────────────────────────────────────
if is_done 7; then skip_msg 7; else
    echo "=== [7] Genome coverage (individual) ==="
    bash "${SCRIPT_DIR}/genomecoverage_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig"
    mv "${OUTPUT_DIR}/bigwig/"*.bedGraph.gz       "${OUTPUT_DIR}/bedGraph/"     2>/dev/null || true
    mv "${OUTPUT_DIR}/bigwig/"*Snorm*.bedGraph.gz "${OUTPUT_DIR}/NormBedGraph/" 2>/dev/null || true
    mark_done 7
fi

# ── Step 8: Replicate merging ─────────────────────────────────────────────────
if is_done 8; then skip_msg 8; else
    echo "=== [8] Replicate merging + merged tracks ==="
    for GENOME in hg38 mm39; do
        grep -q ",$GENOME," "$SAMPLESHEET" && \
            bash "${SCRIPT_DIR}/merge_replicates.sh" \
                "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig_merged" "$GENOME"
    done
    mark_done 8
fi

# ── Step 9: MACS2 ─────────────────────────────────────────────────────────────
if is_done 9; then skip_msg 9; else
    echo "=== [9] MACS2 peak calling (narrow + broad) ==="
    bash "${SCRIPT_DIR}/macs2_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks"
    mark_done 9
fi

# ── Step 10: ChIPQC ───────────────────────────────────────────────────────────
if is_done 10; then skip_msg 10; else
    echo "=== [10] ChIPQC ==="
    for GENOME in hg38 mm39; do
        if grep -q ",$GENOME," "$SAMPLESHEET"; then
            "${R_BIN}" "${SCRIPT_DIR}/run_chipqc.R" \
                "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks" \
                "${OUTPUT_DIR}/chipqc" "$GENOME" "${THREADS_CHIPQC}" "narrow" "$F2T_CONFIG"
#             "${R_BIN}" "${SCRIPT_DIR}/run_chipqc.R" \
#                 "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks" \
#                 "${OUTPUT_DIR}/chipqc" "$GENOME" "${THREADS_CHIPQC}" "broad" "$F2T_CONFIG"
        fi
    done
    mark_done 10
fi

# ── Step 11: DiffBind prep ────────────────────────────────────────────────────
if is_done 11; then skip_msg 11; else
    echo "=== [11] DiffBind samplesheet preparation ==="
    for GENOME in hg38 mm39; do
        grep -q ",$GENOME," "$SAMPLESHEET" && \
            "${R_BIN}" "${SCRIPT_DIR}/prepare_diffbind.R" \
                "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks" \
                "${OUTPUT_DIR}/diffbind" "$GENOME"
    done
    mark_done 11
fi

# ── Step 12: UCSC tracks ──────────────────────────────────────────────────────
if is_done 12; then skip_msg 12; else
    echo "=== [12] UCSC tracks ==="
    bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
        "${OUTPUT_DIR}/bigwig"        "http://your-server.com/data"
    bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
        "${OUTPUT_DIR}/bigwig_merged" "http://your-server.com/data/merged"
    mark_done 12
fi

# ── Step 13: Report ───────────────────────────────────────────────────────────
if is_done 13; then skip_msg 13; else
    echo "=== [13] Pipeline report ==="
    bash "${SCRIPT_DIR}/generate_pipeline_report.sh" "$OUTPUT_DIR" \
        "${OUTPUT_DIR}/reports/pipeline_report_$(date +%Y%m%d)" html
    mark_done 13
fi

# ── Cleanup ───────────────────────────────────────────────────────────────────
[[ "${KEEP_INTERMEDIATE_BAMS}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/bams/"*.bam "${OUTPUT_DIR}/bams/"*.bai 2>/dev/null || true
[[ "${KEEP_TRIMMED_FASTQ}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/trimmedFastq/"*.gz 2>/dev/null || true

echo ""
echo "=== fastq2tracks v3.0.4 complete ==="
echo "  Checkpoints   : ${CHECKPOINT_DIR}/"
echo "  BigWig        : ${OUTPUT_DIR}/bigwig/"
echo "  BigWig merged : ${OUTPUT_DIR}/bigwig_merged/"
echo "  Peaks narrow  : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/narrow/"
echo "  Peaks broad   : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/broad/"
echo "  ChIPQC        : ${OUTPUT_DIR}/chipqc/"
echo "  DiffBind      : ${OUTPUT_DIR}/diffbind/"
echo "  Report        : ${OUTPUT_DIR}/reports/"
