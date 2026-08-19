#!/bin/bash
# =============================================================================
# ATACseq2tracks v4.0.0
# ATAC-seq / chromatin profiling track-generation and QC workflow
#
# Usage:
#   bash /path/to/ATACseq2tracks/atacseq2tracks.sh --config /path/to/config.conf
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
INPUT_SANITIZER="${SCRIPT_DIR}/sanitize_text_inputs.py"
F2T_CONFIG=""
PIPELINE_VERSION="$(tr -d '[:space:]' < "${INSTALL_DIR}/VERSION" 2>/dev/null || echo 4.0.0)"

usage() {
    echo "Usage: bash atacseq2tracks.sh --config /absolute/path/to/config.conf"
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
[[ -f "$INPUT_SANITIZER" ]] || { echo "ERROR: input sanitizer not found: $INPUT_SANITIZER" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required to validate input text files" >&2; exit 1; }

# Normalize the config before Bash sources it. Affected files are backed up
# beside the original; clean UTF-8/LF files are left byte-for-byte unchanged.
python3 "$INPUT_SANITIZER" "$F2T_CONFIG"

export F2T_CONFIG
source "$F2T_CONFIG"

[[ -z "${SAMPLESHEET:-}" ]] && { echo "ERROR: SAMPLESHEET not set in $F2T_CONFIG" >&2; exit 1; }
[[ -z "${OUTPUT_DIR:-}"  ]] && { echo "ERROR: OUTPUT_DIR not set in $F2T_CONFIG"  >&2; exit 1; }
[[ -f "$SAMPLESHEET"     ]] || { echo "ERROR: SAMPLESHEET not found: $SAMPLESHEET" >&2; exit 1; }
python3 "$INPUT_SANITIZER" "$SAMPLESHEET"

echo "============================================================"
echo " ATACseq2tracks v${PIPELINE_VERSION}"
echo " Install dir : $INSTALL_DIR"
echo " Config      : $F2T_CONFIG"
echo " Samplesheet : $SAMPLESHEET"
echo " Output      : $OUTPUT_DIR"
echo " Max jobs    : ${THREADS_PARALLEL_JOBS}"
echo " QC jobs     : ${QC_SAMPLE_PARALLEL_JOBS:-4}"
echo " ataqv jobs  : ${ATAQV_PARALLEL_JOBS:-4}"
echo " Track jobs  : ${TRACK_PARALLEL_JOBS:-2}"
echo " Pooled MACS : ${POOLED_MACS_PARALLEL_JOBS:-2}"
echo " Merge jobs  : ${MERGE_PARALLEL_JOBS:-2}"
echo "============================================================"

# ── Checkpoint helpers ────────────────────────────────────────────────────────
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "$CHECKPOINT_DIR"
RUN_SIGNATURE="$(sha256sum "$SAMPLESHEET" "$F2T_CONFIG" "${INSTALL_DIR}/VERSION" | sha256sum | awk '{print $1}')"

is_done()   { [[ -f "${CHECKPOINT_DIR}/step${1}.done" ]] && [[ "$(cat "${CHECKPOINT_DIR}/step${1}.done")" == "$RUN_SIGNATURE" ]]; }
mark_done() {
    printf '%s\n' "$RUN_SIGNATURE" > "${CHECKPOINT_DIR}/step${1}.done"
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
    "${OUTPUT_DIR}/bigwig_deseq2_consensus" \
    "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm" \
    "${OUTPUT_DIR}/bigwig_merged" \
    "${OUTPUT_DIR}/peaks/per_replicate" \
    "${OUTPUT_DIR}/peaks/pooled" \
    "${OUTPUT_DIR}/chipqc" \
    "${OUTPUT_DIR}/qc_post_alignment" \
    "${OUTPUT_DIR}/diffbind" \
    "${OUTPUT_DIR}/diffbind_results" \
    "${OUTPUT_DIR}/deseq2atac" \
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
        "$RAW_FASTQ_DIR" "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" "${THREADS_PARALLEL_JOBS}" samplesheet
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
        "${OUTPUT_DIR}/trimmedFastq" "${OUTPUT_DIR}/fastQC/fastQC_trimmed" "${THREADS_PARALLEL_JOBS}" directory
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
    multiqc "${OUTPUT_DIR}/dedupBams" "${OUTPUT_DIR}/logs/picard" -n multiQC_deduplication \
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
    echo "=== [7] Fragment/read CPM coverage (individual) ==="
    bash "${SCRIPT_DIR}/genomecoverage_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig"
    mark_done 7
fi

# ── Step 8: Replicate merging ─────────────────────────────────────────────────
if is_done 8; then skip_msg 8; else
    echo "=== [8] Replicate merging + merged CPM tracks ==="
    for GENOME in hg38 mm39; do
        grep -q ",$GENOME," "$SAMPLESHEET" && \
            bash "${SCRIPT_DIR}/merge_replicates.sh" \
                "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig_merged" "$GENOME"
    done
    mark_done 8
fi

# ── Step 9: MACS3 (legacy script name retained) ──────────────────────────────
if is_done 9; then skip_msg 9; else
    echo "=== [9] MACS3 peak calling (sample-sheet mode) ==="
    bash "${SCRIPT_DIR}/macs2_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks"
    mark_done 9
fi
# ── Step 10: Post-alignment QC (deepTools — replaces ChIPQC) ─────────────────
if is_done 10; then skip_msg 10; else
    echo "=== [10] Post-alignment QC (deepTools) ==="
    mkdir -p "${OUTPUT_DIR}/qc_post_alignment"
    bash "${SCRIPT_DIR}/post_alignment_qc_batch.sh" \
        "$SAMPLESHEET" \
        "${OUTPUT_DIR}/filteredBams" \
        "${OUTPUT_DIR}/peaks" \
        "${OUTPUT_DIR}/bigwig" \
        "${OUTPUT_DIR}/qc_post_alignment"
    if [[ "${RUN_ATAQV_QC:-true}" == "true" ]]; then
        echo "=== [10] ATAC-specific QC (TSS enrichment and fragment periodicity) ==="
        bash "${SCRIPT_DIR}/ataqv_qc_batch.sh" \
            "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks" \
            "${OUTPUT_DIR}/qc_post_alignment/atac_qc"
    fi
    RUN_GENOME="$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f5 | tr -d '"')"
    if [[ "${RUN_PEAK_ANNOTATION:-false}" == "true" || "${RUN_MOTIF_ENRICHMENT:-false}" == "true" ]]; then
        bash "${SCRIPT_DIR}/peak_interpretation.sh" \
            "${OUTPUT_DIR}/qc_post_alignment/peak_sets/consensus_peaks.bed" \
            "$RUN_GENOME" "${OUTPUT_DIR}/peak_interpretation"
    fi
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

# ── Step 12: DiffBind differential analysis ──────────────────────────────────
DIFFERENTIAL_ANALYSIS_FAILURES=0
if is_done 12; then skip_msg 12; else
    echo "=== [12] DiffBind differential analysis ==="
    if bash "${SCRIPT_DIR}/diffbind_analysis.sh" \
        "${OUTPUT_DIR}/diffbind" "${OUTPUT_DIR}/diffbind_results"; then
        mark_done 12
    else
        echo "ERROR: DiffBind failed; DESeq2ATAC will still be attempted" >&2
        DIFFERENTIAL_ANALYSIS_FAILURES=$((DIFFERENTIAL_ANALYSIS_FAILURES + 1))
    fi
fi

# Independent peer analysis. The alphanumeric checkpoint preserves the
# established Step 13/14 checkpoint names and resumes separately from DiffBind.
if [[ "${RUN_DESEQ2ATAC:-true}" == "true" ]]; then
    if is_done 12a; then skip_msg 12a; else
        echo "=== [12a] DESeq2ATAC differential accessibility analysis ==="
        if bash "${SCRIPT_DIR}/deseq2atac_analysis.sh" \
            "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks" \
            "${OUTPUT_DIR}/deseq2atac"; then
            mark_done 12a
        else
            echo "ERROR: DESeq2ATAC failed; existing DiffBind results are retained" >&2
            DIFFERENTIAL_ANALYSIS_FAILURES=$((DIFFERENTIAL_ANALYSIS_FAILURES + 1))
        fi
    fi
else
    echo "=== [12a] DESeq2ATAC disabled by RUN_DESEQ2ATAC=false ==="
fi

(( DIFFERENTIAL_ANALYSIS_FAILURES == 0 )) || {
    echo "WARNING: ${DIFFERENTIAL_ANALYSIS_FAILURES} differential-accessibility module(s) failed; continuing through reporting" >&2
    rm -f "${CHECKPOINT_DIR}/step14.done"
}

# ── Step 13: UCSC tracks ──────────────────────────────────────────────────────
if is_done 13; then skip_msg 13; else
    echo "=== [13] UCSC tracks ==="
    bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
        "${OUTPUT_DIR}/bigwig" "${UCSC_BIGDATA_URL_BASE:-}" "${UCSC_TRACK_PREFIX:-ATAC-seq} CPM"
    if find "${OUTPUT_DIR}/bigwig_deseq2_consensus" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        DESEQ2_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/deseq2_consensus}"
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig_deseq2_consensus" "$DESEQ2_URL" "${UCSC_TRACK_PREFIX:-ATAC-seq} DESeq2 consensus"
    fi
    if find "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        ROBUST_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/deseq2_robust_cpm}"
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm" "$ROBUST_URL" "${UCSC_TRACK_PREFIX:-ATAC-seq} DESeq2 robust CPM"
    fi
    if find "${OUTPUT_DIR}/bigwig_merged" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        MERGED_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/merged}"
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig_merged" "$MERGED_URL" "${UCSC_TRACK_PREFIX:-ATAC-seq} merged CPM"
    fi
    mark_done 13
fi

# ── Step 14: Report ───────────────────────────────────────────────────────────
if is_done 14; then skip_msg 14; else
    echo "=== [14] Pipeline report ==="
    bash "${SCRIPT_DIR}/generate_pipeline_report.sh" "$OUTPUT_DIR" \
        "${OUTPUT_DIR}/reports/pipeline_report_$(date +%Y%m%d)" html
    if (( DIFFERENTIAL_ANALYSIS_FAILURES == 0 )); then
        mark_done 14
    else
        echo "WARNING: Step 14 report was written but not checkpointed because differential analysis failed" >&2
    fi
fi

(( DIFFERENTIAL_ANALYSIS_FAILURES == 0 )) || {
    echo "ERROR: ${DIFFERENTIAL_ANALYSIS_FAILURES} differential-accessibility module(s) failed" >&2
    echo "ERROR: reports were generated and automatic cleanup was suppressed" >&2
    exit 1
}

# ── Cleanup ───────────────────────────────────────────────────────────────────
if [[ "${ENABLE_AUTOMATIC_CLEANUP:-true}" == "true" ]]; then
[[ "${KEEP_INTERMEDIATE_BAMS:-false}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/bams/"*.bam "${OUTPUT_DIR}/bams/"*.bai 2>/dev/null || true
[[ "${KEEP_TRIMMED_FASTQ:-false}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/trimmedFastq/"*.gz 2>/dev/null || true
[[ "${KEEP_DEDUP_BAMS:-false}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/dedupBams/"*.bam "${OUTPUT_DIR}/dedupBams/"*.bai 2>/dev/null || true
[[ "${KEEP_FILTERED_BAMS:-true}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/filteredBams/"*.bam "${OUTPUT_DIR}/filteredBams/"*.bai 2>/dev/null || true
[[ "${KEEP_RAW_BEDGRAPH:-false}" == "false" ]] && \
    rm -f "${OUTPUT_DIR}/bedGraph/"*.bedGraph.gz 2>/dev/null || true
fi

echo ""
echo "=== ATACseq2tracks v${PIPELINE_VERSION} complete ==="
echo "  Checkpoints   : ${CHECKPOINT_DIR}/"
echo "  BigWig CPM           : ${OUTPUT_DIR}/bigwig/"
echo "  DESeq2 tracks        : ${OUTPUT_DIR}/bigwig_deseq2_consensus/"
echo "  DESeq2 robust CPM    : ${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/"
echo "  BigWig merged : ${OUTPUT_DIR}/bigwig_merged/"
echo "  Peaks narrow  : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/narrow/"
echo "  Peaks broad   : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/broad/"
echo "  QC deepTools     : ${OUTPUT_DIR}/qc_post_alignment/"
echo "  QC TSS/periodicity: ${OUTPUT_DIR}/qc_post_alignment/atac_qc/"
echo "  DiffBind samples : ${OUTPUT_DIR}/diffbind/"
echo "  DiffBind results : ${OUTPUT_DIR}/diffbind_results/"
echo "  DESeq2ATAC       : ${OUTPUT_DIR}/deseq2atac/"
echo "  Report           : ${OUTPUT_DIR}/reports/"
