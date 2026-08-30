#!/bin/bash
# =============================================================================
# ATACseq2tracks v4.3.0
# ATAC-seq / chromatin profiling track-generation and QC workflow
#
# Preferred usage after shared installation:
#   atacseq2tracks --config /path/to/config.conf
#
# Checkpoint system: completed steps are skipped on re-run.
# To rerun a named stage and everything after it:
#   atacseq2tracks --config /path/to/config.conf --from-stage peaks
# =============================================================================
set -euo pipefail

INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="${INSTALL_DIR}/scripts"
INPUT_SANITIZER="${SCRIPT_DIR}/sanitize_text_inputs.py"
CONFIG_RESOLVER="${SCRIPT_DIR}/resolve_config.py"
METADATA_WRITER="${SCRIPT_DIR}/prepare_run_metadata.py"
USER_CONFIG=""
F2T_CONFIG=""
PLAN_ONLY=false
PREFLIGHT_ONLY=false
FROM_STAGE=""
STOP_AFTER=""
PIPELINE_VERSION="$(tr -d '[:space:]' < "${INSTALL_DIR}/VERSION" 2>/dev/null || echo 4.3.0)"

declare -A STAGE_STEPS=(
    [fastqc_raw]=1 [trim]=2 [fastqc_trimmed]=3 [alignment]=4 [dedup]=5
    [filtering]=6 [spikein]=6s [coverage]=7 [merged_tracks]=8 [peaks]=9
    [qc]=10 [filtering_sensitivity]=10b [diffbind_prep]=11 [diffbind]=12
    [deseq2atac]=12a [browser]=13 [report]=14
)
STAGE_ORDER=(fastqc_raw trim fastqc_trimmed alignment dedup filtering spikein coverage merged_tracks peaks qc filtering_sensitivity diffbind_prep diffbind deseq2atac browser report)

usage() {
    cat <<'EOF'
Usage: atacseq2tracks --config /absolute/path/config.conf [options]

Options:
  --plan                 Validate inputs and write/print the execution plan only
  --preflight-only       Run complete pre-flight validation, then exit
  --from-stage NAME      Re-run NAME and every later stage
  --stop-after NAME      Stop successfully after NAME completes
  --version              Print the workflow version
  -h, --help             Show this help

Stage names: fastqc_raw, trim, fastqc_trimmed, alignment, dedup, filtering,
spikein, coverage, merged_tracks, peaks, qc, filtering_sensitivity,
diffbind_prep, diffbind, deseq2atac, browser, report.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) USER_CONFIG="$(realpath "${2:?missing value for --config}")"; shift 2 ;;
        --plan) PLAN_ONLY=true; shift ;;
        --preflight-only) PREFLIGHT_ONLY=true; shift ;;
        --from-stage) FROM_STAGE="${2:?missing value for --from-stage}"; shift 2 ;;
        --stop-after) STOP_AFTER="${2:?missing value for --stop-after}"; shift 2 ;;
        --version) echo "$PIPELINE_VERSION"; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$USER_CONFIG" ]] || { echo "ERROR: --config is required" >&2; usage >&2; exit 2; }
[[ -f "$USER_CONFIG" ]] || { echo "ERROR: config file not found: $USER_CONFIG" >&2; exit 1; }
[[ -f "$INPUT_SANITIZER" ]] || { echo "ERROR: input sanitizer not found: $INPUT_SANITIZER" >&2; exit 1; }
[[ -f "$CONFIG_RESOLVER" ]] || { echo "ERROR: config resolver not found: $CONFIG_RESOLVER" >&2; exit 1; }
[[ -f "$METADATA_WRITER" ]] || { echo "ERROR: metadata writer not found: $METADATA_WRITER" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required to validate input text files" >&2; exit 1; }

for requested_stage in "$FROM_STAGE" "$STOP_AFTER"; do
    [[ -z "$requested_stage" || -n "${STAGE_STEPS[$requested_stage]:-}" ]] || {
        echo "ERROR: unknown stage: $requested_stage" >&2
        usage >&2
        exit 2
    }
done

# Normalize the config before Bash sources it. Affected files are backed up
# beside the original; clean UTF-8/LF files are left byte-for-byte unchanged.
python3 "$INPUT_SANITIZER" "$USER_CONFIG"

RUNTIME_DIR="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-config.XXXXXX")"
trap 'status=$?; rm -rf -- "$RUNTIME_DIR"; exit "$status"' EXIT
python3 "$CONFIG_RESOLVER" \
    --config "$USER_CONFIG" --template "${INSTALL_DIR}/config/config.conf" \
    --shell-output "${RUNTIME_DIR}/resolved_config.conf" \
    --tsv-output "${RUNTIME_DIR}/resolved_config.tsv"
F2T_CONFIG="${RUNTIME_DIR}/resolved_config.conf"

export F2T_CONFIG
source "$F2T_CONFIG"

[[ -z "${SAMPLESHEET:-}" ]] && { echo "ERROR: SAMPLESHEET not set in $F2T_CONFIG" >&2; exit 1; }
[[ -z "${OUTPUT_DIR:-}"  ]] && { echo "ERROR: OUTPUT_DIR not set in $F2T_CONFIG"  >&2; exit 1; }
[[ -f "$SAMPLESHEET"     ]] || { echo "ERROR: SAMPLESHEET not found: $SAMPLESHEET" >&2; exit 1; }
python3 "$INPUT_SANITIZER" "$SAMPLESHEET"

mkdir -p "${OUTPUT_DIR}/metadata"
cp "${RUNTIME_DIR}/resolved_config.conf" "${OUTPUT_DIR}/metadata/resolved_config.conf"
cp "${RUNTIME_DIR}/resolved_config.tsv" "${OUTPUT_DIR}/metadata/resolved_config.tsv"
printf '%s\n' "$USER_CONFIG" > "${OUTPUT_DIR}/metadata/user_config_path.txt"
python3 "${SCRIPT_DIR}/validate_samplesheet.py" "$SAMPLESHEET"
python3 "$METADATA_WRITER" --samplesheet "$SAMPLESHEET" \
    --resolved-config-tsv "${OUTPUT_DIR}/metadata/resolved_config.tsv" \
    --output-dir "${OUTPUT_DIR}/metadata"

RESOURCE_EXCESS="$(awk -F '\t' 'NR > 1 && $8 == "yes" {print $1 ":" $6}' "${OUTPUT_DIR}/metadata/resource_budget.tsv" | paste -sd, -)"
if [[ -n "$RESOURCE_EXCESS" ]]; then
    case "${RESOURCE_CHECK_MODE:-warn}" in
        error) echo "ERROR: configured thread demand exceeds TOTAL_CPU_BUDGET: $RESOURCE_EXCESS" >&2; exit 1 ;;
        warn) echo "WARNING: configured thread demand exceeds TOTAL_CPU_BUDGET: $RESOURCE_EXCESS" >&2 ;;
        off) ;;
        *) echo "ERROR: RESOURCE_CHECK_MODE must be warn, error or off" >&2; exit 1 ;;
    esac
fi

export ATACSEQ2TRACKS_USER_CONFIG="$USER_CONFIG"

if [[ "$PLAN_ONLY" == "true" ]]; then
    echo "ATACseq2tracks v${PIPELINE_VERSION} execution plan"
    echo "Config: $USER_CONFIG"
    echo "Samplesheet: $SAMPLESHEET"
    echo "Output: $OUTPUT_DIR"
    column -t -s $'\t' "${OUTPUT_DIR}/metadata/planned_stages.tsv" 2>/dev/null || \
        cat "${OUTPUT_DIR}/metadata/planned_stages.tsv"
    echo "Technical-replicate audit: ${OUTPUT_DIR}/metadata/technical_merge_audit.tsv"
    echo "Resource budget: ${OUTPUT_DIR}/metadata/resource_budget.tsv"
    exit 0
fi

echo "============================================================"
echo " ATACseq2tracks v${PIPELINE_VERSION}"
echo " Install dir : $INSTALL_DIR"
echo " Config      : $USER_CONFIG"
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
WORKFLOW_SIGNATURE="$(find "$INSTALL_DIR" -maxdepth 2 -type f \
    \( -name '*.sh' -o -name '*.py' -o -name '*.R' -o -name 'VERSION' \) \
    -print0 | sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')"
RUN_SIGNATURE="$(sha256sum "$SAMPLESHEET" "$F2T_CONFIG" | { cat; printf '%s  workflow\n' "$WORKFLOW_SIGNATURE"; } | sha256sum | awk '{print $1}')"
EVENT_LOG="${OUTPUT_DIR}/metadata/workflow_events.tsv"
[[ -s "$EVENT_LOG" ]] || printf 'timestamp\tevent\tstep\tstage\tdetail\n' > "$EVENT_LOG"
CURRENT_STEP="0"
CURRENT_STAGE="initialization"
STAGE_STARTED="$(date +%s)"

stage_index() {
    local target="$1" index
    for index in "${!STAGE_ORDER[@]}"; do
        [[ "${STAGE_ORDER[$index]}" == "$target" ]] && { echo "$index"; return 0; }
    done
    return 1
}
stage_for_step() {
    local step="$1" stage
    for stage in "${!STAGE_STEPS[@]}"; do
        [[ "${STAGE_STEPS[$stage]}" == "$step" ]] && { echo "$stage"; return 0; }
    done
    echo "unknown"
}
log_event() {
    printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" "$2" "$3" "${4:-}" >> "$EVENT_LOG"
}
begin_stage() {
    CURRENT_STEP="$1"; CURRENT_STAGE="$2"; STAGE_STARTED="$(date +%s)"
    log_event start "$CURRENT_STEP" "$CURRENT_STAGE" ""
    echo "=== [${CURRENT_STEP}] $3 ==="
}
on_exit() {
    local status="$1"
    if (( status != 0 )); then
        log_event failure "$CURRENT_STEP" "$CURRENT_STAGE" "exit_status=${status}"
    fi
}
trap 'status=$?; on_exit "$status"; rm -rf -- "$RUNTIME_DIR"; exit "$status"' EXIT
log_event run_start 0 initialization "version=${PIPELINE_VERSION};signature=${RUN_SIGNATURE}"

is_done() {
    local step="$1" stage checkpoint_matches=false
    stage="$(stage_for_step "$step")"
    if [[ -f "${CHECKPOINT_DIR}/step${step}.done" ]] && [[ "$(cat "${CHECKPOINT_DIR}/step${step}.done")" == "$RUN_SIGNATURE" ]]; then
        checkpoint_matches=true
    fi
    if [[ -n "$FROM_STAGE" ]]; then
        if (( $(stage_index "$stage") < $(stage_index "$FROM_STAGE") )); then
            [[ "$checkpoint_matches" == "true" ]] || {
                echo "ERROR: --from-stage ${FROM_STAGE} requires a matching earlier checkpoint for ${stage}" >&2
                exit 1
            }
            return 0
        fi
        return 1
    fi
    [[ "$checkpoint_matches" == "true" ]]
}
mark_done() {
    printf '%s\n' "$RUN_SIGNATURE" > "${CHECKPOINT_DIR}/step${1}.done"
    local stage elapsed
    stage="$(stage_for_step "$1")"
    elapsed=$(( $(date +%s) - STAGE_STARTED ))
    log_event complete "$1" "$stage" "elapsed_seconds=${elapsed}"
    printf '{"step":"%s","stage":"%s","signature":"%s","completed_utc":"%s","elapsed_seconds":%s}\n' \
        "$1" "$stage" "$RUN_SIGNATURE" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$elapsed" \
        > "${CHECKPOINT_DIR}/step${1}.json"
    echo "[CHECKPOINT] Step ${1} (${stage}) complete -- to re-run: atacseq2tracks --config ${USER_CONFIG} --from-stage ${stage}"
    if [[ "$STOP_AFTER" == "$stage" ]]; then
        log_event stop_after "$1" "$stage" "requested"
        echo "Requested stop after stage '${stage}'."
        exit 0
    fi
}
skip_msg()  {
    local stage
    stage="$(stage_for_step "$1")"
    log_event skip "$1" "$stage" checkpoint
    echo "=== [${1}] ${stage} SKIPPED (matching checkpoint) ==="
    if [[ "$STOP_AFTER" == "$stage" ]]; then
        log_event stop_after "$1" "$stage" "matching_checkpoint"
        echo "Requested stop boundary '${stage}' is already complete."
        exit 0
    fi
}

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
    "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/permissive" \
    "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/intermediate" \
    "${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/stringent" \
    "${OUTPUT_DIR}/coverage_filtering_sensitivity" \
    "${OUTPUT_DIR}/coverage_filtering_policy_bams" \
    "${OUTPUT_DIR}/bigwig_spikein/stringent" \
    "${OUTPUT_DIR}/bigwig_spikein/dm6_control" \
    "${OUTPUT_DIR}/spikein/composite_bams" \
    "${OUTPUT_DIR}/spikein/dedup_bams" \
    "${OUTPUT_DIR}/spikein/stringent_host_bams" \
    "${OUTPUT_DIR}/spikein/stringent_dm6_bams" \
    "${OUTPUT_DIR}/spikein/logs" \
    "${OUTPUT_DIR}/spikein/tables" \
    "${OUTPUT_DIR}/bigwig_merged" \
    "${OUTPUT_DIR}/peaks/per_replicate" \
    "${OUTPUT_DIR}/peaks/pooled" \
    "${OUTPUT_DIR}/chipqc" \
    "${OUTPUT_DIR}/qc_post_alignment" \
    "${OUTPUT_DIR}/diffbind" \
    "${OUTPUT_DIR}/diffbind_results" \
    "${OUTPUT_DIR}/deseq2atac" \
    "${OUTPUT_DIR}/reports"

RAW_FASTQ_DIR=$(dirname "$(tail -n +2 "$SAMPLESHEET" | head -1 | cut -d',' -f2 | tr -d '"')")

# ── Step 0: Pre-flight ────────────────────────────────────────────────────────
begin_stage 0 preflight "Pre-flight checks"
bash "${SCRIPT_DIR}/smoke_test.sh" "$SAMPLESHEET" "$F2T_CONFIG"
bash "${SCRIPT_DIR}/capture_provenance.sh" "${OUTPUT_DIR}/metadata" \
    "${OUTPUT_DIR}/metadata/resolved_config.tsv"
log_event complete 0 preflight ""
if [[ "$PREFLIGHT_ONLY" == "true" ]]; then
    echo "Pre-flight validation completed successfully; workflow stages were not launched."
    exit 0
fi

# ── Step 1: FastQC raw ────────────────────────────────────────────────────────
if is_done 1; then skip_msg 1; else
    begin_stage 1 fastqc_raw "FastQC (raw)"
    bash "${SCRIPT_DIR}/fastqc_batch.sh" \
        "$RAW_FASTQ_DIR" "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" "${THREADS_PARALLEL_JOBS}" samplesheet
    multiqc "${OUTPUT_DIR}/fastQC/fastQC_unTrimmed" -n multiQC_unTrimmed \
        -o "${OUTPUT_DIR}/multiQC/multiQC_unTrimmed" --data-format tsv --export
    mark_done 1
fi

# ── Step 2: TrimGalore ────────────────────────────────────────────────────────
if is_done 2; then skip_msg 2; else
    begin_stage 2 trim "TrimGalore"
    bash "${SCRIPT_DIR}/trimgalore_batch.sh" \
        "$SAMPLESHEET" "$RAW_FASTQ_DIR" "${OUTPUT_DIR}/trimmedFastq"
    mark_done 2
fi

# ── Step 3: FastQC trimmed ────────────────────────────────────────────────────
if is_done 3; then skip_msg 3; else
    begin_stage 3 fastqc_trimmed "FastQC (trimmed)"
    bash "${SCRIPT_DIR}/fastqc_batch.sh" \
        "${OUTPUT_DIR}/trimmedFastq" "${OUTPUT_DIR}/fastQC/fastQC_trimmed" "${THREADS_PARALLEL_JOBS}" directory
    multiqc "${OUTPUT_DIR}/fastQC/fastQC_trimmed" -n multiQC_trimmed \
        -o "${OUTPUT_DIR}/multiQC/multiQC_trimmed" --data-format tsv --export
    mark_done 3
fi

# ── Step 4: Alignment ─────────────────────────────────────────────────────────
if is_done 4; then skip_msg 4; else
    begin_stage 4 alignment "Bowtie2 alignment"
    bash "${SCRIPT_DIR}/bowtie2_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/trimmedFastq" "${OUTPUT_DIR}/bams"
    multiqc "${OUTPUT_DIR}/bams" -n multiQC_alignments \
        -o "${OUTPUT_DIR}/multiQC/multiQC_alignments" --data-format tsv --export
    mark_done 4
fi

# ── Step 5: Deduplication ─────────────────────────────────────────────────────
if is_done 5; then skip_msg 5; else
    begin_stage 5 dedup "Picard deduplication"
    bash "${SCRIPT_DIR}/picard_dedup_batch.sh" \
        "${OUTPUT_DIR}/bams" "${OUTPUT_DIR}/dedupBams" "${THREADS_PARALLEL_JOBS}"
    multiqc "${OUTPUT_DIR}/dedupBams" "${OUTPUT_DIR}/logs/picard" -n multiQC_deduplication \
        -o "${OUTPUT_DIR}/multiQC/multiQC_deduplication" --data-format tsv --export
    mark_done 5
fi

# ── Step 6: Blacklist filtering ───────────────────────────────────────────────
if is_done 6; then skip_msg 6; else
    begin_stage 6 filtering "Blacklist filtering"
    bash "${SCRIPT_DIR}/blacklist_filter_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/dedupBams" "${OUTPUT_DIR}/filteredBams"
    mark_done 6
fi

# Step 6s is independent of the existing host-only branches. It uses a
# competitive host+dm6 composite alignment and creates only spike-in tracks.
if [[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
    if is_done 6s; then skip_msg 6s; else
        begin_stage 6s spikein "Drosophila spike-in stringent coverage"
        bash "${SCRIPT_DIR}/drosophila_spikein_tracks.sh" \
            "$SAMPLESHEET" "${OUTPUT_DIR}/trimmedFastq" "$OUTPUT_DIR"
        mark_done 6s
    fi
else
    echo "=== [6s] Drosophila spike-in coverage disabled ==="
fi

# ── Step 7: Genome coverage ───────────────────────────────────────────────────
if is_done 7; then skip_msg 7; else
    begin_stage 7 coverage "Fragment/read CPM coverage (individual)"
    if [[ "${GENERATE_CPM_TRACKS:-true}" == "true" ]]; then
        bash "${SCRIPT_DIR}/genomecoverage_batch.sh" \
            "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig"
    else
        echo "=== [7] CPM coverage disabled by GENERATE_CPM_TRACKS=false ==="
    fi
    mark_done 7
fi

# ── Step 8: Replicate merging ─────────────────────────────────────────────────
if is_done 8; then skip_msg 8; else
    begin_stage 8 merged_tracks "Replicate merging + merged CPM tracks"
    if [[ "${GENERATE_CPM_TRACKS:-true}" == "true" ]]; then
        for GENOME in hg38 mm39; do
            grep -q ",$GENOME," "$SAMPLESHEET" && \
                bash "${SCRIPT_DIR}/merge_replicates.sh" \
                    "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/bigwig_merged" "$GENOME"
        done
    else
        echo "=== [8] Merged CPM coverage disabled with CPM tracks ==="
    fi
    mark_done 8
fi

# ── Step 9: MACS3 (legacy script name retained) ──────────────────────────────
if is_done 9; then skip_msg 9; else
    begin_stage 9 peaks "MACS3 peak calling (sample-sheet mode)"
    bash "${SCRIPT_DIR}/macs2_batch.sh" \
        "$SAMPLESHEET" "${OUTPUT_DIR}/filteredBams" "${OUTPUT_DIR}/peaks"
    mark_done 9
fi
# ── Step 10: Post-alignment QC (deepTools — replaces ChIPQC) ─────────────────
if is_done 10; then skip_msg 10; else
    begin_stage 10 qc "Post-alignment QC (deepTools)"
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
# Step 10b: additional robust-CPM policies use the fixed Step 10 consensus universe.
# Their own count matrices and DESeq2 factors isolate duplicate and MAPQ effects.
if [[ "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" ]]; then
    if is_done 10b; then skip_msg 10b; else
        begin_stage 10b filtering_sensitivity "Read-filtering sensitivity robust-CPM tracks"
        bash "${SCRIPT_DIR}/generate_filtering_sensitivity_tracks.sh" \
            "$SAMPLESHEET" "$OUTPUT_DIR" \
            "${OUTPUT_DIR}/qc_post_alignment/peak_sets/consensus_peaks.bed"
        mark_done 10b
    fi
else
    echo "=== [10b] Permissive/intermediate robust-CPM tracks disabled ==="
fi

# Step 11: DiffBind samplesheet preparation
if is_done 11; then skip_msg 11; else
    begin_stage 11 diffbind_prep "DiffBind samplesheet preparation"
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
    begin_stage 12 diffbind "DiffBind differential analysis"
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
        begin_stage 12a deseq2atac "DESeq2ATAC differential accessibility analysis"
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
    begin_stage 13 browser "UCSC tracks"
    if find "${OUTPUT_DIR}/bigwig" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig" "${UCSC_BIGDATA_URL_BASE:-}" "${UCSC_TRACK_PREFIX:-ATAC-seq} CPM"
    fi
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
    for POLICY in permissive intermediate stringent; do
        POLICY_DIR="${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/${POLICY}"
        if find "$POLICY_DIR" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
            POLICY_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/deseq2_robust_cpm/${POLICY}}"
            bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" "$POLICY_DIR" "$POLICY_URL" \
                "${UCSC_TRACK_PREFIX:-ATAC-seq} DESeq2 robust CPM ${POLICY}"
        fi
    done
    if find "${OUTPUT_DIR}/bigwig_spikein/stringent" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        SPIKEIN_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/spikein/stringent}"
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig_spikein/stringent" "$SPIKEIN_URL" \
            "${UCSC_TRACK_PREFIX:-ATAC-seq} dm6 spike-in stringent"
    fi
    if find "${OUTPUT_DIR}/bigwig_spikein/dm6_control" -maxdepth 1 -name '*.bw' -print -quit | grep -q .; then
        DM6_CONTROL_URL="${UCSC_BIGDATA_URL_BASE:+${UCSC_BIGDATA_URL_BASE%/}/spikein/dm6_control}"
        bash "${SCRIPT_DIR}/create_ucsc_tracks.sh" \
            "${OUTPUT_DIR}/bigwig_spikein/dm6_control" "$DM6_CONTROL_URL" \
            "${UCSC_TRACK_PREFIX:-ATAC-seq} dm6 calibration control"
        mv "${OUTPUT_DIR}/bigwig_spikein/dm6_control/ucsc_tracks.txt" \
            "${OUTPUT_DIR}/bigwig_spikein/dm6_control/ucsc_tracks_dm6.txt"
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
    begin_stage 14 report "Pipeline report"
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
    bash "${SCRIPT_DIR}/cleanup_intermediates.sh" "$OUTPUT_DIR"
fi

log_event run_complete 14 report "version=${PIPELINE_VERSION}"

echo ""
echo "=== ATACseq2tracks v${PIPELINE_VERSION} complete ==="
echo "  Checkpoints   : ${CHECKPOINT_DIR}/"
echo "  CPM coverage         : ${OUTPUT_DIR}/bigwig/"
echo "  DESeq2 tracks        : ${OUTPUT_DIR}/bigwig_deseq2_consensus/"
echo "  DESeq2 robust CPM    : ${OUTPUT_DIR}/bigwig_deseq2_robust_cpm/"
echo "  Filtering sensitivity: ${OUTPUT_DIR}/coverage_filtering_sensitivity/"
echo "  dm6 spike-in tracks  : ${OUTPUT_DIR}/bigwig_spikein/stringent/"
echo "  dm6 control tracks   : ${OUTPUT_DIR}/bigwig_spikein/dm6_control/"
echo "  dm6 spike-in QC      : ${OUTPUT_DIR}/spikein/tables/"
echo "  Merged CPM    : ${OUTPUT_DIR}/bigwig_merged/"
echo "  Peaks narrow  : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/narrow/"
echo "  Peaks broad   : ${OUTPUT_DIR}/peaks/per_replicate/<sample>/broad/"
echo "  QC deepTools     : ${OUTPUT_DIR}/qc_post_alignment/"
echo "  QC TSS/periodicity: ${OUTPUT_DIR}/qc_post_alignment/atac_qc/"
echo "  DiffBind samples : ${OUTPUT_DIR}/diffbind/"
echo "  DiffBind results : ${OUTPUT_DIR}/diffbind_results/"
echo "  DESeq2ATAC       : ${OUTPUT_DIR}/deseq2atac/"
echo "  Report           : ${OUTPUT_DIR}/reports/"
