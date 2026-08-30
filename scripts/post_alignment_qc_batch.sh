#!/bin/bash
# =============================================================================
# ATACseq2tracks v4.1.0 - post-alignment QC module (deepTools-based)
# Replaces ChIPQC with a robust, crash-proof deepTools + samtools/bedtools QC.
#
# Compatible assays: ChIP-seq, CUT&RUN, CUT&Tag, ChIPmentation, ATAC-seq
#
# Features:
#   - Per-sample + multi-sample ChIPQC-style karyogram plots (100 kb bins)
#     via plot_chrom_coverage.py (must live alongside this script in scripts/)
#   - ATAC-seq mitochondrial read fraction QC metric
#   - FRiP per sample (narrow + broad) and over consensus peaks
#   - deepTools: plotFingerprint, multiBamSummary, plotCorrelation, plotPCA
#   - deepTools: computeMatrix, plotHeatmap, plotProfile over consensus peaks
#   - Robust: zero-peak samples stay in all BAM-level analyses; never dropped
#
# Usage (called by atacseq2tracks.sh as Step 10):
#   bash scripts/post_alignment_qc_batch.sh \
#       <SAMPLESHEET> <BAM_DIR> <PEAKS_DIR> <BIGWIG_DIR> <OUT_QC_DIR>
# =============================================================================
set -uo pipefail   # NOTE: -e intentionally omitted — errors handled per-step

# ── Config loading ─────────────────────────────────────────────────────────────
_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.conf not found. Export F2T_CONFIG or pass --config." >&2
            exit 1
        }
    fi
}
_load_config

# ── Arguments ──────────────────────────────────────────────────────────────────
SAMPLESHEET="${1:?Usage: post_alignment_qc_batch.sh <samplesheet> <bam_dir> <peaks_dir> <bigwig_dir> <out_dir>}"
BAM_DIR="${2:?BAM_DIR required}"
PEAKS_DIR="${3:?PEAKS_DIR required}"
BIGWIG_DIR="${4:?BIGWIG_DIR required}"
OUT_DIR="${5:?OUT_DIR required}"

THREADS="${THREADS_DEEPTOOLS:-8}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QC_HELPERS="${SCRIPT_DIR}/qc_table_helpers.sh"
[[ -f "$QC_HELPERS" ]] || { echo "ERROR: missing QC helper: $QC_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$QC_HELPERS"
TRACK_HELPERS="${SCRIPT_DIR}/track_normalization_helpers.sh"
[[ -f "$TRACK_HELPERS" ]] || { echo "ERROR: missing track helper: $TRACK_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$TRACK_HELPERS"
PARALLEL_HELPERS="${SCRIPT_DIR}/parallel_job_helpers.sh"
[[ -f "$PARALLEL_HELPERS" ]] || { echo "ERROR: missing parallel helper: $PARALLEL_HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$PARALLEL_HELPERS"
QC_JOBS="${QC_SAMPLE_PARALLEL_JOBS:-4}"
TRACK_JOBS="${TRACK_PARALLEL_JOBS:-2}"
parallel_require_positive_integer QC_SAMPLE_PARALLEL_JOBS "$QC_JOBS" || exit 1
parallel_require_positive_integer TRACK_PARALLEL_JOBS "$TRACK_JOBS" || exit 1

# ── Output directory structure ─────────────────────────────────────────────────
QC_TABLES="${OUT_DIR}/tables"
QC_PLOTS="${OUT_DIR}/plots"
QC_CHRPLOTS="${OUT_DIR}/plots/chromosome_coverage"
QC_MATRICES="${OUT_DIR}/matrices"
QC_LOGS="${OUT_DIR}/logs"
QC_PEAKS="${OUT_DIR}/peak_sets"
QC_DT="${OUT_DIR}/deeptools"
BIGWIG_PEAKNORM_DIR="$(dirname "$BIGWIG_DIR")/bigwig_deseq2_consensus"
BIGWIG_ROBUST_DIR="$(dirname "$BIGWIG_DIR")/bigwig_deseq2_robust_cpm"
BIGWIG_ROBUST_STRINGENT_DIR="${BIGWIG_ROBUST_DIR}/stringent"
mkdir -p "$QC_TABLES" "$QC_PLOTS" "$QC_CHRPLOTS" "$QC_MATRICES" \
         "$QC_LOGS" "$QC_PEAKS" "$QC_DT" "$BIGWIG_PEAKNORM_DIR" \
         "$BIGWIG_ROBUST_DIR" "$BIGWIG_ROBUST_STRINGENT_DIR"
QC_WORK_ROOT="$(mktemp -d "${OUT_DIR}/.parallel-workers.XXXXXX")"
trap 'rm -rf "$QC_WORK_ROOT"' EXIT

MAIN_LOG="${QC_LOGS}/post_qc_${TIMESTAMP}.log"
SUMMARY_TSV="${QC_TABLES}/qc_summary.tsv"
WARNINGS_TSV="${QC_TABLES}/qc_warnings.tsv"
PARALLEL_TIMING_TSV="${QC_TABLES}/parallel_job_timing.tsv"

# ── Logging ────────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $1" | tee -a "$MAIN_LOG"; }

log "=== Post-alignment QC module started ==="
log "Samplesheet : $SAMPLESHEET"
log "BAM dir     : $BAM_DIR"
log "Peaks dir   : $PEAKS_DIR"
log "BigWig dir  : $BIGWIG_DIR"
log "Output dir  : $OUT_DIR"
log "Threads     : $THREADS"
log "QC jobs     : $QC_JOBS"
log "Track jobs  : $TRACK_JOBS"

# ── Tool checks ────────────────────────────────────────────────────────────────
MISSING_TOOLS=()
for tool in samtools bedtools bamCoverage multiBamSummary multiBigwigSummary \
            plotCorrelation plotPCA plotFingerprint computeMatrix \
            plotHeatmap plotProfile; do
    command -v "$tool" &>/dev/null || MISSING_TOOLS+=("$tool")
done
[[ ${#MISSING_TOOLS[@]} -gt 0 ]] && \
    warn "Missing tools (some steps may skip): ${MISSING_TOOLS[*]}"

# Check for karyogram Python script
KARYOGRAM_PY="${SCRIPT_DIR}/plot_chrom_coverage.py"
if [[ ! -f "$KARYOGRAM_PY" ]]; then
    warn "plot_chrom_coverage.py not found at ${KARYOGRAM_PY} — karyogram plots will be skipped"
    KARYOGRAM_PY=""
fi

# Check matplotlib is available if the script is present
if [[ -n "$KARYOGRAM_PY" ]]; then
    python3 -c "import matplotlib, numpy" 2>/dev/null \
        || { warn "matplotlib/numpy not available — karyogram plots will be skipped"; KARYOGRAM_PY=""; }
fi

# ── TSV header ─────────────────────────────────────────────────────────────────
echo -e "sample_id\tassay\tbam_file\ttotal_reads\tduplication_pct\tbl_filtered_reads\tmito_reads\tmito_pct\tn_narrow_peaks\tn_broad_peaks\tpeak_width_median\tpeak_max_signal\tfrip_narrow\tfrip_broad\tqc_status\tnotes" \
    > "$SUMMARY_TSV"
echo -e "sample_id\twarning_type\tvalue\tthreshold\tmessage" > "$WARNINGS_TSV"
printf 'scope\tlabel\tstart_epoch\tend_epoch\telapsed_seconds\tparallel_jobs\tthreads_per_job\tstatus\n' \
    > "$PARALLEL_TIMING_TSV"

# ── Helpers ────────────────────────────────────────────────────────────────────
add_warning() {
    local sid="$1" wtype="$2" val="$3" thr="$4" msg="$5"
    echo -e "${sid}\t${wtype}\t${val}\t${thr}\t${msg}" >> "$WARNINGS_TSV"
    warn "$sid — $wtype: $msg"
}

generate_deseq2_tracks() {
    local bam="$1" key="$2" layout="$3" size_factor="$4" robust_scale="$5"
    local consensus_count_sum="$6" cohort_geometric_mean="$7" metadata_row="$8"
    local genome="${SAMPLE_GENOME[$key]:-hg38}" sample tmp_dir std_bam
    local signal_count signal_unit consensus_scale final_output temp_output mapq_bins consensus_sha
    local consensus_bw consensus_bg robust_bw robust_bg legacy_bw legacy_bg
    local -a std_chr=()
    sample=$(basename "$bam" .bam)

    for value in "$size_factor" "$robust_scale" "$consensus_count_sum" "$cohort_geometric_mean"; do
        awk -v value="$value" 'BEGIN {
            valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value > 0
            exit(valid ? 0 : 1)
        }' || { warn "Invalid normalization value for $key: $value"; return 1; }
    done
    [[ -f "$bam" ]] || { warn "BAM not found for DESeq2 tracks: $bam"; return 1; }
    [[ "${TRACK_STANDARD_CHROMS_ONLY:-true}" == "true" ]] || {
        warn "TRACK_STANDARD_CHROMS_ONLY=false is incompatible with the unified canonical track universe"
        return 1
    }

    layout="$(resolve_signal_layout "$bam" "$layout")" || return 1
    signal_unit="$(signal_unit_for_layout "$layout")"
    samtools index -@ "${THREADS_BIGWIG:-2}" "$bam" || { warn "BAM indexing failed for $key"; return 1; }
    mapfile -t std_chr < <(canonical_contigs_from_bam "$bam" "$genome")
    (( ${#std_chr[@]} > 0 )) || { warn "No canonical chromosomes found in $bam for $genome"; return 1; }

    tmp_dir="$(mktemp -d "${BIGWIG_PEAKNORM_DIR}/.${sample}.tracks.XXXXXX")"
    std_bam="${tmp_dir}/${sample}.canonical.bam"
    if ! samtools view -@ "${THREADS_BIGWIG:-2}" -b -o "$std_bam" "$bam" "${std_chr[@]}"; then
        warn "Canonical-chromosome BAM creation failed for $key"
        rm -rf "$tmp_dir"
        return 1
    fi
    samtools index -@ "${THREADS_BIGWIG:-2}" "$std_bam" || {
        warn "Canonical-chromosome BAM indexing failed for $key"
        rm -rf "$tmp_dir"
        return 1
    }
    signal_count="$(signal_count_for_bam "$std_bam" "$layout" 2>/dev/null || echo 0)"
    if [[ "$signal_count" -le 0 ]]; then
        warn "Zero $signal_unit records in canonical BAM for $key"
        rm -rf "$tmp_dir"
        return 1
    fi
    mapq_bins="$(signal_mapq_bin_counts "$std_bam" "$layout")"
    consensus_sha="$(sha256sum "$CONSENSUS_PEAK" | awk '{print $1}')"

    consensus_scale=$(awk "BEGIN {printf \"%.12g\", 1/$size_factor}")
    write_track() {
        local destination="$1" format="$2" scale="$3" label="$4"
        temp_output="${tmp_dir}/$(basename "$destination")"
        if ! write_scaled_coverage_track "$std_bam" "$temp_output" "$format" "$scale" "$layout" \
            "${TRACK_BIN_SIZE:-10}" "${THREADS_BIGWIG:-2}"; then
            warn "bamCoverage failed for $label: $key"
            return 1
        fi
        mv -f "$temp_output" "$destination"
    }

    if [[ "${GENERATE_DESEQ2_CONSENSUS_TRACKS:-true}" == "true" ]]; then
        if [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" ]]; then
            consensus_bw="${BIGWIG_PEAKNORM_DIR}/${sample}_DESeq2Consensus.bw"
            write_track "$consensus_bw" bigwig "$consensus_scale" DESeq2-consensus || { rm -rf "$tmp_dir"; return 1; }
        fi
        if [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" ]]; then
            consensus_bg="${BIGWIG_PEAKNORM_DIR}/${sample}_DESeq2Consensus.bedGraph"
            write_track "$consensus_bg" bedgraph "$consensus_scale" DESeq2-consensus || { rm -rf "$tmp_dir"; return 1; }
        fi
    fi

    if [[ "${GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true}" == "true" ]]; then
        if [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" ]]; then
            robust_bw="${BIGWIG_ROBUST_STRINGENT_DIR}/${key}_DESeq2RobustCPM_Stringent.bw"
            legacy_bw="${BIGWIG_ROBUST_DIR}/${sample}_DESeq2RobustCPM.bw"
            write_track "$robust_bw" bigwig "$robust_scale" DESeq2-robust-CPM-stringent || { rm -rf "$tmp_dir"; return 1; }
            cp -f "$robust_bw" "$legacy_bw"
        fi
        if [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" ]]; then
            robust_bg="${BIGWIG_ROBUST_STRINGENT_DIR}/${key}_DESeq2RobustCPM_Stringent.bedGraph"
            legacy_bg="${BIGWIG_ROBUST_DIR}/${sample}_DESeq2RobustCPM.bedGraph"
            write_track "$robust_bg" bedgraph "$robust_scale" DESeq2-robust-CPM-stringent || { rm -rf "$tmp_dir"; return 1; }
            cp -f "$robust_bg" "$legacy_bg"
        fi
    fi
    rm -rf "$tmp_dir"
    printf '%s\tstringent\tfilteredBams\tremoved_by_Picard\tone_primary_default\t%s\texcluded\texcluded\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$key" "${MIN_MAPQ:-30}" "$genome" "$layout" "$signal_unit" "$signal_count" \
        $mapq_bins "$consensus_sha" "$consensus_count_sum" "$cohort_geometric_mean" \
        "$size_factor" "$consensus_scale" "$robust_scale" \
        > "$metadata_row"
    log "  Enabled DESeq2 consensus/stringent robust-CPM tracks generated for $key"
}

compute_frip() {
    local bam="$1" peak="$2"
    if [[ ! -f "$bam" || ! -f "$peak" || ! -s "$peak" ]]; then echo "NA"; return; fi
    local total rip
    total=$(samtools view -c -F 4 "$bam" 2>/dev/null || echo 0)
    [[ "$total" -eq 0 ]] && echo "NA" && return
    rip=$(bedtools sort -i "$peak" 2>/dev/null \
        | bedtools merge -i stdin 2>/dev/null \
        | bedtools intersect -u -a "$bam" -b stdin -ubam 2>/dev/null \
        | samtools view -c 2>/dev/null || echo 0)
    awk "BEGIN {printf \"%.4f\", $rip / $total}"
}

compute_mito_fraction() {
    local key="$1"
    # Prefer pre-blacklist-filter dedup BAM (still has chrM)
    local bam_mito="${BAM_DIR}/../dedupBams/${key}_dedup.bam"
    [[ ! -f "$bam_mito" ]] && bam_mito="$2"  # fall back to blFilt BAM arg $2
    if [[ ! -f "$bam_mito" ]]; then echo "NA NA"; return; fi
    local mito total
    mito=$(samtools idxstats "$bam_mito" 2>/dev/null \
        | awk '$1=="chrM"||$1=="MT"||$1=="chrMT"{sum+=$3}END{print sum+0}')
    total=$(samtools idxstats "$bam_mito" 2>/dev/null \
        | awk '{sum+=$3}END{print sum+0}')
    [[ "$total" -eq 0 ]] && echo "NA NA" && return
    awk "BEGIN{printf \"%d %.2f\", $mito, ($mito/$total)*100}"
}

# ── Phase 1: Per-sample metrics ────────────────────────────────────────────────
log "=== Phase 1: Per-sample metrics ==="

declare -A seen_keys
declare -a SAMPLE_KEYS SAMPLE_BAMS SAMPLE_NARROW SAMPLE_BROAD SAMPLE_BIGWIGS
declare -A SAMPLE_GENOME SAMPLE_ASSAY SAMPLE_LAYOUT

while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment \
    cell_type rep tech_rep is_ctrl ctrl_id macs2_mode blacklist rest; do

    [[ "$sid" == "sample_id" ]] && continue
    sid="${sid//\"/}"; rep="${rep//\"/}"; is_ctrl="${is_ctrl//\"/}"
    assay="${assay//\"/}"; genome="${genome//\"/}"; layout="${layout//\"/}"
    [[ "${is_ctrl,,}" == "true" || "$is_ctrl" == "1" ]] && continue

    KEY="${sid}_bioR${rep}"
    [[ -n "${seen_keys[$KEY]+x}" ]] && continue
    seen_keys["$KEY"]=1

    BAM="${BAM_DIR}/${KEY}_dedup_blFilt.bam"
    if [[ ! -f "$BAM" ]]; then
        warn "$KEY: required filtered BAM not found"
        exit 1
    fi

    SAMPLE_GENOME["$KEY"]="$genome"
    SAMPLE_ASSAY["$KEY"]="${assay,,}"
    SAMPLE_LAYOUT["$KEY"]="${layout^^}"
    NARROW_PEAK="${PEAKS_DIR}/per_replicate/${KEY}/narrow/${KEY}_peaks.narrowPeak"
    BROAD_PEAK="${PEAKS_DIR}/per_replicate/${KEY}/broad/${KEY}_peaks.broadPeak"
    BW="${BIGWIG_DIR}/${KEY}_dedup_blFilt_CPM.bw"
    [[ ! -f "$BW" ]] && BW="${BIGWIG_DIR}/${KEY}_CPM.bw"

    SAMPLE_KEYS+=("$KEY")
    SAMPLE_BAMS+=("$BAM")
    SAMPLE_NARROW+=("$NARROW_PEAK")
    SAMPLE_BROAD+=("$BROAD_PEAK")
    SAMPLE_BIGWIGS+=("$BW")

done < <(tail -n +2 "$SAMPLESHEET")

PHASE1_DIR="${QC_WORK_ROOT}/phase1"
mkdir -p "$PHASE1_DIR"

metric_worker_warning() {
    local file="$1" sid="$2" wtype="$3" val="$4" threshold="$5" message="$6"
    printf '%s\t%s\t%s\t%s\t%s\n' "$sid" "$wtype" "$val" "$threshold" "$message" >> "$file"
}

run_metric_worker_impl() {
    local i="$1" key bam assay narrow broad result warnings info
    local dup_metrics dup_pct total_reads dup_int mito_reads mito_pct mito_int
    local n_narrow=0 n_broad=0 peak_width_med=NA peak_max_sig=NA
    local frip_narrow=NA frip_broad=NA frip_pct qc_status
    local -a mito_info
    key="${SAMPLE_KEYS[$i]}"; bam="${SAMPLE_BAMS[$i]}"; assay="${SAMPLE_ASSAY[$key]}"
    narrow="${SAMPLE_NARROW[$i]}"; broad="${SAMPLE_BROAD[$i]}"
    result="${PHASE1_DIR}/${i}.summary"; warnings="${PHASE1_DIR}/${i}.warnings"; info="${PHASE1_DIR}/${i}.info"
    : > "$warnings"

    dup_metrics=$(find "${OUT_DIR}/../logs/picard" "${BAM_DIR}/../logs/picard" \
        "${BAM_DIR}/../dedupBams" -name "${key}_dup_metrics.txt" 2>/dev/null | head -1 || true)
    dup_pct=NA
    [[ -f "$dup_metrics" ]] && dup_pct=$(picard_duplication_pct "$dup_metrics" 2>/dev/null || echo NA)
    total_reads=$(samtools view -c -F 4 "$bam" 2>/dev/null) || return 1
    [[ "$total_reads" -lt 1000000 ]] && metric_worker_warning "$warnings" "$key" LOW_READS \
        "$total_reads" 1000000 "Very low aligned reads (<1M)"
    if [[ "$dup_pct" != NA ]]; then
        dup_int=${dup_pct%.*}
        [[ "$dup_int" -gt 80 ]] && metric_worker_warning "$warnings" "$key" HIGH_DUPLICATION \
            "$dup_pct" 80 "Duplication >80%"
    fi

    mito_info=($(compute_mito_fraction "$key" "$bam"))
    mito_reads="${mito_info[0]:-NA}"; mito_pct="${mito_info[1]:-NA}"
    if [[ "$assay" =~ ^atac && "$mito_pct" != NA ]]; then
        mito_int=${mito_pct%.*}
        [[ "$mito_int" -ge 25 ]] && metric_worker_warning "$warnings" "$key" HIGH_MITO_ATAC \
            "$mito_pct" 25 "ATAC-seq mitochondrial reads >=25%"
        [[ "$mito_int" -ge 10 && "$mito_int" -lt 25 ]] && metric_worker_warning "$warnings" "$key" ELEVATED_MITO_ATAC \
            "$mito_pct" 10 "ATAC-seq mitochondrial reads >=10%"
    fi

    if [[ -s "$narrow" ]]; then
        n_narrow=$(wc -l < "$narrow")
        peak_width_med=$(awk '{print $3-$2}' "$narrow" | sort -n \
            | awk '{a[NR]=$0}END{print (NR%2==0)?(a[NR/2]+a[NR/2+1])/2:a[int(NR/2)+1]}') || return 1
        peak_max_sig=$(awk 'BEGIN{max=0}{if($7>max)max=$7}END{print max}' "$narrow") || return 1
    fi
    [[ -s "$broad" ]] && n_broad=$(wc -l < "$broad")
    if [[ "$n_narrow" -eq 0 && "$n_broad" -eq 0 ]]; then
        metric_worker_warning "$warnings" "$key" ZERO_PEAKS 0 1 "No peaks called (narrow or broad)"
        qc_status=NO_PEAKS
    elif [[ "$n_narrow" -lt 200 && "$n_broad" -lt 200 ]]; then
        metric_worker_warning "$warnings" "$key" FEW_PEAKS "$n_narrow" 200 "Very few peaks (<200)"
        qc_status=POOR
    else
        qc_status=PASS
    fi

    if [[ "$n_narrow" -gt 0 ]]; then
        frip_narrow=$(compute_frip "$bam" "$narrow") || return 1
        if [[ "$frip_narrow" != NA ]]; then
            frip_pct=$(awk "BEGIN{printf \"%.0f\", $frip_narrow * 100}")
            [[ "$frip_pct" -lt 1 ]] && metric_worker_warning "$warnings" "$key" LOW_FRIP_NARROW \
                "$frip_narrow" 0.01 "FRiP (narrow) <1%"
        fi
    fi
    [[ "$n_broad" -gt 0 ]] && frip_broad=$(compute_frip "$bam" "$broad")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' \
        "$key" "$assay" "$bam" "$total_reads" "$dup_pct" "$total_reads" "$mito_reads" "$mito_pct" \
        "$n_narrow" "$n_broad" "$peak_width_med" "$peak_max_sig" "$frip_narrow" "$frip_broad" "$qc_status" \
        > "$result"
    printf '%s: reads=%s dup=%s mito=%s%% narrow=%s broad=%s frip_n=%s\n' \
        "$key" "$total_reads" "$dup_pct" "$mito_pct" "$n_narrow" "$n_broad" "$frip_narrow" > "$info"
}

run_metric_worker() {
    local i="$1" key start end status=0
    key="${SAMPLE_KEYS[$i]}"; start="$(date +%s)"
    run_metric_worker_impl "$i" > "${PHASE1_DIR}/${i}.worker.log" 2>&1 || status=$?
    end="$(date +%s)"
    parallel_write_timing_row "${PHASE1_DIR}/${i}.timing" qc_metrics "$key" \
        "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" "$QC_JOBS" 1
    return "$status"
}

parallel_pool_init "$QC_JOBS"
for i in "${!SAMPLE_KEYS[@]}"; do
    log "QUEUE metrics: ${SAMPLE_KEYS[$i]}"
    parallel_pool_submit "${SAMPLE_KEYS[$i]}" run_metric_worker "$i"
done
phase1_status=0
parallel_pool_wait_all || phase1_status=$?
(( phase1_status == 0 )) || { warn "Per-sample QC job(s) failed: $(parallel_failed_labels_csv)"; exit 1; }

for i in "${!SAMPLE_KEYS[@]}"; do
    cat "${PHASE1_DIR}/${i}.summary" >> "$SUMMARY_TSV"
    while IFS=$'\t' read -r warning_sid warning_type warning_value warning_threshold warning_message; do
        [[ -n "$warning_sid" ]] || continue
        add_warning "$warning_sid" "$warning_type" "$warning_value" "$warning_threshold" "$warning_message"
    done < "${PHASE1_DIR}/${i}.warnings"
    log "  $(cat "${PHASE1_DIR}/${i}.info")"
    cat "${PHASE1_DIR}/${i}.timing" >> "$PARALLEL_TIMING_TSV"
done

N_SAMPLES=${#SAMPLE_KEYS[@]}
log "Collected $N_SAMPLES IP samples"

if [[ $N_SAMPLES -eq 0 ]]; then
    warn "No valid samples — skipping deepTools cross-sample steps"
    exit 0
fi

# One multiBamSummary invocation must use one signal definition. Mixing PE
# fragments and SE reads in the same normalization cohort is not comparable.
declare -A NORMALIZATION_LAYOUTS=()
for KEY in "${SAMPLE_KEYS[@]}"; do
    NORMALIZATION_LAYOUTS["${SAMPLE_LAYOUT[$KEY]}"]=1
done
if (( ${#NORMALIZATION_LAYOUTS[@]} != 1 )); then
    warn "Mixed PE/SE layouts are not supported in one DESeq2 normalization cohort"
    exit 1
fi
NORMALIZATION_LAYOUT="${SAMPLE_LAYOUT[${SAMPLE_KEYS[0]}]}"
declare -a PEAK_SIGNAL_ARGS=()
deeptools_signal_args "$NORMALIZATION_LAYOUT" PEAK_SIGNAL_ARGS || {
    warn "Unsupported normalization layout: $NORMALIZATION_LAYOUT"
    exit 1
}

# =============================================================================
# Phase 2: Chromosome-wide karyogram plots (ChIPQC-style)
# Strategy: bamCoverage --binSize 100000 → bedGraph per sample
#           → plot_chrom_coverage.py → per-sample + multi-sample PNG
# =============================================================================
log "=== Phase 2: Chromosome-wide karyogram plots (100 kb bins) ==="

KARYOGRAM_BG_LIST=()    # bedGraph paths for karyogram Python script
KARYOGRAM_LABELS=()     # corresponding labels

PHASE2_DIR="${QC_WORK_ROOT}/phase2"
mkdir -p "$PHASE2_DIR"

run_karyogram_worker_impl() {
    local i="$1" key bam genome bg_out chr_tsv karyo_tmp karyo_bam total
    local -a karyo_contigs=() karyo_signal_args=()
    key="${SAMPLE_KEYS[$i]}"; bam="${SAMPLE_BAMS[$i]}"; genome="${SAMPLE_GENOME[$key]:-hg38}"
    bg_out="${QC_CHRPLOTS}/${key}_100kb.bedGraph"
    chr_tsv="${QC_CHRPLOTS}/${key}_per_chrom.tsv"
    : > "${PHASE2_DIR}/${i}.warning"

    if [[ ! -s "$bg_out" ]]; then
        mapfile -t karyo_contigs < <(canonical_contigs_from_bam "$bam" "$genome")
        deeptools_signal_args "${SAMPLE_LAYOUT[$key]}" karyo_signal_args || {
            echo "unsupported layout for karyogram: $key" > "${PHASE2_DIR}/${i}.warning"; return 2;
        }
        if (( ${#karyo_contigs[@]} == 0 )); then
            echo "no canonical chromosomes for karyogram: $key" > "${PHASE2_DIR}/${i}.warning"
            return 2
        fi
        karyo_tmp="$(mktemp -d "${QC_CHRPLOTS}/.${key}.canonical.XXXXXX")"
        karyo_bam="${karyo_tmp}/${key}.canonical.bam"
        if ! samtools view -@ "$THREADS" -b -o "$karyo_bam" "$bam" "${karyo_contigs[@]}" \
            || ! samtools index -@ "$THREADS" "$karyo_bam"; then
            echo "canonical BAM preparation failed for karyogram: $key" > "${PHASE2_DIR}/${i}.warning"
            rm -rf "$karyo_tmp"; return 2
        fi
        if ! bamCoverage -b "$karyo_bam" --binSize 100000 --normalizeUsing RPKM \
            --skipNonCoveredRegions --outFileFormat bedgraph "${karyo_signal_args[@]}" \
            -p "$THREADS" -o "$bg_out"; then
            rm -f "$bg_out"
            echo "bamCoverage failed for karyogram: $key" > "${PHASE2_DIR}/${i}.warning"
            rm -rf "$karyo_tmp"; return 2
        fi
        rm -rf "$karyo_tmp"
    fi

    if [[ ! -s "$chr_tsv" ]]; then
        total=$(samtools view -c -F 4 "$bam" 2>/dev/null) || return 1
        {
            echo -e "chrom\tchrom_reads\ttotal_reads\treads_per_mb"
            samtools idxstats "$bam" 2>/dev/null \
                | awk -v tot="$total" '$1~/^chr[0-9XY]+$/ || $1=="chrM" || $1=="MT" {
                    rpm=($2>0)?($3/$2*1e6):0
                    printf "%s\t%d\t%d\t%.4f\n",$1,$3,tot,rpm
                }'
        } > "$chr_tsv"
    fi
}

run_karyogram_worker() {
    local i="$1" key start end status=0 reported_status
    key="${SAMPLE_KEYS[$i]}"; start="$(date +%s)"
    run_karyogram_worker_impl "$i" > "${PHASE2_DIR}/${i}.worker.log" 2>&1 || status=$?
    end="$(date +%s)"
    if (( status == 0 )); then reported_status=SUCCESS
    elif (( status == 2 )); then reported_status=WARNING
    else reported_status=FAILED; fi
    parallel_write_timing_row "${PHASE2_DIR}/${i}.timing" karyogram "$key" \
        "$start" "$end" "$reported_status" "$QC_JOBS" "$THREADS"
    (( status == 2 )) && return 0
    return "$status"
}

parallel_pool_init "$QC_JOBS"
for i in "${!SAMPLE_KEYS[@]}"; do
    log "QUEUE karyogram: ${SAMPLE_KEYS[$i]}"
    parallel_pool_submit "${SAMPLE_KEYS[$i]}" run_karyogram_worker "$i"
done
phase2_status=0
parallel_pool_wait_all || phase2_status=$?
(( phase2_status == 0 )) || { warn "Karyogram worker(s) failed: $(parallel_failed_labels_csv)"; exit 1; }
for i in "${!SAMPLE_KEYS[@]}"; do
    KEY="${SAMPLE_KEYS[$i]}"; BG_OUT="${QC_CHRPLOTS}/${KEY}_100kb.bedGraph"
    [[ -s "${PHASE2_DIR}/${i}.warning" ]] && warn "  $(cat "${PHASE2_DIR}/${i}.warning")"
    [[ -s "$BG_OUT" ]] && KARYOGRAM_BG_LIST+=("$BG_OUT") && KARYOGRAM_LABELS+=("$KEY")
    cat "${PHASE2_DIR}/${i}.timing" >> "$PARALLEL_TIMING_TSV"
done

# Consolidated per-chromosome table (all samples)
CHR_SUMMARY_TSV="${QC_TABLES}/per_chromosome_reads.tsv"
{
    echo -e "sample_id\tchrom\tchrom_reads\ttotal_reads\treads_per_mb"
    for KEY in "${SAMPLE_KEYS[@]}"; do
        tsv="${QC_CHRPLOTS}/${KEY}_per_chrom.tsv"
        [[ -f "$tsv" ]] && tail -n +2 "$tsv" | awk -v k="$KEY" '{print k"\t"$0}'
    done
} > "$CHR_SUMMARY_TSV"
log "Per-chrom table: $CHR_SUMMARY_TSV"

# Run karyogram Python plotter
if [[ -n "$KARYOGRAM_PY" && ${#KARYOGRAM_BG_LIST[@]} -gt 0 ]]; then
    FIRST_KEY="${SAMPLE_KEYS[0]}"
    FIRST_GENOME="${SAMPLE_GENOME[$FIRST_KEY]:-hg38}"
    if [[ "${FIRST_GENOME,,}" == "hg38" ]]; then
        CS="${CHROM_SIZES_HUMAN:-}"
    else
        CS="${CHROM_SIZES_MOUSE:-}"
    fi

    if [[ -n "$CS" && -f "$CS" ]]; then
        log "  Generating karyogram plots (${#KARYOGRAM_BG_LIST[@]} samples)..."
        python3 "$KARYOGRAM_PY" \
            --bedgraph "${KARYOGRAM_BG_LIST[@]}" \
            --labels   "${KARYOGRAM_LABELS[@]}" \
            --genome   "$FIRST_GENOME" \
            --chrom-sizes "$CS" \
            --outdir   "$QC_CHRPLOTS" \
            >> "$MAIN_LOG" 2>&1 \
        && log "  Karyogram plots complete → $QC_CHRPLOTS/" \
        || warn "  Karyogram plots FAILED (check log)"
    else
        warn "CHROM_SIZES not set or file missing — skipping karyogram plots"
    fi
else
    warn "Skipping karyogram: no bedGraphs generated or plot_chrom_coverage.py missing"
fi

# =============================================================================
# Phase 3: deepTools genome-wide QC
# =============================================================================
log "=== Phase 3: deepTools genome-wide QC ==="

BAM_LIST=("${SAMPLE_BAMS[@]}")
LABELS=("${SAMPLE_KEYS[@]}")

# plotFingerprint + multiBamSummary bins (parallel)
(
    multiBamSummary bins \
        -b "${BAM_LIST[@]}" \
        --labels "${LABELS[@]}" \
        "${PEAK_SIGNAL_ARGS[@]}" \
        -p "$THREADS" \
        -o "${QC_DT}/multiBamSummary_bins.npz" \
        --outRawCounts "${QC_MATRICES}/multiBamSummary_bins.tab" \
        >> "$MAIN_LOG" 2>&1 \
    && log "multiBamSummary (bins) complete" \
    || warn "multiBamSummary (bins) failed"
) &
FP_PID=$!

plotFingerprint \
    -b "${BAM_LIST[@]}" \
    --labels "${LABELS[@]}" \
    -p "$THREADS" \
    --minMappingQuality 20 \
    --skipZeros \
    --plotFile "${QC_PLOTS}/fingerprint.png" \
    --outRawCounts "${QC_MATRICES}/fingerprint.tab" \
    --outQualityMetrics "${QC_TABLES}/fingerprint_metrics.tsv" \
    >> "$MAIN_LOG" 2>&1 \
    && log "plotFingerprint complete" \
    || warn "plotFingerprint failed"

wait $FP_PID

if [[ -f "${QC_DT}/multiBamSummary_bins.npz" ]]; then
    for method in pearson spearman; do
        plotCorrelation \
            -in "${QC_DT}/multiBamSummary_bins.npz" \
            --corMethod "$method" --skipZeros \
            --whatToPlot heatmap --colorMap RdYlBu \
            --plotTitle "Sample ${method} correlation (genome-wide bins)" \
            --plotFile "${QC_PLOTS}/correlation_heatmap_${method}.png" \
            --outFileCorMatrix "${QC_MATRICES}/correlation_matrix_${method}.tsv" \
            >> "$MAIN_LOG" 2>&1 || warn "plotCorrelation ($method) failed"
    done
    plotPCA \
        -in "${QC_DT}/multiBamSummary_bins.npz" \
        --plotTitle "PCA — genome-wide bins" \
        --plotFile "${QC_PLOTS}/pca_bins.png" \
        --outFileNameData "${QC_MATRICES}/pca_bins_data.tsv" \
        >> "$MAIN_LOG" 2>&1 || warn "plotPCA (bins) failed"
fi

# =============================================================================
# Phase 4: Consensus peaks + peak-centric deepTools
# =============================================================================
log "=== Phase 4: Consensus peaks + peak-centric QC ==="

NARROW_WITH_PEAKS=(); BROAD_WITH_PEAKS=()
CANONICAL_PEAK_INPUT_DIR="${QC_PEAKS}/canonical_inputs"
mkdir -p "$CANONICAL_PEAK_INPUT_DIR"
for i in "${!SAMPLE_KEYS[@]}"; do
    KEY="${SAMPLE_KEYS[$i]}"; GENOME="${SAMPLE_GENOME[$KEY]}"
    pk="${SAMPLE_NARROW[$i]}"
    if [[ -f "$pk" && -s "$pk" ]]; then
        canonical_pk="${CANONICAL_PEAK_INPUT_DIR}/${KEY}.narrow.canonical.bed"
        canonicalize_peak_file "$pk" "$canonical_pk" "$GENOME" || { warn "Canonical narrow-peak filtering failed for $KEY"; exit 1; }
        [[ -s "$canonical_pk" ]] && NARROW_WITH_PEAKS+=("$canonical_pk")
    fi
    pk="${SAMPLE_BROAD[$i]}"
    if [[ -f "$pk" && -s "$pk" ]]; then
        canonical_pk="${CANONICAL_PEAK_INPUT_DIR}/${KEY}.broad.canonical.bed"
        canonicalize_peak_file "$pk" "$canonical_pk" "$GENOME" || { warn "Canonical broad-peak filtering failed for $KEY"; exit 1; }
        [[ -s "$canonical_pk" ]] && BROAD_WITH_PEAKS+=("$canonical_pk")
    fi
done

MERGED_NARROW="${QC_PEAKS}/merged_narrow.bed"
MERGED_BROAD="${QC_PEAKS}/merged_broad.bed"
CONSENSUS_PEAK="${QC_PEAKS}/consensus_peaks.bed"
CONSENSUS_SUPPORTED="${QC_PEAKS}/consensus_peaks_supported.bed"

[[ ${#NARROW_WITH_PEAKS[@]} -gt 0 ]] && \
    cat "${NARROW_WITH_PEAKS[@]}" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > "$MERGED_NARROW" 2>>"$MAIN_LOG" \
    && log "Merged narrow: $(wc -l < "$MERGED_NARROW") regions"

[[ ${#BROAD_WITH_PEAKS[@]} -gt 0 ]] && \
    cat "${BROAD_WITH_PEAKS[@]}" | cut -f1-3 | sort -k1,1 -k2,2n \
    | bedtools merge -i stdin > "$MERGED_BROAD" 2>>"$MAIN_LOG" \
    && log "Merged broad: $(wc -l < "$MERGED_BROAD") regions"

CONSENSUS_MIN_SAMPLES="${CONSENSUS_MIN_SAMPLES:-2}"
CONSENSUS_PEAK_TYPE="${CONSENSUS_PEAK_TYPE:-narrow}"
if [[ "$CONSENSUS_PEAK_TYPE" == "broad" ]]; then
    CONSENSUS_INPUTS=("${BROAD_WITH_PEAKS[@]}")
else
    CONSENSUS_INPUTS=("${NARROW_WITH_PEAKS[@]}")
fi

if (( ${#CONSENSUS_INPUTS[@]} >= CONSENSUS_MIN_SAMPLES )); then
    log "Building ${CONSENSUS_PEAK_TYPE} consensus: support >=${CONSENSUS_MIN_SAMPLES} biological samples"
    bedtools multiinter -i "${CONSENSUS_INPUTS[@]}" \
        | awk -v minimum="$CONSENSUS_MIN_SAMPLES" '$4 >= minimum {print $1"\t"$2"\t"$3}' \
        | bedtools sort -i stdin | bedtools merge -i stdin > "$CONSENSUS_SUPPORTED" 2>>"$MAIN_LOG"
    [[ -s "$CONSENSUS_SUPPORTED" ]] || { warn "No consensus peaks meet support >=${CONSENSUS_MIN_SAMPLES}"; rm -f "$CONSENSUS_SUPPORTED"; }
elif [[ "${ALLOW_SINGLE_SAMPLE_CONSENSUS:-false}" == "true" && ${#CONSENSUS_INPUTS[@]} -eq 1 ]]; then
    cut -f1-3 "${CONSENSUS_INPUTS[0]}" | bedtools sort -i stdin | bedtools merge -i stdin > "$CONSENSUS_SUPPORTED"
    warn "Single-sample consensus enabled explicitly; this is not replicate-supported"
else
    warn "Consensus not built: ${#CONSENSUS_INPUTS[@]} peak files, minimum=${CONSENSUS_MIN_SAMPLES}"
fi

if [[ -s "$CONSENSUS_SUPPORTED" ]]; then
    cp "$CONSENSUS_SUPPORTED" "$CONSENSUS_PEAK"
    log "Consensus peaks: $(wc -l < "$CONSENSUS_PEAK") regions"
fi

CURRENT_DESEQ2_TRACKS_REQUIRED=false
if [[ "${GENERATE_DESEQ2_CONSENSUS_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true}" == "true" ]]; then
    CURRENT_DESEQ2_TRACKS_REQUIRED=true
fi
ANY_DESEQ2_TRACKS_REQUIRED="$CURRENT_DESEQ2_TRACKS_REQUIRED"
if [[ "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" ]]; then
    ANY_DESEQ2_TRACKS_REQUIRED=true
fi

if [[ "$ANY_DESEQ2_TRACKS_REQUIRED" == "true" && ! -s "$CONSENSUS_PEAK" ]]; then
    warn "Enabled DESeq2 coverage families require a non-empty consensus peak set"
    exit 1
fi

if [[ -f "$CONSENSUS_PEAK" && -s "$CONSENSUS_PEAK" ]]; then
    multiBamSummary BED-file \
        --BED "$CONSENSUS_PEAK" \
        -b "${BAM_LIST[@]}" \
        --labels "${LABELS[@]}" \
        "${PEAK_SIGNAL_ARGS[@]}" \
        -p "$THREADS" \
        -o "${QC_DT}/multiBamSummary_peaks.npz" \
        --outRawCounts "${QC_MATRICES}/multiBamSummary_peaks.tab" \
        >> "$MAIN_LOG" 2>&1 \
        && log "multiBamSummary (peaks) complete" \
        || { warn "multiBamSummary (peaks) failed"; exit 1; }

    if [[ -f "${QC_DT}/multiBamSummary_peaks.npz" ]]; then
        plotCorrelation \
            -in "${QC_DT}/multiBamSummary_peaks.npz" \
            --corMethod pearson --skipZeros \
            --whatToPlot heatmap --colorMap RdYlBu \
            --plotTitle "Sample correlation (consensus peaks)" \
            --plotFile "${QC_PLOTS}/correlation_heatmap_pearson_peaks.png" \
            --outFileCorMatrix "${QC_MATRICES}/correlation_matrix_pearson_peaks.tsv" \
            >> "$MAIN_LOG" 2>&1 || warn "plotCorrelation (peaks) failed"

        plotPCA \
            -in "${QC_DT}/multiBamSummary_peaks.npz" \
            --plotTitle "PCA — consensus peaks" \
            --plotFile "${QC_PLOTS}/pca_peaks.png" \
            --outFileNameData "${QC_MATRICES}/pca_peaks_data.tsv" \
            >> "$MAIN_LOG" 2>&1 || warn "plotPCA (peaks) failed"
    fi
fi

# computeMatrix + plotHeatmap + plotProfile
BW_EXIST=(); BW_LABELS=()
for i in "${!SAMPLE_KEYS[@]}"; do
    bw="${SAMPLE_BIGWIGS[$i]}"
    [[ -f "$bw" ]] && BW_EXIST+=("$bw") && BW_LABELS+=("${SAMPLE_KEYS[$i]}")
done

if [[ ${#BW_EXIST[@]} -gt 0 && -f "$CONSENSUS_PEAK" && -s "$CONSENSUS_PEAK" ]]; then
    MATRIX="${QC_MATRICES}/signal_over_peaks.gz"
    computeMatrix reference-point \
        -S "${BW_EXIST[@]}" \
        -R "$CONSENSUS_PEAK" \
        --referencePoint center \
        -b 2000 -a 2000 \
        --skipZeros \
        --samplesLabel "${BW_LABELS[@]}" \
        -p "$THREADS" \
        -o "$MATRIX" \
        >> "$MAIN_LOG" 2>&1 && log "computeMatrix complete" || warn "computeMatrix failed"

    if [[ -f "$MATRIX" ]]; then
        plotHeatmap \
            -m "$MATRIX" \
            --colorMap RdYlBu_r \
            --whatToShow "heatmap and colorbar" \
            --plotTitle "Signal over consensus peaks (±2kb)" \
            --heatmapHeight 15 \
            -o "${QC_PLOTS}/heatmap_signal_over_peaks.png" \
            --outFileSortedRegions "${QC_PEAKS}/heatmap_sorted_regions.bed" \
            >> "$MAIN_LOG" 2>&1 || warn "plotHeatmap failed"

        plotProfile \
            -m "$MATRIX" \
            --plotTitle "Average signal over consensus peaks (±2kb)" \
            --perGroup \
            -o "${QC_PLOTS}/profile_signal_over_peaks.png" \
            --outFileNameData "${QC_MATRICES}/profile_signal_over_peaks.tsv" \
            >> "$MAIN_LOG" 2>&1 || warn "plotProfile failed"
    fi
fi

# FRiP over consensus peaks
if [[ -f "$CONSENSUS_PEAK" && -s "$CONSENSUS_PEAK" ]]; then
    FRIP_CONSENSUS_TSV="${QC_TABLES}/frip_consensus.tsv"
    echo -e "sample_id\ttotal_reads\treads_in_consensus_peaks\tfrip_consensus" > "$FRIP_CONSENSUS_TSV"
    FRIP_WORK_DIR="${QC_WORK_ROOT}/frip_consensus"
    mkdir -p "$FRIP_WORK_DIR"

    frip_consensus_worker_impl() {
        local i="$1" key bam total rip frip_c
        key="${SAMPLE_KEYS[$i]}"; bam="${SAMPLE_BAMS[$i]}"
        [[ -f "$bam" ]] || { printf '%s\tNA\tNA\tNA\n' "$key" > "${FRIP_WORK_DIR}/${i}.row"; return 0; }
        total=$(samtools view -c -F 4 "$bam" 2>/dev/null) || return 1
        rip=$(bedtools intersect -u -a "$bam" -b "$CONSENSUS_PEAK" -ubam 2>/dev/null \
            | samtools view -c 2>/dev/null) || return 1
        frip_c=$(awk "BEGIN{if($total>0){printf \"%.4f\",$rip/$total}else{print \"NA\"}}")
        printf '%s\t%s\t%s\t%s\n' "$key" "$total" "$rip" "$frip_c" > "${FRIP_WORK_DIR}/${i}.row"
    }

    frip_consensus_worker() {
        local i="$1" key start end status=0
        key="${SAMPLE_KEYS[$i]}"; start="$(date +%s)"
        frip_consensus_worker_impl "$i" > "${FRIP_WORK_DIR}/${i}.worker.log" 2>&1 || status=$?
        end="$(date +%s)"
        parallel_write_timing_row "${FRIP_WORK_DIR}/${i}.timing" consensus_frip "$key" \
            "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" "$QC_JOBS" 1
        return "$status"
    }

    parallel_pool_init "$QC_JOBS"
    for i in "${!SAMPLE_KEYS[@]}"; do
        log "QUEUE consensus FRiP: ${SAMPLE_KEYS[$i]}"
        parallel_pool_submit "${SAMPLE_KEYS[$i]}" frip_consensus_worker "$i"
    done
    frip_status=0
    parallel_pool_wait_all || frip_status=$?
    (( frip_status == 0 )) || { warn "Consensus FRiP worker(s) failed: $(parallel_failed_labels_csv)"; exit 1; }
    for i in "${!SAMPLE_KEYS[@]}"; do
        cat "${FRIP_WORK_DIR}/${i}.row" >> "$FRIP_CONSENSUS_TSV"
        cat "${FRIP_WORK_DIR}/${i}.timing" >> "$PARALLEL_TIMING_TSV"
    done
    log "FRiP consensus: $FRIP_CONSENSUS_TSV"

    CONSENSUS_SIZEFACTORS_TSV="${QC_TABLES}/consensus_sizeFactors.tsv"
    CONSENSUS_COUNTS_TSV="${QC_MATRICES}/consensus_peak_counts.tsv"
    CONSENSUS_NORM_COUNTS_TSV="${QC_MATRICES}/consensus_peak_normCounts.tsv"

    R_SCRIPT="${R_BIN:-Rscript}"
    if [[ "$CURRENT_DESEQ2_TRACKS_REQUIRED" == "true" ]] && \
       command -v "$R_SCRIPT" >/dev/null 2>&1 && [[ -f "${QC_MATRICES}/multiBamSummary_peaks.tab" ]]; then
        "$R_SCRIPT" "${SCRIPT_DIR}/consensus_peak_size_factors.R" \
            "$SAMPLESHEET" "${QC_MATRICES}/multiBamSummary_peaks.tab" "${QC_TABLES}" \
            "$CONSENSUS_COUNTS_TSV" "$CONSENSUS_NORM_COUNTS_TSV" \
            >> "$MAIN_LOG" 2>&1 \
        && log "Consensus peak count matrix + DESeq2 size factors written" \
        || { warn "Consensus peak size factor estimation failed"; exit 1; }
    else
        if [[ "$CURRENT_DESEQ2_TRACKS_REQUIRED" == "true" ]]; then
            warn "Rscript or consensus peak-count matrix is unavailable"
            exit 1
        fi
    fi

    if [[ "$CURRENT_DESEQ2_TRACKS_REQUIRED" == "true" && -s "$CONSENSUS_SIZEFACTORS_TSV" ]] && \
       command -v bamCoverage >/dev/null 2>&1; then
        log "=== Phase 5: DESeq2 consensus and robust CPM bigWig/bedGraph generation ==="
        TRACK_NORMALIZATION_TSV="${QC_TABLES}/track_normalization_metadata.tsv"
        printf 'key\tpolicy\tsource_bam_stage\tduplicate_policy\tbowtie2_reporting\tminimum_mapq\tsecondary_alignments\tsupplementary_alignments\tgenome\tlayout\tsignal_unit\tsignal_count\tmapq_0\tmapq_1_9\tmapq_10_29\tmapq_ge_30\txs_tagged_signal_records\tconsensus_peak_sha256\tconsensus_count_sum\tcohort_geometric_mean_column_sum\tsize_factor\tdeseq2_consensus_scale\tdeseq2_robust_cpm_scale\n' \
            > "$TRACK_NORMALIZATION_TSV"
        TRACK_WORK_DIR="${QC_WORK_ROOT}/tracks"
        mkdir -p "$TRACK_WORK_DIR"

        track_worker_impl() {
            local i="$1" key bam size_factor robust_scale consensus_count_sum cohort_geomean
            key="${SAMPLE_KEYS[$i]}"; bam="${SAMPLE_BAMS[$i]}"
            size_factor=$(consensus_size_factor_for_key "$CONSENSUS_SIZEFACTORS_TSV" "$key" 2>/dev/null || true)
            [[ -n "$size_factor" ]] || { echo "No valid DESeq2 size factor found for $key" >&2; return 1; }
            robust_scale=$(consensus_table_value_for_key "$CONSENSUS_SIZEFACTORS_TSV" "$key" deseq2_robust_cpm_scale 2>/dev/null || true)
            consensus_count_sum=$(consensus_table_value_for_key "$CONSENSUS_SIZEFACTORS_TSV" "$key" consensus_count_sum 2>/dev/null || true)
            cohort_geomean=$(consensus_table_value_for_key "$CONSENSUS_SIZEFACTORS_TSV" "$key" cohort_geometric_mean_column_sum 2>/dev/null || true)
            [[ -n "$robust_scale" && -n "$consensus_count_sum" && -n "$cohort_geomean" ]] || {
                echo "Incomplete robust CPM metadata for $key" >&2; return 1;
            }
            generate_deseq2_tracks "$bam" "$key" "${SAMPLE_LAYOUT[$key]}" "$size_factor" \
                "$robust_scale" "$consensus_count_sum" "$cohort_geomean" "${TRACK_WORK_DIR}/${i}.metadata"
        }

        track_worker() {
            local i="$1" key start end status=0
            key="${SAMPLE_KEYS[$i]}"; start="$(date +%s)"
            track_worker_impl "$i" > "${QC_LOGS}/${key}.deseq2_tracks.log" 2>&1 || status=$?
            end="$(date +%s)"
            parallel_write_timing_row "${TRACK_WORK_DIR}/${i}.timing" deseq2_tracks "$key" \
                "$start" "$end" "$([[ $status -eq 0 ]] && echo SUCCESS || echo FAILED)" \
                "$TRACK_JOBS" "${THREADS_BIGWIG:-2}"
            return "$status"
        }

        parallel_pool_init "$TRACK_JOBS"
        for i in "${!SAMPLE_KEYS[@]}"; do
            log "QUEUE DESeq2 tracks: ${SAMPLE_KEYS[$i]}"
            parallel_pool_submit "${SAMPLE_KEYS[$i]}" track_worker "$i"
        done
        track_status=0
        parallel_pool_wait_all || track_status=$?
        (( track_status == 0 )) || { warn "DESeq2 track worker(s) failed: $(parallel_failed_labels_csv)"; exit 1; }
        for i in "${!SAMPLE_KEYS[@]}"; do
            cat "${TRACK_WORK_DIR}/${i}.metadata" >> "$TRACK_NORMALIZATION_TSV"
            cat "${TRACK_WORK_DIR}/${i}.timing" >> "$PARALLEL_TIMING_TSV"
        done
        log "DESeq2 consensus tracks: $BIGWIG_PEAKNORM_DIR"
        log "DESeq2 stringent robust CPM tracks: $BIGWIG_ROBUST_STRINGENT_DIR"
        log "Track normalization metadata: $TRACK_NORMALIZATION_TSV"
        log "Parallel job timing: $PARALLEL_TIMING_TSV"
    else
        if [[ "$CURRENT_DESEQ2_TRACKS_REQUIRED" == "true" ]]; then
            warn "Missing size factors or bamCoverage for required DESeq2 tracks"
            exit 1
        fi
    fi
fi

# =============================================================================
log ""
log "=== Post-alignment QC complete ==="
log "Summary      : $SUMMARY_TSV"
log "Warnings     : $WARNINGS_TSV"
log "Karyograms   : $QC_CHRPLOTS/"
log "Plots        : $QC_PLOTS/"
log "Matrices     : $QC_MATRICES/"
log ""
log "Output layout:"
log "  ${OUT_DIR}/"
log "    tables/    — qc_summary.tsv, qc_warnings.tsv, fingerprint_metrics.tsv"
log "               frip_consensus.tsv, per_chromosome_reads.tsv"
log "    plots/"
log "      chromosome_coverage/  — <KEY>_karyogram.png (per sample)"
log "                              karyogram_all_samples.png (multi-sample)"
log "                              <KEY>_100kb.bedGraph, <KEY>_per_chrom.tsv"
log "      fingerprint.png, correlation heatmaps, PCA, signal heatmap, profile"
log "    matrices/  — raw count matrices, PCA data"
log "    logs/      — main log"
log "    peak_sets/ — merged_narrow.bed, merged_broad.bed, consensus_peaks.bed"
log "    deeptools/ — .npz files (bins, peaks)"
