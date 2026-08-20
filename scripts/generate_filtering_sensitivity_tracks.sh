#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - permissive/intermediate robust-CPM sensitivity tracks
# Usage: generate_filtering_sensitivity_tracks.sh <samplesheet.csv> <output_dir> <consensus_peaks.bed>
set -euo pipefail

[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || {
    echo "ERROR: F2T_CONFIG is not set" >&2
    exit 1
}
# shellcheck disable=SC1090
source "$F2T_CONFIG"

SAMPLESHEET="${1:?samplesheet required}"
OUTPUT_DIR="${2:?workflow output directory required}"
CONSENSUS_PEAK="${3:?consensus peak BED required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRACK_HELPERS="${SCRIPT_DIR}/track_normalization_helpers.sh"
QC_HELPERS="${SCRIPT_DIR}/qc_table_helpers.sh"
PARALLEL_HELPERS="${SCRIPT_DIR}/parallel_job_helpers.sh"
# shellcheck disable=SC1090
source "$TRACK_HELPERS"
# shellcheck disable=SC1090
source "$QC_HELPERS"
# shellcheck disable=SC1090
source "$PARALLEL_HELPERS"

[[ -s "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet missing or empty: $SAMPLESHEET" >&2; exit 1; }
[[ -s "$CONSENSUS_PEAK" ]] || { echo "ERROR: consensus peak BED missing or empty: $CONSENSUS_PEAK" >&2; exit 1; }

POLICY_ROOT="${OUTPUT_DIR}/coverage_filtering_sensitivity"
POLICY_BAM_ROOT="${OUTPUT_DIR}/coverage_filtering_policy_bams"
ROBUST_TRACK_ROOT="${OUTPUT_DIR}/bigwig_deseq2_robust_cpm"
LOG_DIR="${POLICY_ROOT}/logs"
mkdir -p "$POLICY_ROOT" "$POLICY_BAM_ROOT" "$ROBUST_TRACK_ROOT" "$LOG_DIR"
MAIN_LOG="${LOG_DIR}/filtering_sensitivity.log"
TRACK_JOBS="${TRACK_PARALLEL_JOBS:-2}"
FILTER_JOBS="${COVERAGE_FILTER_PARALLEL_JOBS:-${TRACK_PARALLEL_JOBS:-2}}"
THREADS="${THREADS_DEEPTOOLS:-8}"
parallel_require_positive_integer TRACK_PARALLEL_JOBS "$TRACK_JOBS"
parallel_require_positive_integer COVERAGE_FILTER_PARALLEL_JOBS "$FILTER_JOBS"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$1" | tee -a "$MAIN_LOG"; }

run_policy() {
    local policy="$1" title source_dir source_suffix minimum_mapq duplicate_policy
    local policy_dir matrix_dir table_dir track_dir bam_dir count_tab raw_counts norm_counts sf_table
    local normalization_layout genome consensus_sha
    local -a signal_args=()

    case "$policy" in
        permissive)
            title="Permissive"
            source_dir="${OUTPUT_DIR}/bams"
            source_suffix=""
            minimum_mapq="${PERMISSIVE_MIN_MAPQ:-0}"
            duplicate_policy="retained"
            ;;
        intermediate)
            title="Intermediate"
            source_dir="${OUTPUT_DIR}/dedupBams"
            source_suffix="_dedup"
            minimum_mapq="${INTERMEDIATE_MIN_MAPQ:-0}"
            duplicate_policy="removed_by_Picard"
            ;;
        *) echo "ERROR: unsupported policy: $policy" >&2; return 1 ;;
    esac

    policy_dir="${POLICY_ROOT}/${policy}"
    matrix_dir="${policy_dir}/matrices"
    table_dir="${policy_dir}/tables"
    track_dir="${ROBUST_TRACK_ROOT}/${policy}"
    bam_dir="${POLICY_BAM_ROOT}/${policy}"
    mkdir -p "$matrix_dir" "$table_dir" "$track_dir" "$bam_dir" "${LOG_DIR}/${policy}"

    POLICY_KEYS=()
    POLICY_BAMS=()
    POLICY_LAYOUTS=()
    POLICY_GENOMES=()
    POLICY_BLACKLISTS=()
    declare -A seen=()
    while IFS=',' read -r sid fq1 fq2 layout sample_genome assay factor condition treatment \
        cell_type rep tech_rep is_control control_id macs2_mode blacklist rest; do
        [[ "$sid" == "sample_id" ]] && continue
        for name in sid layout sample_genome rep is_control blacklist; do
            value="${!name}"; printf -v "$name" '%s' "${value//\"/}"
        done
        [[ "${is_control,,}" =~ ^(true|1|yes)$ ]] && continue
        key="${sid}_bioR${rep}"
        [[ -n "${seen[$key]+x}" ]] && continue
        seen["$key"]=1
        POLICY_KEYS+=("$key")
        POLICY_LAYOUTS+=("${layout^^}")
        POLICY_GENOMES+=("$sample_genome")
        POLICY_BLACKLISTS+=("$blacklist")
        POLICY_BAMS+=("${bam_dir}/${key}_${policy}.bam")
    done < "$SAMPLESHEET"
    (( ${#POLICY_KEYS[@]} > 0 )) || { echo "ERROR: no biological samples for $policy" >&2; return 1; }

    normalization_layout="${POLICY_LAYOUTS[0]}"
    genome="${POLICY_GENOMES[0]}"
    for i in "${!POLICY_KEYS[@]}"; do
        [[ "${POLICY_LAYOUTS[$i]}" == "$normalization_layout" ]] || {
            echo "ERROR: mixed PE/SE layouts cannot share a normalization cohort" >&2; return 1;
        }
        [[ "${POLICY_GENOMES[$i]}" == "$genome" ]] || {
            echo "ERROR: mixed genome builds cannot share a normalization cohort" >&2; return 1;
        }
    done
    deeptools_signal_args "$normalization_layout" signal_args

    filter_worker() {
        local i="$1" key input output
        key="${POLICY_KEYS[$i]}"
        input="${source_dir}/${key}${source_suffix}.bam"
        output="${POLICY_BAMS[$i]}"
        [[ -s "$input" ]] || { echo "ERROR: source BAM missing: $input" >&2; return 1; }
        [[ -s "${POLICY_BLACKLISTS[$i]}" ]] || {
            echo "ERROR: blacklist missing for $key: ${POLICY_BLACKLISTS[$i]}" >&2; return 1;
        }
        if [[ -s "$output" ]] && samtools quickcheck "$output" 2>/dev/null; then
            echo "SKIP existing policy BAM: $output"
            return 0
        fi
        bash "${SCRIPT_DIR}/filter_bam_for_coverage_policy.sh" \
            "$input" "${POLICY_BLACKLISTS[$i]}" "$output" "${POLICY_LAYOUTS[$i]}" \
            "${POLICY_GENOMES[$i]}" "$policy" "$minimum_mapq"
    }

    log "Filtering $policy policy BAMs with MAPQ >=${minimum_mapq}"
    parallel_pool_init "$FILTER_JOBS"
    for i in "${!POLICY_KEYS[@]}"; do
        parallel_pool_submit "${POLICY_KEYS[$i]}" filter_worker "$i" \
            > "${LOG_DIR}/${policy}/${POLICY_KEYS[$i]}.filter.log" 2>&1
    done
    parallel_pool_wait_all || {
        echo "ERROR: $policy BAM filtering failed: $(parallel_failed_labels_csv)" >&2
        return 1
    }

    count_tab="${matrix_dir}/multiBamSummary_consensus_peaks.tab"
    multiBamSummary BED-file --BED "$CONSENSUS_PEAK" \
        -b "${POLICY_BAMS[@]}" --labels "${POLICY_KEYS[@]}" "${signal_args[@]}" \
        -p "$THREADS" -o "${matrix_dir}/multiBamSummary_consensus_peaks.npz" \
        --outRawCounts "$count_tab" >> "$MAIN_LOG" 2>&1

    raw_counts="${matrix_dir}/consensus_peak_counts.tsv"
    norm_counts="${matrix_dir}/consensus_peak_normCounts.tsv"
    "${R_BIN:-Rscript}" "${SCRIPT_DIR}/consensus_peak_size_factors.R" \
        "$SAMPLESHEET" "$count_tab" "$table_dir" "$raw_counts" "$norm_counts" \
        >> "$MAIN_LOG" 2>&1
    sf_table="${table_dir}/consensus_sizeFactors.tsv"
    [[ -s "$sf_table" ]] || { echo "ERROR: size-factor table missing: $sf_table" >&2; return 1; }

    consensus_sha="$(sha256sum "$CONSENSUS_PEAK" | awk '{print $1}')"
    metadata="${table_dir}/track_normalization_metadata.tsv"
    printf 'key\tpolicy\tsource_bam_stage\tduplicate_policy\tbowtie2_reporting\tminimum_mapq\tsecondary_alignments\tsupplementary_alignments\tgenome\tlayout\tsignal_unit\tsignal_count\tmapq_0\tmapq_1_9\tmapq_10_29\tmapq_ge_30\txs_tagged_signal_records\tconsensus_peak_sha256\tconsensus_count_sum\tcohort_geometric_mean_column_sum\tsize_factor\tdeseq2_consensus_scale\tdeseq2_robust_cpm_scale\n' > "$metadata"

    track_worker() {
        local i="$1" key bam layout signal_unit signal_count sf consensus_scale robust_scale count_sum geomean bins
        local bw bg tmp_dir
        key="${POLICY_KEYS[$i]}"; bam="${POLICY_BAMS[$i]}"; layout="${POLICY_LAYOUTS[$i]}"
        signal_unit="$(signal_unit_for_layout "$layout")"
        signal_count="$(signal_count_for_bam "$bam" "$layout")"
        sf="$(consensus_size_factor_for_key "$sf_table" "$key")"
        consensus_scale="$(awk -v value="$sf" 'BEGIN{printf "%.12g", 1/value}')"
        robust_scale="$(consensus_table_value_for_key "$sf_table" "$key" deseq2_robust_cpm_scale)"
        count_sum="$(consensus_table_value_for_key "$sf_table" "$key" consensus_count_sum)"
        geomean="$(consensus_table_value_for_key "$sf_table" "$key" cohort_geometric_mean_column_sum)"
        tmp_dir="$(mktemp -d "${track_dir}/.${key}.XXXXXX")"
        if [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" ]]; then
            bw="${track_dir}/${key}_DESeq2RobustCPM_${title}.bw"
            write_scaled_coverage_track "$bam" "${tmp_dir}/$(basename "$bw")" bigwig \
                "$robust_scale" "$layout" "${TRACK_BIN_SIZE:-10}" "${THREADS_BIGWIG:-2}"
            mv "${tmp_dir}/$(basename "$bw")" "$bw"
        fi
        if [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" ]]; then
            bg="${track_dir}/${key}_DESeq2RobustCPM_${title}.bedGraph"
            write_scaled_coverage_track "$bam" "${tmp_dir}/$(basename "$bg")" bedgraph \
                "$robust_scale" "$layout" "${TRACK_BIN_SIZE:-10}" "${THREADS_BIGWIG:-2}"
            mv "${tmp_dir}/$(basename "$bg")" "$bg"
        fi
        rmdir "$tmp_dir"
        bins="$(signal_mapq_bin_counts "$bam" "$layout")"
        printf '%s\t%s\t%s\t%s\t%s\t%s\texcluded\texcluded\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$key" "$policy" "$source_dir" "$duplicate_policy" "one_primary_default" \
            "$minimum_mapq" "${POLICY_GENOMES[$i]}" "$layout" "$signal_unit" "$signal_count" \
            $bins "$consensus_sha" "$count_sum" "$geomean" "$sf" "$consensus_scale" "$robust_scale" \
            > "${policy_dir}/.${i}.metadata"
    }

    log "Generating $policy robust-CPM tracks"
    parallel_pool_init "$TRACK_JOBS"
    for i in "${!POLICY_KEYS[@]}"; do
        parallel_pool_submit "${POLICY_KEYS[$i]}" track_worker "$i" \
            > "${LOG_DIR}/${policy}/${POLICY_KEYS[$i]}.tracks.log" 2>&1
    done
    parallel_pool_wait_all || {
        echo "ERROR: $policy track generation failed: $(parallel_failed_labels_csv)" >&2
        return 1
    }
    for i in "${!POLICY_KEYS[@]}"; do
        cat "${policy_dir}/.${i}.metadata" >> "$metadata"
        rm -f "${policy_dir}/.${i}.metadata"
    done
    log "Completed $policy policy: tracks=$track_dir metadata=$metadata"
}

if [[ "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" ]]; then
    run_policy permissive
else
    log "Permissive robust-CPM tracks disabled"
fi

if [[ "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" ]]; then
    run_policy intermediate
else
    log "Intermediate robust-CPM tracks disabled"
fi

COMBINED_METADATA="${POLICY_ROOT}/track_normalization_metadata.tsv"
: > "$COMBINED_METADATA"
METADATA_TABLES=()
[[ "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" ]] && \
    METADATA_TABLES+=("${POLICY_ROOT}/permissive/tables/track_normalization_metadata.tsv")
[[ "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" ]] && \
    METADATA_TABLES+=("${POLICY_ROOT}/intermediate/tables/track_normalization_metadata.tsv")
[[ "${GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true}" == "true" ]] && \
    METADATA_TABLES+=("${OUTPUT_DIR}/qc_post_alignment/tables/track_normalization_metadata.tsv")
for table in "${METADATA_TABLES[@]}"; do
    [[ -s "$table" ]] || continue
    if [[ ! -s "$COMBINED_METADATA" ]]; then
        cat "$table" > "$COMBINED_METADATA"
    else
        tail -n +2 "$table" >> "$COMBINED_METADATA"
    fi
done
[[ -s "$COMBINED_METADATA" ]] || {
    echo "ERROR: no filtering-policy normalization metadata were generated" >&2
    exit 1
}
log "Combined policy metadata: $COMBINED_METADATA"
log "Coverage filtering sensitivity module complete"
