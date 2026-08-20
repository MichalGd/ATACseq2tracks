#!/usr/bin/env bash
# ATACseq2tracks v4.2.0 - stringent dm6 spike-in calibrated host tracks
# Usage: drosophila_spikein_tracks.sh <samplesheet.csv> <trimmed_fastq_dir> <output_dir>
set -euo pipefail

[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"

SAMPLESHEET="${1:?samplesheet required}"
TRIM_DIR="${2:?trimmed FASTQ directory required}"
OUTPUT_DIR="${3:?workflow output directory required}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "${SCRIPT_DIR}/track_normalization_helpers.sh"
# shellcheck disable=SC1090
source "${SCRIPT_DIR}/parallel_job_helpers.sh"
# shellcheck disable=SC1090
source "${SCRIPT_DIR}/qc_table_helpers.sh"

[[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]] || {
    echo "Drosophila spike-in tracks disabled"
    exit 0
}

SPIKE_ROOT="${OUTPUT_DIR}/spikein"
COMPOSITE_BAM_DIR="${SPIKE_ROOT}/composite_bams"
DEDUP_BAM_DIR="${SPIKE_ROOT}/dedup_bams"
HOST_BAM_DIR="${SPIKE_ROOT}/stringent_host_bams"
DM6_BAM_DIR="${SPIKE_ROOT}/stringent_dm6_bams"
LOG_DIR="${SPIKE_ROOT}/logs"
TABLE_DIR="${SPIKE_ROOT}/tables"
PER_SAMPLE_TABLE_DIR="${TABLE_DIR}/per_sample"
TRACK_DIR="${OUTPUT_DIR}/bigwig_spikein/stringent"
mkdir -p "$COMPOSITE_BAM_DIR" "$DEDUP_BAM_DIR" "$HOST_BAM_DIR" "$DM6_BAM_DIR" \
    "$LOG_DIR" "$TABLE_DIR" "$PER_SAMPLE_TABLE_DIR" "$TRACK_DIR"

SAMPLE_TABLE="${TABLE_DIR}/declared_spikein_samples.tsv"
python3 "${SCRIPT_DIR}/extract_spikein_samples.py" "$SAMPLESHEET" > "$SAMPLE_TABLE"
[[ -s "$SAMPLE_TABLE" ]] || { echo "ERROR: no declared dm6 spike-in samples" >&2; exit 1; }

THREADS_ALIGN_VALUE="${THREADS_ALIGN:-2}"
THREADS_SAMTOOLS_VALUE="${THREADS_SAMTOOLS:-2}"
THREADS_BIGWIG_VALUE="${THREADS_BIGWIG:-2}"
PARALLEL_JOBS="${SPIKEIN_PARALLEL_JOBS:-2}"
MINIMUM_MAPQ="${SPIKEIN_MIN_MAPQ:-30}"
SCALE_TARGET="${SPIKEIN_SCALE_TARGET:-1000000}"
MIN_FAIL="${SPIKEIN_MIN_FRAGMENTS_FAIL:-1000}"
MIN_WARN="${SPIKEIN_MIN_FRAGMENTS_WARN:-10000}"
LOW_FRACTION="${SPIKEIN_WARN_LOW_FRACTION:-0.001}"
HIGH_FRACTION="${SPIKEIN_WARN_HIGH_FRACTION:-0.20}"

parallel_require_positive_integer SPIKEIN_PARALLEL_JOBS "$PARALLEL_JOBS"
for numeric in "$MINIMUM_MAPQ" "$MIN_FAIL" "$MIN_WARN"; do
    [[ "$numeric" =~ ^[0-9]+$ ]] || { echo "ERROR: spike-in count/MAPQ settings must be integers" >&2; exit 1; }
done
awk -v value="$SCALE_TARGET" 'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value > 0)}' \
    || { echo "ERROR: SPIKEIN_SCALE_TARGET must be positive" >&2; exit 1; }

KEYS=(); SAMPLE_IDS=(); REPLICATES=(); LAYOUTS=(); HOST_GENOMES=(); HOST_BLACKLISTS=()
RATIOS=(); STAGES=(); CONTROL_STATES=(); COMPOSITE_INDICES=()
while IFS=$'\t' read -r key sample_id replicate layout host_genome host_blacklist spike_genome ratio stage is_control; do
    [[ "$key" == "key" ]] && continue
    case "$host_genome" in
        hg38) composite_index="${INDEX_HG38_DM6:-}" ;;
        mm39) composite_index="${INDEX_MM39_DM6:-}" ;;
        *) echo "ERROR: unsupported host genome for spike-in: $host_genome" >&2; exit 1 ;;
    esac
    [[ -n "$composite_index" ]] || { echo "ERROR: composite Bowtie2 index is unset for $host_genome" >&2; exit 1; }
    [[ -s "$host_blacklist" ]] || { echo "ERROR: host blacklist missing for $key: $host_blacklist" >&2; exit 1; }
    KEYS+=("$key"); SAMPLE_IDS+=("$sample_id"); REPLICATES+=("$replicate"); LAYOUTS+=("$layout")
    HOST_GENOMES+=("$host_genome"); HOST_BLACKLISTS+=("$host_blacklist"); RATIOS+=("$ratio")
    STAGES+=("$stage"); CONTROL_STATES+=("$is_control"); COMPOSITE_INDICES+=("$composite_index")
done < "$SAMPLE_TABLE"
(( ${#KEYS[@]} > 0 )) || { echo "ERROR: no valid dm6 spike-in samples" >&2; exit 1; }
[[ -s "${BLACKLIST_DM6:-}" ]] || { echo "ERROR: BLACKLIST_DM6 is missing: ${BLACKLIST_DM6:-unset}" >&2; exit 1; }
[[ -s "${CHROM_SIZES_DM6:-}" ]] || { echo "ERROR: CHROM_SIZES_DM6 is missing: ${CHROM_SIZES_DM6:-unset}" >&2; exit 1; }

check_index() {
    local prefix="$1" extension component missing
    for extension in bt2 bt2l; do
        missing=0
        for component in 1 2 3 4; do [[ -s "${prefix}.${component}.${extension}" ]] || missing=$((missing + 1)); done
        for component in 1 2; do [[ -s "${prefix}.rev.${component}.${extension}" ]] || missing=$((missing + 1)); done
        (( missing == 0 )) && return 0
    done
    return 1
}
for index in "${COMPOSITE_INDICES[@]}"; do
    check_index "$index" || { echo "ERROR: incomplete composite Bowtie2 index: $index" >&2; exit 1; }
done

WARNINGS_FILE="${TABLE_DIR}/spikein_warnings.tsv"
printf 'key\twarning_type\tvalue\tthreshold\tmessage\n' > "$WARNINGS_FILE"

process_sample() {
    local i="$1" key layout host_genome host_blacklist ratio index composite_bam dedup_bam
    local host_bam dm6_bam metrics_file row_file warning_file log_file r1 r2 staging
    local host_count dm6_count combined_count fly_fraction scale duplicate_pct q0 q30 low_mapq
    local bw bg tmp_track_dir index_sha host_blacklist_sha dm6_blacklist_sha outputs_complete
    key="${KEYS[$i]}"; layout="${LAYOUTS[$i]}"; host_genome="${HOST_GENOMES[$i]}"
    host_blacklist="${HOST_BLACKLISTS[$i]}"; ratio="${RATIOS[$i]}"; index="${COMPOSITE_INDICES[$i]}"
    composite_bam="${COMPOSITE_BAM_DIR}/${key}_host_dm6.bam"
    dedup_bam="${DEDUP_BAM_DIR}/${key}_host_dm6_dedup.bam"
    host_bam="${HOST_BAM_DIR}/${key}_host_stringent.bam"
    dm6_bam="${DM6_BAM_DIR}/${key}_dm6_stringent.bam"
    metrics_file="${LOG_DIR}/${key}.picard_dup_metrics.txt"
    row_file="${PER_SAMPLE_TABLE_DIR}/${key}.normalization.tsv"
    warning_file="${PER_SAMPLE_TABLE_DIR}/${key}.warnings.tsv"
    log_file="${LOG_DIR}/${key}.log"
    bw="${TRACK_DIR}/${key}_SpikeInDM6_Stringent.bw"
    bg="${TRACK_DIR}/${key}_SpikeInDM6_Stringent.bedGraph"
    : > "$warning_file"

    outputs_complete=true
    [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" && ! -s "$bw" ]] && outputs_complete=false
    [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" && ! -s "$bg" ]] && outputs_complete=false
    [[ ! -s "$row_file" ]] && outputs_complete=false
    if [[ "$outputs_complete" == "true" ]]; then
        echo "SKIP existing spike-in tracks: $key"
        return 0
    fi

    if [[ ! -s "$composite_bam" ]] || ! samtools quickcheck "$composite_bam" 2>/dev/null; then
        staging="${COMPOSITE_BAM_DIR}/.${key}.staging.$$.bam"
        # Keep competitive assignment reproducible and emit only Bowtie2's one
        # best primary placement. In particular, do not inherit a host branch
        # that might contain -k/-a and report multiple placements.
        extra_args=(--very-sensitive)
        common=( -x "$index" -p "$THREADS_ALIGN_VALUE" --rg-id "$key" --rg "SM:${key}" --rg 'PL:ILLUMINA' )
        if [[ "$layout" == "PE" ]]; then
            r1="${TRIM_DIR}/${key}_1_val_1.fq.gz"; r2="${TRIM_DIR}/${key}_2_val_2.fq.gz"
            [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: missing trimmed PE FASTQs for $key" >&2; return 1; }
            bowtie2 "${common[@]}" "${extra_args[@]}" --no-mixed --no-discordant --dovetail \
                -1 "$r1" -2 "$r2" 2> "$log_file" \
                | samtools sort -@ "$THREADS_SAMTOOLS_VALUE" -o "$staging"
        else
            r1="${TRIM_DIR}/${key}_trimmed.fq.gz"
            [[ -s "$r1" ]] || { echo "ERROR: missing trimmed SE FASTQ for $key" >&2; return 1; }
            bowtie2 "${common[@]}" "${extra_args[@]}" -U "$r1" 2> "$log_file" \
                | samtools sort -@ "$THREADS_SAMTOOLS_VALUE" -o "$staging"
        fi
        samtools quickcheck "$staging"; samtools index -@ "$THREADS_SAMTOOLS_VALUE" "$staging"
        mv "$staging" "$composite_bam"; mv "${staging}.bai" "${composite_bam}.bai"
    fi

    if [[ ! -s "$dedup_bam" ]] || ! samtools quickcheck "$dedup_bam" 2>/dev/null; then
        staging="${DEDUP_BAM_DIR}/.${key}.dedup.$$.bam"
        java "${PICARD_XMX:--Xmx8g}" -jar "$PICARD_JAR" MarkDuplicates \
            INPUT="$composite_bam" OUTPUT="$staging" METRICS_FILE="$metrics_file" \
            REMOVE_DUPLICATES=true OPTICAL_DUPLICATE_PIXEL_DISTANCE="${PICARD_OPTICAL_DISTANCE:-100}" \
            TMP_DIR="${PICARD_TMP:-/tmp}" ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT \
            >> "$log_file" 2>&1
        samtools quickcheck "$staging"; samtools index -@ "$THREADS_SAMTOOLS_VALUE" "$staging"
        mv "$staging" "$dedup_bam"; mv "${staging}.bai" "${dedup_bam}.bai"
    fi

    if [[ ! -s "$host_bam" || ! -s "$dm6_bam" ]] || \
       ! samtools quickcheck "$host_bam" "$dm6_bam" 2>/dev/null; then
        bash "${SCRIPT_DIR}/filter_composite_spikein_bam.sh" "$dedup_bam" "$host_genome" \
            "$host_blacklist" "$BLACKLIST_DM6" "$host_bam" "$dm6_bam" "$layout" \
            >> "$log_file" 2>&1
    fi

    host_count="$(signal_count_for_bam "$host_bam" "$layout")"
    dm6_count="$(signal_count_for_bam "$dm6_bam" "$layout")"
    (( dm6_count >= MIN_FAIL )) || {
        echo "ERROR: $key retained only $dm6_count dm6 observations; minimum is $MIN_FAIL" >&2
        return 1
    }
    combined_count=$((host_count + dm6_count))
    fly_fraction="$(awk -v fly="$dm6_count" -v total="$combined_count" 'BEGIN{printf "%.12g", fly/total}')"
    scale="$(spikein_scale_factor "$SCALE_TARGET" "$ratio" "$dm6_count")" || {
        echo "ERROR: invalid spike-in scale inputs for $key" >&2
        return 1
    }

    if [[ "$layout" == "PE" ]]; then
        q0="$(samtools view -c -f 66 -F 3852 "$dedup_bam")"
        q30="$(samtools view -c -q "$MINIMUM_MAPQ" -f 66 -F 3852 "$dedup_bam")"
    else
        q0="$(samtools view -c -F 3844 "$dedup_bam")"
        q30="$(samtools view -c -q "$MINIMUM_MAPQ" -F 3844 "$dedup_bam")"
    fi
    low_mapq=$((q0 - q30))
    duplicate_pct="$(picard_duplication_pct "$metrics_file" 2>/dev/null || echo NA)"

    if (( dm6_count < MIN_WARN )); then
        printf '%s\tlow_dm6_count\t%s\t%s\tRetained dm6 count has elevated sampling uncertainty\n' \
            "$key" "$dm6_count" "$MIN_WARN" >> "$warning_file"
    fi
    if awk -v value="$fly_fraction" -v threshold="$LOW_FRACTION" 'BEGIN{exit !(value < threshold)}'; then
        printf '%s\tlow_dm6_fraction\t%s\t%s\tDrosophila fraction is below the configured warning threshold\n' \
            "$key" "$fly_fraction" "$LOW_FRACTION" >> "$warning_file"
    fi
    if awk -v value="$fly_fraction" -v threshold="$HIGH_FRACTION" 'BEGIN{exit !(value > threshold)}'; then
        printf '%s\thigh_dm6_fraction\t%s\t%s\tDrosophila fraction may consume excessive sequencing capacity\n' \
            "$key" "$fly_fraction" "$HIGH_FRACTION" >> "$warning_file"
    fi

    tmp_track_dir="$(mktemp -d "${TRACK_DIR}/.${key}.XXXXXX")"
    trap 'rm -rf -- "$tmp_track_dir"' EXIT
    if [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" == "true" ]]; then
        write_scaled_coverage_track "$host_bam" "${tmp_track_dir}/$(basename "$bw")" bigwig \
            "$scale" "$layout" "${TRACK_BIN_SIZE:-10}" "$THREADS_BIGWIG_VALUE"
        mv "${tmp_track_dir}/$(basename "$bw")" "$bw"
    fi
    if [[ "${GENERATE_COVERAGE_BEDGRAPHS:-true}" == "true" ]]; then
        write_scaled_coverage_track "$host_bam" "${tmp_track_dir}/$(basename "$bg")" bedgraph \
            "$scale" "$layout" "${TRACK_BIN_SIZE:-10}" "$THREADS_BIGWIG_VALUE"
        mv "${tmp_track_dir}/$(basename "$bg")" "$bg"
    fi
    rmdir "$tmp_track_dir"; trap - EXIT

    index_sha="$({
        for suffix in 1 2 3 4; do
            sha256sum "${index}.${suffix}.bt2" 2>/dev/null || sha256sum "${index}.${suffix}.bt2l"
        done
        for suffix in 1 2; do
            sha256sum "${index}.rev.${suffix}.bt2" 2>/dev/null || sha256sum "${index}.rev.${suffix}.bt2l"
        done
    } | sha256sum | awk '{print $1}')"
    host_blacklist_sha="$(sha256sum "$host_blacklist" | awk '{print $1}')"
    dm6_blacklist_sha="$(sha256sum "$BLACKLIST_DM6" | awk '{print $1}')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$host_genome" dm6 "$layout" "$(signal_unit_for_layout "$layout")" \
        "${CONTROL_STATES[$i]}" "${STAGES[$i]}" "$ratio" "$MINIMUM_MAPQ" true \
        "$host_count" "$dm6_count" "$fly_fraction" "$q0" "$q30" "$low_mapq" \
        "$duplicate_pct" "$SCALE_TARGET" "$scale" "$index" "$index_sha" \
        "${host_blacklist_sha}:${dm6_blacklist_sha}" > "$row_file"
    echo "Spike-in tracks complete: key=$key host=$host_count dm6=$dm6_count scale=$scale"
}

parallel_pool_init "$PARALLEL_JOBS"
for i in "${!KEYS[@]}"; do
    parallel_pool_submit "${KEYS[$i]}" process_sample "$i" \
        > "${LOG_DIR}/${KEYS[$i]}.worker.log" 2>&1
done
parallel_pool_wait_all || {
    echo "ERROR: Drosophila spike-in processing failed: $(parallel_failed_labels_csv)" >&2
    exit 1
}

NORMALIZATION_TABLE="${TABLE_DIR}/spikein_normalization.tsv"
printf 'key\thost_genome\tspikein_genome\tlayout\tsignal_unit\tis_control\tspikein_stage\tspikein_to_host_ratio\tminimum_mapq\tduplicates_removed\thost_signal_count\tdm6_signal_count\tdm6_fraction\tcomposite_primary_q0\tcomposite_primary_q_threshold\tlow_mapq_primary\tpercent_duplication\tscale_target\tapplied_scale\tcomposite_index\tcomposite_index_sha256\tblacklist_sha256_host_dm6\n' \
    > "$NORMALIZATION_TABLE"
for key in "${KEYS[@]}"; do
    cat "${PER_SAMPLE_TABLE_DIR}/${key}.normalization.tsv" >> "$NORMALIZATION_TABLE"
    [[ -s "${PER_SAMPLE_TABLE_DIR}/${key}.warnings.tsv" ]] && \
        cat "${PER_SAMPLE_TABLE_DIR}/${key}.warnings.tsv" >> "$WARNINGS_FILE"
done

PROVENANCE="${TABLE_DIR}/spikein_provenance.tsv"
printf 'field\tvalue\nworkflow_version\t4.2.0\nformula\traw_host_coverage_x_scale_target_x_spikein_to_host_ratio_div_dm6_signal_count\ndm6_reference\tRelease_6_dm6\ndm6_chrom_sizes\t%s\ndm6_blacklist\t%s\nbowtie2_version\t%s\nsamtools_version\t%s\n' \
    "$CHROM_SIZES_DM6" "$BLACKLIST_DM6" \
    "$(bowtie2 --version | head -1 | tr '\t' ' ')" "$(samtools --version | head -1 | tr '\t' ' ')" \
    > "$PROVENANCE"

echo "Drosophila spike-in module complete"
echo "  Tracks    : $TRACK_DIR"
echo "  Metadata  : $NORMALIZATION_TABLE"
echo "  Warnings  : $WARNINGS_FILE"
