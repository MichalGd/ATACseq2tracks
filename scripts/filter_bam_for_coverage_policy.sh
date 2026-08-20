#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - policy-specific BAM filtering for coverage sensitivity tracks
# Usage: filter_bam_for_coverage_policy.sh <input.bam> <blacklist.bed> <output.bam> <PE|SE> <genome> <policy> <min_mapq>
set -euo pipefail

[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || {
    echo "ERROR: F2T_CONFIG is not set" >&2
    exit 1
}
# shellcheck disable=SC1090
source "$F2T_CONFIG"

INPUT_BAM="${1:?input BAM required}"
BLACKLIST_BED="${2:?blacklist BED required}"
OUTPUT_BAM="${3:?output BAM required}"
LAYOUT="${4:?PE or SE required}"
GENOME="${5:?genome required}"
POLICY="${6:?policy required}"
MINIMUM_MAPQ="${7:?minimum MAPQ required}"
THREADS="${THREADS_SAMTOOLS:-2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "${SCRIPT_DIR}/track_normalization_helpers.sh"

case "$POLICY" in
    permissive|intermediate) ;;
    *) echo "ERROR: unsupported coverage filtering policy: $POLICY" >&2; exit 1 ;;
esac
[[ "$LAYOUT" == "PE" || "$LAYOUT" == "SE" ]] || {
    echo "ERROR: layout must be PE or SE: $LAYOUT" >&2
    exit 1
}
[[ "$MINIMUM_MAPQ" =~ ^[0-9]+$ ]] || {
    echo "ERROR: minimum MAPQ must be a non-negative integer: $MINIMUM_MAPQ" >&2
    exit 1
}
[[ -s "$INPUT_BAM" ]] || { echo "ERROR: input BAM missing or empty: $INPUT_BAM" >&2; exit 1; }
[[ -s "$BLACKLIST_BED" ]] || { echo "ERROR: blacklist missing or empty: $BLACKLIST_BED" >&2; exit 1; }

mkdir -p "$(dirname "$OUTPUT_BAM")"
samtools quickcheck "$INPUT_BAM"
[[ -s "${INPUT_BAM}.bai" ]] || samtools index -@ "$THREADS" "$INPUT_BAM"

tmp_dir="$(mktemp -d "$(dirname "$OUTPUT_BAM")/.coverage-policy.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT
prefilter="${tmp_dir}/prefilter.bam"
canonical="${tmp_dir}/canonical.bam"
no_blacklist="${tmp_dir}/no_blacklist.bam"

input_records="$(samtools view -@ "$THREADS" -c "$INPUT_BAM")"
if [[ "$LAYOUT" == "PE" ]]; then
    # Require a proper pair; exclude unmapped, mate-unmapped, secondary,
    # QC-failed and supplementary records. Duplicate records are deliberately
    # not excluded here: the source BAM defines whether duplicates were removed.
    samtools view -@ "$THREADS" -b -q "$MINIMUM_MAPQ" -f 2 -F 2828 \
        -o "$prefilter" "$INPUT_BAM"
else
    samtools view -@ "$THREADS" -b -q "$MINIMUM_MAPQ" -F 2820 \
        -o "$prefilter" "$INPUT_BAM"
fi
after_flags="$(samtools view -@ "$THREADS" -c "$prefilter")"

samtools index -@ "$THREADS" "$prefilter"
mapfile -t contigs < <(canonical_contigs_from_bam "$prefilter" "$GENOME")
(( ${#contigs[@]} > 0 )) || {
    echo "ERROR: no canonical chromosomes remain for $OUTPUT_BAM" >&2
    exit 1
}
samtools view -@ "$THREADS" -b -o "$canonical" "$prefilter" "${contigs[@]}"
after_canonical="$(samtools view -@ "$THREADS" -c "$canonical")"

bedtools intersect -v -abam "$canonical" -b "$BLACKLIST_BED" > "$no_blacklist"
if [[ "$LAYOUT" == "PE" ]]; then
    samtools sort -n -@ "$THREADS" -o "${tmp_dir}/name.bam" "$no_blacklist"
    samtools fixmate -r -@ "$THREADS" "${tmp_dir}/name.bam" "${tmp_dir}/fixmate.bam"
    samtools view -@ "$THREADS" -b -f 2 -o "${tmp_dir}/paired.bam" "${tmp_dir}/fixmate.bam"
    samtools sort -@ "$THREADS" -o "${tmp_dir}/final.bam" "${tmp_dir}/paired.bam"
else
    samtools sort -@ "$THREADS" -o "${tmp_dir}/final.bam" "$no_blacklist"
fi

final_records="$(samtools view -@ "$THREADS" -c "${tmp_dir}/final.bam")"
if (( final_records == 0 )) && [[ "${ALLOW_EMPTY_FILTERED_BAM:-false}" != "true" ]]; then
    echo "ERROR: $POLICY filtering removed every alignment from $INPUT_BAM" >&2
    exit 1
fi
samtools quickcheck "${tmp_dir}/final.bam"
samtools index -@ "$THREADS" "${tmp_dir}/final.bam"
mv "${tmp_dir}/final.bam" "$OUTPUT_BAM"
mv "${tmp_dir}/final.bam.bai" "${OUTPUT_BAM}.bai"

signal_count="$(signal_count_for_bam "$OUTPUT_BAM" "$LAYOUT")"
metrics="${OUTPUT_BAM%.bam}.filtering_metadata.tsv"
printf 'policy\tsource_bam\toutput_bam\tgenome\tlayout\tminimum_mapq\tduplicates_removed_upstream\tinclude_flags\texclude_flags\tinput_alignment_records\tafter_primary_pair_mapq\tafter_canonical_chromosomes\tafter_blacklist_pair_cleanup\tsignal_count\n' > "$metrics"
if [[ "$POLICY" == "intermediate" ]]; then duplicates_removed=true; else duplicates_removed=false; fi
if [[ "$LAYOUT" == "PE" ]]; then include_flags=2; exclude_flags=2828; else include_flags=0; exclude_flags=2820; fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$POLICY" "$INPUT_BAM" "$OUTPUT_BAM" "$GENOME" "$LAYOUT" "$MINIMUM_MAPQ" \
    "$duplicates_removed" "$include_flags" "$exclude_flags" "$input_records" \
    "$after_flags" "$after_canonical" "$final_records" "$signal_count" >> "$metrics"
echo "Coverage policy BAM: policy=$POLICY signal_count=$signal_count output=$OUTPUT_BAM"
