#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - canonical-chromosome fragment/read CPM tracks
# Usage: genomecoverage_single.sh <filtered.bam> <genome> <output_dir> [PE|SE]
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
BAM="${1:?BAM required}"; GENOME="${2:?genome required}"; OUT="${3:?output directory required}"
REQUESTED_LAYOUT="${4:-}"
SAMPLE="$(basename "$BAM" .bam)"; THREADS="${THREADS_BIGWIG:-2}"; BIN_SIZE="${TRACK_BIN_SIZE:-10}"
HELPERS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/track_normalization_helpers.sh"
[[ -f "$HELPERS" ]] || { echo "ERROR: missing track helper: $HELPERS" >&2; exit 1; }
# shellcheck disable=SC1090
source "$HELPERS"
mkdir -p "$OUT"; samtools quickcheck "$BAM"
[[ "${TRACK_STANDARD_CHROMS_ONLY:-true}" == "true" ]] || {
    echo "ERROR: TRACK_STANDARD_CHROMS_ONLY=false is incompatible with the unified canonical track universe" >&2
    exit 1
}
LAYOUT="$(resolve_signal_layout "$BAM" "$REQUESTED_LAYOUT")"
SIGNAL_UNIT="$(signal_unit_for_layout "$LAYOUT")"
declare -a SIGNAL_ARGS=()
deeptools_signal_args "$LAYOUT" SIGNAL_ARGS
tmp_dir="$(mktemp -d "${OUT}/.${SAMPLE}.cpm.XXXXXX")"; trap 'rm -rf "$tmp_dir"' EXIT
mapfile -t contigs < <(canonical_contigs_from_bam "$BAM" "$GENOME")
(( ${#contigs[@]} > 0 )) || { echo "ERROR: no canonical chromosomes found in $BAM for $GENOME" >&2; exit 1; }
track_bam="${tmp_dir}/${SAMPLE}.canonical.bam"
samtools view -@ "$THREADS" -b -o "$track_bam" "$BAM" "${contigs[@]}"
samtools index -@ "$THREADS" "$track_bam"
SIGNAL_COUNT="$(signal_count_for_bam "$track_bam" "$LAYOUT")"
(( SIGNAL_COUNT > 0 )) || { echo "ERROR: no canonical $SIGNAL_UNIT records found in $BAM" >&2; exit 1; }
GENERATE_BW="${GENERATE_COVERAGE_BIGWIGS:-true}"
GENERATE_BG="${GENERATE_COVERAGE_BEDGRAPHS:-true}"
[[ "$GENERATE_BW" == "true" || "$GENERATE_BG" == "true" ]] || {
    echo "ERROR: both GENERATE_COVERAGE_BIGWIGS and GENERATE_COVERAGE_BEDGRAPHS are false" >&2
    exit 1
}
out_bw="NA"; out_bg="NA"
if [[ "$GENERATE_BW" == "true" ]]; then
    out_bw="${OUT}/${SAMPLE}_CPM.bw"
    bamCoverage --bam "$track_bam" --outFileName "$out_bw" --outFileFormat bigwig \
        --normalizeUsing CPM --exactScaling --binSize "$BIN_SIZE" --numberOfProcessors "$THREADS" \
        "${SIGNAL_ARGS[@]}" ${BAMCOVERAGE_COMMON_ARGS:-}
    [[ -s "$out_bw" ]] || { echo "ERROR: CPM bigWig was not created: $out_bw" >&2; exit 1; }
fi
if [[ "$GENERATE_BG" == "true" ]]; then
    out_bg="${OUT}/${SAMPLE}_CPM.bedGraph"
    bamCoverage --bam "$track_bam" --outFileName "$out_bg" --outFileFormat bedgraph \
        --normalizeUsing CPM --exactScaling --binSize "$BIN_SIZE" --numberOfProcessors "$THREADS" \
        "${SIGNAL_ARGS[@]}" ${BAMCOVERAGE_COMMON_ARGS:-}
    [[ -s "$out_bg" ]] || { echo "ERROR: CPM bedGraph was not created: $out_bg" >&2; exit 1; }
fi
printf 'sample\tgenome\tlayout\tsignal_unit\tnormalization\tnormalization_count\tbin_size\tcanonical_contigs\tbigwig_file\tbedgraph_file\n%s\t%s\t%s\t%s\tCPM\t%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$GENOME" "$LAYOUT" "$SIGNAL_UNIT" "$SIGNAL_COUNT" "$BIN_SIZE" \
    "$(IFS=,; echo "${contigs[*]}")" "$out_bw" "$out_bg" > "${OUT}/${SAMPLE}_CPM.metadata.tsv"
echo "CPM tracks ($SIGNAL_UNIT-based): bigWig=$out_bw bedGraph=$out_bg"
