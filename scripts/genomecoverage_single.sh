#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - canonical-chromosome fragment/read CPM bigWig
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
out_bw="${OUT}/${SAMPLE}_CPM.bw"
bamCoverage --bam "$track_bam" --outFileName "$out_bw" --outFileFormat bigwig \
    --normalizeUsing CPM --exactScaling --binSize "$BIN_SIZE" --numberOfProcessors "$THREADS" \
    "${SIGNAL_ARGS[@]}" \
    ${BAMCOVERAGE_COMMON_ARGS:-}
[[ -s "$out_bw" ]] || { echo "ERROR: CPM bigWig was not created: $out_bw" >&2; exit 1; }
printf 'sample\tgenome\tlayout\tsignal_unit\tnormalization\tnormalization_count\tbin_size\tcanonical_contigs\tfile\n%s\t%s\t%s\t%s\tCPM\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$GENOME" "$LAYOUT" "$SIGNAL_UNIT" "$SIGNAL_COUNT" "$BIN_SIZE" \
    "$(IFS=,; echo "${contigs[*]}")" "$out_bw" > "${OUT}/${SAMPLE}_CPM.metadata.tsv"
echo "CPM bigWig ($SIGNAL_UNIT-based): $out_bw"
