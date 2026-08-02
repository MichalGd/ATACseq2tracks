#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - fragment-aware RPM/CPM bigWig for UCSC/IGV
# Usage: genomecoverage_single.sh <filtered.bam> <genome> <output_dir>
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
BAM="${1:?BAM required}"; GENOME="${2:?genome required}"; OUT="${3:?output directory required}"
SAMPLE="$(basename "$BAM" .bam)"; THREADS="${THREADS_BIGWIG:-2}"; BIN_SIZE="${TRACK_BIN_SIZE:-10}"
mkdir -p "$OUT"; samtools quickcheck "$BAM"
tmp_dir="$(mktemp -d "${OUT}/.${SAMPLE}.rpm.XXXXXX")"; trap 'rm -rf "$tmp_dir"' EXIT
track_bam="$BAM"
if [[ "${TRACK_STANDARD_CHROMS_ONLY:-true}" == "true" ]]; then
    mapfile -t contigs < <(samtools idxstats "$BAM" | awk '$1~/^chr([0-9]+|X|Y)$/ {print $1}')
    (( ${#contigs[@]} > 0 )) || { echo "ERROR: no standard chromosomes found in $BAM" >&2; exit 1; }
    track_bam="${tmp_dir}/${SAMPLE}.standard.bam"
    samtools view -@ "$THREADS" -b -o "$track_bam" "$BAM" "${contigs[@]}"; samtools index -@ "$THREADS" "$track_bam"
fi
out_bw="${OUT}/${SAMPLE}_RPM.bw"
bamCoverage --bam "$track_bam" --outFileName "$out_bw" --outFileFormat bigwig \
    --normalizeUsing CPM --binSize "$BIN_SIZE" --numberOfProcessors "$THREADS" \
    ${BAMCOVERAGE_COMMON_ARGS:-}
[[ -s "$out_bw" ]] || { echo "ERROR: RPM bigWig was not created: $out_bw" >&2; exit 1; }
printf 'sample\tgenome\tnormalization\tbin_size\tfile\n%s\t%s\tCPM/RPM\t%s\t%s\n' \
    "$SAMPLE" "$GENOME" "$BIN_SIZE" "$out_bw" > "${OUT}/${SAMPLE}_RPM.metadata.tsv"
echo "RPM bigWig: $out_bw"
