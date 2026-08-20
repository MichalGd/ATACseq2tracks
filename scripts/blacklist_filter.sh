#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - usable-read filtering plus pair-safe blacklist removal
# Usage: blacklist_filter.sh <input.bam> <blacklist.bed> <out_dir> <PE|SE> <genome>
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"

INPUT_BAM="${1:?input BAM required}"; BLACKLIST_BED="${2:?blacklist BED required}"; OUT_DIR="${3:?output directory required}"
LAYOUT="${4:?PE or SE required}"; GENOME="${5:?genome required}"
SAMPLE="$(basename "$INPUT_BAM" .bam)"; FILTERED="${OUT_DIR}/${SAMPLE}_blFilt.bam"
THREADS="${THREADS_SAMTOOLS:-2}"; MIN_MAPQ="${MIN_MAPQ:-30}"
mkdir -p "$OUT_DIR"
[[ -s "$INPUT_BAM" ]] || { echo "ERROR: input BAM missing or empty: $INPUT_BAM" >&2; exit 1; }
[[ -s "$BLACKLIST_BED" ]] || { echo "ERROR: blacklist missing or empty: $BLACKLIST_BED" >&2; exit 1; }
samtools quickcheck "$INPUT_BAM"

tmp_dir="$(mktemp -d "${OUT_DIR}/.${SAMPLE}.filter.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT
prefilter="${tmp_dir}/prefilter.bam"; nuclear="${tmp_dir}/nuclear.bam"; no_blacklist="${tmp_dir}/no_blacklist.bam"

before="$(samtools view -@ "$THREADS" -c "$INPUT_BAM")"
if [[ "$LAYOUT" == "PE" ]]; then
    # Proper pairs; exclude unmapped/mate-unmapped, secondary, QC-fail, duplicate and supplementary records.
    samtools view -@ "$THREADS" -b -q "$MIN_MAPQ" -f 2 -F 3852 -o "$prefilter" "$INPUT_BAM"
else
    samtools view -@ "$THREADS" -b -q "$MIN_MAPQ" -F 3844 -o "$prefilter" "$INPUT_BAM"
fi
after_flags="$(samtools view -@ "$THREADS" -c "$prefilter")"

if [[ "${REMOVE_MITO:-true}" == "true" ]]; then
    samtools index -@ "$THREADS" "$prefilter"
    mapfile -t contigs < <(samtools idxstats "$prefilter" | awk '$1!="*" && $1!~/^(chrM|MT|M)$/ {print $1}')
    (( ${#contigs[@]} > 0 )) || { echo "ERROR: no nuclear contigs remain for $SAMPLE" >&2; exit 1; }
    samtools view -@ "$THREADS" -b -o "$nuclear" "$prefilter" "${contigs[@]}"
else
    cp "$prefilter" "$nuclear"
fi
after_mito="$(samtools view -@ "$THREADS" -c "$nuclear")"

bedtools intersect -v -abam "$nuclear" -b "$BLACKLIST_BED" > "$no_blacklist"
if [[ "$LAYOUT" == "PE" ]]; then
    samtools sort -n -@ "$THREADS" -o "${tmp_dir}/name.bam" "$no_blacklist"
    samtools fixmate -r -@ "$THREADS" "${tmp_dir}/name.bam" "${tmp_dir}/fixmate.bam"
    samtools view -@ "$THREADS" -b -f 2 -o "${tmp_dir}/paired.bam" "${tmp_dir}/fixmate.bam"
    samtools sort -@ "$THREADS" -o "${tmp_dir}/final.bam" "${tmp_dir}/paired.bam"
else
    samtools sort -@ "$THREADS" -o "${tmp_dir}/final.bam" "$no_blacklist"
fi

after="$(samtools view -@ "$THREADS" -c "${tmp_dir}/final.bam")"
if (( after == 0 )) && [[ "${ALLOW_EMPTY_FILTERED_BAM:-false}" != "true" ]]; then
    echo "ERROR: filtering removed every alignment for $SAMPLE" >&2; exit 1
fi
samtools quickcheck "${tmp_dir}/final.bam"; samtools index -@ "$THREADS" "${tmp_dir}/final.bam"
mv "${tmp_dir}/final.bam" "$FILTERED"; mv "${tmp_dir}/final.bam.bai" "${FILTERED}.bai"
printf 'sample\tgenome\tlayout\tinput_alignments\tafter_flags_mapq\tafter_mito\tafter_blacklist_pair_cleanup\n' > "${OUT_DIR}/${SAMPLE}_filter_metrics.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$GENOME" "$LAYOUT" "$before" "$after_flags" "$after_mito" "$after" >> "${OUT_DIR}/${SAMPLE}_filter_metrics.tsv"
echo "Filtered: $SAMPLE before=$before flags_mapq=$after_flags mito=$after_mito final=$after"
