#!/usr/bin/env bash
# Split a deduplicated host+dm6 composite BAM and apply the v4.2.0 stringent policy.
# Usage: filter_composite_spikein_bam.sh <dedup.bam> <host_genome> <host_blacklist> <dm6_blacklist> <host_out.bam> <dm6_out.bam> <PE|SE>
set -euo pipefail

[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || {
    echo "ERROR: F2T_CONFIG is not set" >&2
    exit 1
}
# shellcheck disable=SC1090
source "$F2T_CONFIG"

INPUT_BAM="${1:?deduplicated composite BAM required}"
HOST_GENOME="${2:?host genome required}"
HOST_BLACKLIST="${3:?host blacklist required}"
DM6_BLACKLIST="${4:?dm6 blacklist required}"
HOST_OUTPUT="${5:?host output BAM required}"
DM6_OUTPUT="${6:?dm6 output BAM required}"
LAYOUT="${7:?PE or SE required}"
THREADS="${THREADS_SAMTOOLS:-2}"
MINIMUM_MAPQ="${SPIKEIN_MIN_MAPQ:-30}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1090
source "${SCRIPT_DIR}/track_normalization_helpers.sh"

[[ "$LAYOUT" == "PE" || "$LAYOUT" == "SE" ]] || { echo "ERROR: layout must be PE or SE" >&2; exit 1; }
[[ "$MINIMUM_MAPQ" =~ ^[0-9]+$ ]] || { echo "ERROR: SPIKEIN_MIN_MAPQ must be non-negative" >&2; exit 1; }
[[ -s "$INPUT_BAM" ]] || { echo "ERROR: composite BAM missing: $INPUT_BAM" >&2; exit 1; }
[[ -s "$HOST_BLACKLIST" ]] || { echo "ERROR: host blacklist missing: $HOST_BLACKLIST" >&2; exit 1; }
[[ -s "$DM6_BLACKLIST" ]] || { echo "ERROR: dm6 blacklist missing: $DM6_BLACKLIST" >&2; exit 1; }

mkdir -p "$(dirname "$HOST_OUTPUT")" "$(dirname "$DM6_OUTPUT")"
samtools quickcheck "$INPUT_BAM"
[[ -s "${INPUT_BAM}.bai" ]] || samtools index -@ "$THREADS" "$INPUT_BAM"

tmp_dir="$(mktemp -d "$(dirname "$HOST_OUTPUT")/.spikein-filter.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT
prefilter="${tmp_dir}/prefilter.bam"

if [[ "$LAYOUT" == "PE" ]]; then
    samtools view -@ "$THREADS" -b -q "$MINIMUM_MAPQ" -f 2 -F 3852 \
        -o "$prefilter" "$INPUT_BAM"
else
    samtools view -@ "$THREADS" -b -q "$MINIMUM_MAPQ" -F 3844 \
        -o "$prefilter" "$INPUT_BAM"
fi
samtools index -@ "$THREADS" "$prefilter"

mapfile -t host_contigs < <(canonical_contigs_from_bam "$prefilter" "$HOST_GENOME")
(( ${#host_contigs[@]} > 0 )) || { echo "ERROR: no canonical host contigs remain" >&2; exit 1; }

IFS=',' read -r -a configured_dm6 <<< "${SPIKEIN_CANONICAL_CONTIGS:-2L,2R,3L,3R,4,X}"
dm6_contigs=()
header_contigs="$(samtools view -H "$prefilter" | awk -F '\t' '$1=="@SQ" {for(i=1;i<=NF;i++) if($i~/^SN:/){sub(/^SN:/,"",$i); print $i}}')"
for contig in "${configured_dm6[@]}"; do
    contig="${contig//[[:space:]]/}"
    [[ -n "$contig" ]] || continue
    prefixed="dm6__${contig}"
    [[ "$contig" == chr* ]] || prefixed="dm6__chr${contig}"
    if grep -Fxq "$prefixed" <<< "$header_contigs"; then
        dm6_contigs+=("$prefixed")
    else
        echo "ERROR: configured dm6 contig is absent from composite BAM: $prefixed" >&2
        exit 1
    fi
done
(( ${#dm6_contigs[@]} > 0 )) || { echo "ERROR: no dm6 canonical contigs configured" >&2; exit 1; }

samtools view -@ "$THREADS" -b -o "${tmp_dir}/host.canonical.bam" "$prefilter" "${host_contigs[@]}"
samtools view -@ "$THREADS" -b -o "${tmp_dir}/dm6.prefixed.bam" "$prefilter" "${dm6_contigs[@]}"

printf '%s\n' "${host_contigs[@]}" > "${tmp_dir}/host.keep.txt"
printf '%s\n' "${dm6_contigs[@]}" > "${tmp_dir}/dm6.keep.txt"
samtools view -H "${tmp_dir}/host.canonical.bam" \
    | awk -F '\t' 'BEGIN{OFS="\t"}
        NR==FNR {keep[$1]=1; next}
        $1=="@SQ" {
            name=""
            for(i=1;i<=NF;i++) if($i~/^SN:/){name=$i; sub(/^SN:/,"",name)}
            if(keep[name]) print
            next
        }
        {print}
    ' "${tmp_dir}/host.keep.txt" - > "${tmp_dir}/host.header.sam"
samtools reheader "${tmp_dir}/host.header.sam" "${tmp_dir}/host.canonical.bam" \
    > "${tmp_dir}/host.header_filtered.bam"

samtools view -H "${tmp_dir}/dm6.prefixed.bam" \
    | awk -F '\t' 'BEGIN{OFS="\t"}
        NR==FNR {keep[$1]=1; next}
        $1=="@SQ" {
            name=""
            for(i=1;i<=NF;i++) if($i~/^SN:/){name=$i; sub(/^SN:/,"",name)}
            if(keep[name]) {
                for(i=1;i<=NF;i++) sub(/^SN:dm6__/,"SN:",$i)
                print
            }
            next
        }
        {print}
    ' "${tmp_dir}/dm6.keep.txt" - > "${tmp_dir}/dm6.header.sam"
samtools reheader "${tmp_dir}/dm6.header.sam" "${tmp_dir}/dm6.prefixed.bam" \
    > "${tmp_dir}/dm6.canonical.bam"

pair_safe_blacklist_filter() {
    local input="$1" blacklist="$2" output="$3" label="$4"
    bedtools intersect -v -abam "$input" -b "$blacklist" > "${tmp_dir}/${label}.no_blacklist.bam"
    if [[ "$LAYOUT" == "PE" ]]; then
        samtools sort -n -@ "$THREADS" -o "${tmp_dir}/${label}.name.bam" "${tmp_dir}/${label}.no_blacklist.bam"
        samtools fixmate -r -@ "$THREADS" "${tmp_dir}/${label}.name.bam" "${tmp_dir}/${label}.fixmate.bam"
        samtools view -@ "$THREADS" -b -f 2 -o "${tmp_dir}/${label}.paired.bam" "${tmp_dir}/${label}.fixmate.bam"
        samtools sort -@ "$THREADS" -o "$output" "${tmp_dir}/${label}.paired.bam"
    else
        samtools sort -@ "$THREADS" -o "$output" "${tmp_dir}/${label}.no_blacklist.bam"
    fi
}

pair_safe_blacklist_filter "${tmp_dir}/host.header_filtered.bam" "$HOST_BLACKLIST" "${tmp_dir}/host.final.bam" host
pair_safe_blacklist_filter "${tmp_dir}/dm6.canonical.bam" "$DM6_BLACKLIST" "${tmp_dir}/dm6.final.bam" dm6

for bam in "${tmp_dir}/host.final.bam" "${tmp_dir}/dm6.final.bam"; do
    samtools quickcheck "$bam"
    samtools index -@ "$THREADS" "$bam"
    count="$(signal_count_for_bam "$bam" "$LAYOUT")"
    (( count > 0 )) || { echo "ERROR: zero retained $(basename "$bam") observations" >&2; exit 1; }
done

mv "${tmp_dir}/host.final.bam" "$HOST_OUTPUT"
mv "${tmp_dir}/host.final.bam.bai" "${HOST_OUTPUT}.bai"
mv "${tmp_dir}/dm6.final.bam" "$DM6_OUTPUT"
mv "${tmp_dir}/dm6.final.bam.bai" "${DM6_OUTPUT}.bai"
echo "Composite stringent filtering complete: host=$HOST_OUTPUT dm6=$DM6_OUTPUT"
