#!/usr/bin/env bash
# Build one host+dm6 Bowtie2 composite reference without altering source files.
# Usage: prepare_dm6_spikein_references.sh <hg38|mm39> <host_index_prefix> <dm6_fasta[.gz]> <dm6_chrom_sizes> <dm6_blacklist> <shared_reference_root>
set -euo pipefail

HOST_GENOME="${1:?host genome hg38 or mm39 required}"
HOST_INDEX="${2:?host Bowtie2 index prefix required}"
DM6_FASTA="${3:?dm6 FASTA required}"
DM6_CHROM_SIZES="${4:?dm6 chromosome sizes required}"
DM6_BLACKLIST="${5:?dm6 blacklist required}"
REFERENCE_ROOT="${6:?shared reference root required}"

[[ "$HOST_GENOME" == hg38 || "$HOST_GENOME" == mm39 ]] || { echo "ERROR: host must be hg38 or mm39" >&2; exit 1; }
[[ "$REFERENCE_ROOT" == /* && "$REFERENCE_ROOT" != / ]] || {
    echo "ERROR: shared reference root must be an absolute, non-root path" >&2
    exit 1
}
for tool in bowtie2-inspect bowtie2-build sha256sum awk grep sed install; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: tool not in PATH: $tool" >&2; exit 1; }
done
for source in "$DM6_FASTA" "$DM6_CHROM_SIZES" "$DM6_BLACKLIST"; do
    [[ -s "$source" ]] || { echo "ERROR: source missing or empty: $source" >&2; exit 1; }
done

index_complete() {
    local prefix="$1" extension component missing
    for extension in bt2 bt2l; do
        missing=0
        for component in 1 2 3 4; do [[ -s "${prefix}.${component}.${extension}" ]] || missing=$((missing + 1)); done
        for component in 1 2; do [[ -s "${prefix}.rev.${component}.${extension}" ]] || missing=$((missing + 1)); done
        (( missing == 0 )) && return 0
    done
    return 1
}
index_complete "$HOST_INDEX" || { echo "ERROR: incomplete host Bowtie2 index: $HOST_INDEX" >&2; exit 1; }

DM6_ROOT="${REFERENCE_ROOT}/dm6"
COMPOSITE_ROOT="${REFERENCE_ROOT}/composite/${HOST_GENOME}_dm6"
INDEX_DIR="${COMPOSITE_ROOT}/bowtie2"
INDEX_PREFIX="${INDEX_DIR}/${HOST_GENOME}_dm6"
if index_complete "$INDEX_PREFIX" && [[ "${OVERWRITE:-false}" != true ]]; then
    echo "Composite index already complete; leaving unchanged: $INDEX_PREFIX"
    exit 0
fi
if [[ -e "$COMPOSITE_ROOT" && "${OVERWRITE:-false}" != true ]]; then
    echo "ERROR: incomplete destination exists; inspect it or rerun with OVERWRITE=true: $COMPOSITE_ROOT" >&2
    exit 1
fi

mkdir -p "$DM6_ROOT" "$INDEX_DIR"
tmp_dir="$(mktemp -d "${REFERENCE_ROOT}/.dm6-composite.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT
host_fasta="${tmp_dir}/${HOST_GENOME}.fa"
dm6_prefixed="${tmp_dir}/dm6.prefixed.fa"
composite_fasta="${tmp_dir}/${HOST_GENOME}_dm6.fa"

bowtie2-inspect "$HOST_INDEX" > "$host_fasta"
if [[ "$DM6_FASTA" == *.gz ]]; then
    gzip -dc "$DM6_FASTA"
else
    sed -n 'p' "$DM6_FASTA"
fi | awk '
    /^>/ {
        sub(/^>/, "")
        sub(/[[:space:]].*$/, "")
        name=$0
        if (name !~ /^chr/) name="chr" name
        print ">dm6__" name
        next
    }
    {print}
' > "$dm6_prefixed"
for contig in 2L 2R 3L 3R 4 X; do
    grep -q "^>dm6__chr${contig}$" "$dm6_prefixed" || {
        echo "ERROR: dm6 FASTA lacks required canonical contig: chr${contig}" >&2
        exit 1
    }
done

cp "$host_fasta" "$composite_fasta"
printf '\n' >> "$composite_fasta"
sed -n 'p' "$dm6_prefixed" >> "$composite_fasta"

if [[ "${OVERWRITE:-false}" == true ]]; then
    rm -f "${INDEX_PREFIX}."*.bt2 "${INDEX_PREFIX}."*.bt2l "${INDEX_PREFIX}.rev."*.bt2 "${INDEX_PREFIX}.rev."*.bt2l
fi
bowtie2-build --threads "${BOWTIE2_BUILD_THREADS:-8}" "$composite_fasta" "$INDEX_PREFIX"
index_complete "$INDEX_PREFIX" || { echo "ERROR: composite index build is incomplete" >&2; exit 1; }

install -m 0644 "$DM6_CHROM_SIZES" "${DM6_ROOT}/dm6.chrom.sizes"
install -m 0644 "$DM6_BLACKLIST" "${DM6_ROOT}/dm6-blacklist.v2.bed"
install -m 0644 "$composite_fasta" "${COMPOSITE_ROOT}/${HOST_GENOME}_dm6.fa"

MANIFEST="${COMPOSITE_ROOT}/reference_manifest.tsv"
printf 'resource\tpath\tsha256\n' > "$MANIFEST"
for resource in "${COMPOSITE_ROOT}/${HOST_GENOME}_dm6.fa" "${DM6_ROOT}/dm6.chrom.sizes" \
    "${DM6_ROOT}/dm6-blacklist.v2.bed" "${INDEX_PREFIX}."*.bt2 "${INDEX_PREFIX}."*.bt2l; do
    [[ -f "$resource" ]] || continue
    printf '%s\t%s\t%s\n' "$(basename "$resource")" "$resource" "$(sha256sum "$resource" | awk '{print $1}')" >> "$MANIFEST"
done
chmod 0755 "$REFERENCE_ROOT" "$DM6_ROOT" "${REFERENCE_ROOT}/composite" "$COMPOSITE_ROOT" "$INDEX_DIR"
find "$COMPOSITE_ROOT" "$DM6_ROOT" -type f -exec chmod 0644 {} +
echo "Prepared shared composite reference: $INDEX_PREFIX"
echo "Manifest: $MANIFEST"
