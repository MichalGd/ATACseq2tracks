#!/usr/bin/env bash
set -euo pipefail

METADATA_DIR="${1:?usage: capture_provenance.sh METADATA_DIR RESOLVED_CONFIG_TSV}"
CONFIG_TSV="${2:?usage: capture_provenance.sh METADATA_DIR RESOLVED_CONFIG_TSV}"
mkdir -p "$METADATA_DIR"

SOFTWARE="${METADATA_DIR}/software_versions.tsv"
printf 'tool\tpath\tversion\n' > "$SOFTWARE"
for tool in python3 Rscript bowtie2 samtools bedtools trim_galore fastqc macs3 multiqc bamCoverage ataqv mkarv; do
    if command -v "$tool" >/dev/null 2>&1; then
        path="$(command -v "$tool")"
        version="$($tool --version 2>&1 | head -n 1 || true)"
        printf '%s\t%s\t%s\n' "$tool" "$path" "${version//$'\t'/ }" >> "$SOFTWARE"
    else
        printf '%s\t%s\t%s\n' "$tool" "NOT_FOUND" "NOT_FOUND" >> "$SOFTWARE"
    fi
done

REFERENCES="${METADATA_DIR}/reference_manifest.tsv"
printf 'config_key\tconfigured_path\tstatus\n' > "$REFERENCES"
awk -F '\t' 'NR > 1 && $1 ~ /^(INDEX_|GTF_|CHROM_SIZES_|BLACKLIST_|CCRE_BED_|PICARD_JAR|BEDGRAPH_TO_BIGWIG)/ {print $1 "\t" $2}' "$CONFIG_TSV" |
while IFS=$'\t' read -r key path; do
    status="missing"
    [[ -e "$path" ]] && status="present"
    if [[ "$key" == INDEX_* && "$status" == "missing" ]]; then
        for suffix in 1.bt2 1.bt2l; do
            [[ -s "${path}.${suffix}" ]] && { status="index_prefix_present"; break; }
        done
    fi
    printf '%s\t%s\t%s\n' "$key" "$path" "$status" >> "$REFERENCES"
done
