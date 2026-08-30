#!/usr/bin/env bash
set -euo pipefail

OUTPUT_DIR="${1:?usage: cleanup_intermediates.sh OUTPUT_DIR}"
if [[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]]; then
    # F2T_CONFIG is the literal-only file produced by resolve_config.py.
    # shellcheck disable=SC1090
    source "$F2T_CONFIG"
fi
[[ -d "$OUTPUT_DIR" ]] || { echo "ERROR: output directory not found: $OUTPUT_DIR" >&2; exit 1; }
MANIFEST="${OUTPUT_DIR}/metadata/cleanup_manifest.tsv"
mkdir -p "${OUTPUT_DIR}/metadata"
[[ -s "$MANIFEST" ]] || printf 'timestamp\tpath\tsize_bytes\taction\n' > "$MANIFEST"

remove_matches() {
    local directory="$1"; shift
    [[ -d "$directory" ]] || return 0
    local file size
    while IFS= read -r -d '' file; do
        size="$(stat -c '%s' "$file" 2>/dev/null || echo 0)"
        printf '%s\t%s\t%s\tdeleted\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$file" "$size" >> "$MANIFEST"
        rm -f -- "$file"
    done < <(find "$directory" -maxdepth 1 -type f \( "$@" \) -print0)
}

[[ "${KEEP_INTERMEDIATE_BAMS:-false}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/bams" -name '*.bam' -o -name '*.bai'
[[ "${KEEP_TRIMMED_FASTQ:-false}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/trimmedFastq" -name '*.gz'
[[ "${KEEP_DEDUP_BAMS:-false}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/dedupBams" -name '*.bam' -o -name '*.bai'
[[ "${KEEP_FILTERED_BAMS:-true}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/filteredBams" -name '*.bam' -o -name '*.bai'
[[ "${KEEP_NORMALIZATION_POLICY_BAMS:-false}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/coverage_filtering_policy_bams" -name '*.bam' -o -name '*.bai'
if [[ "${KEEP_NORMALIZATION_POLICY_BAMS:-false}" == "false" ]]; then
    for directory in "${OUTPUT_DIR}/coverage_filtering_policy_bams"/*; do
        [[ -d "$directory" ]] && remove_matches "$directory" -name '*.bam' -o -name '*.bai'
    done
fi
if [[ "${KEEP_SPIKEIN_BAMS:-false}" == "false" ]]; then
    for directory in "${OUTPUT_DIR}/spikein" "${OUTPUT_DIR}/spikein"/*; do
        [[ -d "$directory" ]] && remove_matches "$directory" -name '*.bam' -o -name '*.bai'
    done
fi
[[ "${KEEP_RAW_BEDGRAPH:-false}" == "false" ]] && \
    remove_matches "${OUTPUT_DIR}/bedGraph" -name '*.bedGraph.gz'

echo "Cleanup manifest: $MANIFEST"
