#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - UCSC custom-track definitions for bigWig outputs
# Usage: create_ucsc_tracks.sh <bigwig_dir> [public_url_base] [description_prefix]
set -euo pipefail
BIGWIG_DIR="${1:?bigWig directory required}"; URL_BASE="${2:-}"; PREFIX="${3:-ATAC-seq}"
[[ -d "$BIGWIG_DIR" ]] || { echo "ERROR: directory not found: $BIGWIG_DIR" >&2; exit 1; }
TRACK_FILE="${BIGWIG_DIR}/ucsc_tracks.txt"; : > "$TRACK_FILE"; count=0
while IFS= read -r bw; do
    file="$(basename "$bw")"; name="${file%.bw}"; url="$file"
    [[ -n "$URL_BASE" ]] && url="${URL_BASE%/}/${file}"
    if [[ "$name" == *_DESeq2RobustCPM || "$name" == *_DESeq2RobustCPM_* ]]; then
        description="${PREFIX} DESeq2 robust CPM"; color="120,50,160"
    elif [[ "$name" == *_DESeq2Consensus ]]; then
        description="${PREFIX} DESeq2 consensus-peak normalized"; color="180,50,50"
    elif [[ "$name" == *_CPM ]]; then
        description="${PREFIX} counts per million"; color="40,90,180"
    else
        description="${PREFIX} bigWig"; color="80,80,80"
    fi
    printf 'track type=bigWig name="%s" description="%s: %s" bigDataUrl=%s visibility=full autoScale=on alwaysZero=on color=%s\n' \
        "$name" "$description" "$name" "$url" "$color" >> "$TRACK_FILE"
    count=$((count + 1))
done < <(find "$BIGWIG_DIR" -maxdepth 1 -type f -name '*.bw' | sort)
(( count > 0 )) || { echo "ERROR: no bigWig files found in $BIGWIG_DIR" >&2; exit 1; }
echo "UCSC custom tracks: $TRACK_FILE ($count tracks)"
