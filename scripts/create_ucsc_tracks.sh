#!/bin/bash
# ATACseq2tracks v3.0.2 — UCSC BigWig track line generator (no config needed)
# Usage: bash scripts/create_ucsc_tracks.sh <bigwigDir> [url_base]
set -euo pipefail
OUT="$1"; URL_BASE="${2:-http://your-server.com/data}"
TRACK_FILE="${OUT}/ucsc_tracks.txt"
> "$TRACK_FILE"
find "$OUT" -name "*.bw" | sort | while read -r bw; do
    name=$(basename "$bw" .bw)
    echo "track type=bigWig name=\"${name}\" description=\"${name}\" bigDataUrl=${URL_BASE}/$(basename "$bw") visibility=full autoScale=on" >> "$TRACK_FILE"
done
echo "UCSC track lines written to: $TRACK_FILE"
