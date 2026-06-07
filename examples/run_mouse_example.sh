#!/usr/bin/env bash
# Example only. Edit paths before running.
set -euo pipefail

INPUT=/data/project/raw_fastq_mouse
OUTPUT=/data/project/ATACseq2tracks_mouse
MAX_JOBS=4
SPECIES=mouse

../atacseq2tracks.sh "$INPUT" "$OUTPUT" "$MAX_JOBS" "$SPECIES"
