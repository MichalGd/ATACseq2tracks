#!/usr/bin/env bash
# Example only. Edit paths before running.
set -euo pipefail

INPUT=/data/project/raw_fastq_human
OUTPUT=/data/project/ATACseq2tracks_human
MAX_JOBS=4
SPECIES=human

../atacseq2tracks.sh "$INPUT" "$OUTPUT" "$MAX_JOBS" "$SPECIES"
