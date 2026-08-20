#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - optional HOMER annotation and motif enrichment
# Usage: peak_interpretation.sh <consensus_peaks.bed> <hg38|mm39> <output_dir>
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
PEAKS="${1:?consensus BED required}"; GENOME="${2:?genome required}"; OUT="${3:?output directory required}"
[[ -s "$PEAKS" ]] || { echo "ERROR: consensus BED missing or empty: $PEAKS" >&2; exit 1; }
mkdir -p "$OUT"
case "$GENOME" in
    hg38) HOMER_GENOME="${HOMER_GENOME_HG38:-hg38}"; GTF="${GTF_HUMAN:-}" ;;
    mm39) HOMER_GENOME="${HOMER_GENOME_MM39:-mm39}"; GTF="${GTF_MOUSE:-}" ;;
    *) echo "ERROR: unsupported genome: $GENOME" >&2; exit 1 ;;
esac
if [[ "${RUN_PEAK_ANNOTATION:-false}" == "true" ]]; then
    command -v annotatePeaks.pl >/dev/null 2>&1 || { echo "ERROR: annotatePeaks.pl is not installed" >&2; exit 1; }
    annotation_args=(); [[ -s "$GTF" ]] && annotation_args=( -gtf "$GTF" )
    annotatePeaks.pl "$PEAKS" "$HOMER_GENOME" "${annotation_args[@]}" > "${OUT}/consensus_peak_annotation.tsv"
fi
if [[ "${RUN_MOTIF_ENRICHMENT:-false}" == "true" ]]; then
    command -v findMotifsGenome.pl >/dev/null 2>&1 || { echo "ERROR: findMotifsGenome.pl is not installed" >&2; exit 1; }
    findMotifsGenome.pl "$PEAKS" "$HOMER_GENOME" "${OUT}/motifs" -size given -p "${THREADS_MOTIF:-4}"
fi
