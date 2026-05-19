#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — MACS2 peak calling (single sample, narrow or broad)
# Usage: bash scripts/macs2_peaks.sh <ip.bam> <ctrl.bam|none> <outDir> <mode> <genome_key> <sample_name>
#   mode: narrow | broad | both
#   genome_key: hg38 | mm39
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

IP_BAM="$1"
CTRL_BAM="$2"       # use "none" if no control
OUT_DIR="$3"
MODE="${4:-narrow}"
GENOME_KEY="${5:-hg38}"
SAMPLE="${6:-$(basename "$IP_BAM" .bam)}"

mkdir -p "$OUT_DIR"

# genome size flag
[[ "$GENOME_KEY" == "hg38" ]] && GSIZE="$MACS2_GENOME_HG38" || GSIZE="$MACS2_GENOME_MM39"

# control flag
CTRL_FLAG=""
[[ "$CTRL_BAM" != "none" && -f "$CTRL_BAM" ]] && CTRL_FLAG="-c $CTRL_BAM"

run_macs2() {
    local peak_type="$1"  # narrowPeak or broadPeak
    local extra_flags="$2"
    echo "[MACS2] $SAMPLE  mode=$peak_type  control=$([ -z "$CTRL_FLAG" ] && echo none || echo $CTRL_BAM)"
    macs2 callpeak \
        -t "$IP_BAM" \
        $CTRL_FLAG \
        -f BAM \
        -g "$GSIZE" \
        -n "$SAMPLE" \
        --outdir "$OUT_DIR" \
        -q "${MACS2_QVALUE}" \
        $extra_flags \
        2>"${OUT_DIR}/${SAMPLE}_macs2_${peak_type}.log"
}

case "$MODE" in
    narrow)
        run_macs2 "narrowPeak" ""
        ;;
    broad)
        run_macs2 "broadPeak" "--broad --broad-cutoff ${MACS2_BROAD_CUTOFF}"
        ;;
    both)
        run_macs2 "narrowPeak" ""
        run_macs2 "broadPeak"  "--broad --broad-cutoff ${MACS2_BROAD_CUTOFF}"
        ;;
    none)
        echo "MACS2 skipped for $SAMPLE (mode=none)"
        ;;
    *)
        echo "ERROR: Unknown mode '$MODE'" >&2; exit 1
        ;;
esac
