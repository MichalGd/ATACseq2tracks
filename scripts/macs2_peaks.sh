#!/usr/bin/env bash
# ATACseq2tracks v4.1.0 - strict, layout-aware MACS3 peak calling
# Usage: macs2_peaks.sh <sample.bam> <control.bam|none> <out_dir> <narrow|broad|both|none> <genome> <sample_name> <PE|SE>
set -euo pipefail
[[ -n "${F2T_CONFIG:-}" && -f "$F2T_CONFIG" ]] || { echo "ERROR: F2T_CONFIG is not set" >&2; exit 1; }
# shellcheck disable=SC1090
source "$F2T_CONFIG"
IP_BAM="${1:?sample BAM required}"; CTRL_BAM="${2:-none}"; OUT_DIR="${3:?output directory required}"
MODE="${4:-narrow}"; GENOME="${5:?genome required}"; SAMPLE="${6:-$(basename "$IP_BAM" .bam)}"; LAYOUT="${7:?PE or SE required}"
MACS3="${MACS3_COMMAND:-macs3}"
case "$GENOME" in hg38) GSIZE="${MACS3_GENOME_HG38:-${MACS2_GENOME_HG38:-hs}}" ;; mm39) GSIZE="${MACS3_GENOME_MM39:-${MACS2_GENOME_MM39:-mm}}" ;; *) echo "ERROR: unsupported genome: $GENOME" >&2; exit 1 ;; esac
[[ -s "$IP_BAM" ]] || { echo "ERROR: sample BAM missing: $IP_BAM" >&2; exit 1; }
mkdir -p "${OUT_DIR}/narrow" "${OUT_DIR}/broad"
control_args=(); [[ "$CTRL_BAM" != "none" && -s "$CTRL_BAM" ]] && control_args=( -c "$CTRL_BAM" )
format_args=()
if [[ "$LAYOUT" == "PE" ]]; then
    format_args=( -f BAMPE )
else
    format_args=( -f BAM --nomodel --shift "${MACS3_ATAC_SE_SHIFT:--75}" --extsize "${MACS3_ATAC_SE_EXTSIZE:-150}" )
fi

call_one() {
    local peak_type="$1" subdir="$2"; shift 2
    local expected="${OUT_DIR}/${subdir}/${SAMPLE}_peaks.${peak_type}"
    "$MACS3" callpeak -t "$IP_BAM" "${control_args[@]}" "${format_args[@]}" \
        -g "$GSIZE" -n "$SAMPLE" --outdir "${OUT_DIR}/${subdir}" \
        -q "${MACS3_QVALUE:-${MACS2_QVALUE:-0.05}}" --keep-dup all "$@" \
        >"${OUT_DIR}/${subdir}/${SAMPLE}_macs3.log" 2>&1
    if [[ ! -s "$expected" ]]; then
        if [[ "${ALLOW_EMPTY_PEAKS:-false}" == "true" ]]; then
            echo "WARNING: no ${peak_type} peaks for $SAMPLE" >&2; : > "$expected"
        else
            echo "ERROR: MACS3 produced no ${peak_type} peaks for $SAMPLE" >&2; exit 1
        fi
    fi
}

case "${MODE,,}" in
    none) echo "MACS3 skipped: $SAMPLE" ;;
    narrow) call_one narrowPeak narrow ;;
    broad) call_one broadPeak broad --broad --broad-cutoff "${MACS3_BROAD_CUTOFF:-${MACS2_BROAD_CUTOFF:-0.1}}" ;;
    both)
        call_one narrowPeak narrow
        call_one broadPeak broad --broad --broad-cutoff "${MACS3_BROAD_CUTOFF:-${MACS2_BROAD_CUTOFF:-0.1}}"
        ;;
    *) echo "ERROR: invalid peak mode: $MODE" >&2; exit 1 ;;
esac
