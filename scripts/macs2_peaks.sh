#!/bin/bash
# ATACseq2tracks v3.0.4 — MACS2/MACS3 peak calling — always runs BOTH narrow AND broad
# Robust mode: retries with --nomodel on failure; creates empty peak file rather than aborting.
# Usage: bash scripts/macs2_peaks.sh <ip.bam> <ctrl.bam|none> <outDir> <mode> <genome_key> [sample_name]
set -uo pipefail   # NOTE: -e intentionally removed so we handle errors manually

_load_config() {
    if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
        source "${F2T_CONFIG}"
    else
        local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
        local _c="${_d}/../config/config.conf"
        [[ -f "$_c" ]] && source "$_c" || {
            echo "ERROR: config.conf not found. Export F2T_CONFIG or pass --config to atacseq2tracks.sh." >&2
            exit 1
        }
    fi
}
_load_config

IP_BAM="$1"; CTRL_BAM="$2"; OUT_DIR="$3"
MODE="${4:-both}"; GENOME_KEY="${5:-hg38}"; SAMPLE="${6:-$(basename "$IP_BAM" .bam)}"
mkdir -p "${OUT_DIR}/narrow" "${OUT_DIR}/broad"
[[ "$GENOME_KEY" == "hg38" ]] && GSIZE="$MACS2_GENOME_HG38" || GSIZE="$MACS2_GENOME_MM39"
CTRL_FLAG=""; [[ "$CTRL_BAM" != "none" && -f "$CTRL_BAM" ]] && CTRL_FLAG="-c $CTRL_BAM"

# Shared no-peaks registry (appended to by all parallel jobs; created by batch script)
NO_PEAKS_FILE="${OUT_DIR}/../no_peaks_samples.txt"

# ---------------------------------------------------------------------------
# _run_macs3_attempt: single MACS3 callpeak attempt; returns 0 if peak file
# produced with >0 lines, 1 otherwise.
# Usage: _run_macs3_attempt <peak_type> <subdir> <extra_flags...>
# ---------------------------------------------------------------------------
_run_macs3_attempt() {
    local peak_type="$1" sub="$2"; shift 2
    local extra_flags="$*"
    local peak_ext; [[ "$peak_type" == "broadPeak" ]] && peak_ext="broadPeak" || peak_ext="narrowPeak"
    local expected_peak="${OUT_DIR}/${sub}/${SAMPLE}_peaks.${peak_ext}"
    local log_file="${OUT_DIR}/${sub}/${SAMPLE}_macs2_${peak_type}.log"

    macs3 callpeak \
        -t "$IP_BAM" $CTRL_FLAG \
        -f BAM -g "$GSIZE" -n "$SAMPLE" \
        --outdir "${OUT_DIR}/${sub}" \
        -q "${MACS2_QVALUE}" $extra_flags \
        2>"$log_file"

    # Return success only if peak file exists and is non-empty
    if [[ -s "$expected_peak" ]]; then
        return 0
    else
        return 1
    fi
}

# ---------------------------------------------------------------------------
# run_macs2: robust wrapper with --nomodel retry and graceful empty fallback
# ---------------------------------------------------------------------------
run_macs2() {
    local peak_type="$1" sub="$2" extra_flags="$3"
    local peak_ext; [[ "$peak_type" == "broadPeak" ]] && peak_ext="broadPeak" || peak_ext="narrowPeak"
    local expected_peak="${OUT_DIR}/${sub}/${SAMPLE}_peaks.${peak_ext}"

    echo "[MACS2] $SAMPLE  mode=${peak_type}  attempt=1 (model-based)"

    # --- Attempt 1: standard model-based ---
    if _run_macs3_attempt "$peak_type" "$sub" "$extra_flags"; then
        echo "[MACS2] $SAMPLE  mode=${peak_type}  result=OK (model-based)"
        return 0
    fi

    echo "[WARN]  $SAMPLE  mode=${peak_type}  model-based produced no peaks -- retrying with --nomodel"

    # --- Attempt 2: no-model, relaxed cutoff, keep-dup all ---
    local nomodel_flags="--nomodel --extsize 200 --keep-dup all"
    [[ "$peak_type" == "broadPeak" ]] && \
        nomodel_flags="$nomodel_flags" || \
        nomodel_flags="$nomodel_flags"

    # Use relaxed q-value for no-model retry
    local orig_q="${MACS2_QVALUE}"; MACS2_QVALUE=0.1
    if _run_macs3_attempt "$peak_type" "$sub" "$extra_flags $nomodel_flags"; then
        MACS2_QVALUE="$orig_q"
        echo "[WARN]  $SAMPLE  mode=${peak_type}  result=OK (--nomodel fallback, q=0.1)"
        return 0
    fi
    MACS2_QVALUE="$orig_q"

    # --- Both attempts failed: create empty file and register sample ---
    echo "[WARN]  $SAMPLE  mode=${peak_type}  result=NO_PEAKS -- creating empty peak file"
    touch "$expected_peak"
    # Thread-safe append to shared no-peaks registry
    echo "${SAMPLE}" >> "$NO_PEAKS_FILE" 2>/dev/null || true
    return 0   # Always exit 0 so pipeline continues
}

if [[ "${MODE,,}" == "none" ]]; then
    echo "MACS2 skipped for $SAMPLE (mode=none)"; exit 0
fi

[[ "${MODE,,}" == "both" || "${MODE,,}" == "narrow" ]] && \
    run_macs2 "narrowPeak" "narrow" ""
[[ "${MODE,,}" == "both" || "${MODE,,}" == "broad"  ]] && \
    run_macs2 "broadPeak"  "broad"  "--broad --broad-cutoff ${MACS2_BROAD_CUTOFF}"

exit 0
