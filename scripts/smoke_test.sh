#!/bin/bash
# =============================================================================
# ATACseq2tracks v3.0.3 — Pre-flight smoke test
# Usage: bash scripts/smoke_test.sh <samplesheet.csv> <config.conf>
# =============================================================================
set -euo pipefail

SAMPLESHEET="${1:?ERROR: pass samplesheet.csv as arg 1}"
CONFIG="${2:?ERROR: pass config.conf as arg 2}"

PASS=0; FAIL=0; WARN=0

ok()   { echo "[OK]   $1"; ((PASS++))  || true; }
fail() { echo "[FAIL] $1" >&2; ((FAIL++)) || true; }
warn() { echo "[WARN] $1"; ((WARN++)) || true; }

echo "============================================================"
echo " ATACseq2tracks v3.0.3 -- Pre-flight check"
echo " $(date)"
echo "============================================================"

echo ""; echo "--- Config ---"
if [[ -f "$CONFIG" ]]; then
    ok "Config exists: $CONFIG"
else
    fail "Config not found: $CONFIG"; exit 1
fi

if source "$CONFIG" 2>/dev/null; then
    ok "Config sourced OK"
else
    fail "Config failed to source: $CONFIG"; exit 1
fi

if [[ -n "${SAMPLESHEET:-}" ]]; then
    ok "SAMPLESHEET defined: $SAMPLESHEET"
else
    fail "SAMPLESHEET not set in $CONFIG"
fi

if [[ -n "${OUTPUT_DIR:-}" ]]; then
    ok "OUTPUT_DIR defined: $OUTPUT_DIR"
else
    fail "OUTPUT_DIR not set in $CONFIG"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""; echo "--- Scripts ---"
REQUIRED_SCRIPTS=(
    validate_samplesheet.py trimgalore_batch.sh   fastqc_batch.sh
    bowtie2_align.sh        bowtie2_batch.sh
    picard_dedup.sh         picard_dedup_batch.sh
    blacklist_filter.sh     blacklist_filter_batch.sh
    genomecoverage_single.sh genomecoverage_batch.sh
    merge_replicates.sh     macs2_peaks.sh macs2_batch.sh
    run_chipqc.R            prepare_diffbind.R
    create_ucsc_tracks.sh   generate_pipeline_report.sh
)
for s in "${REQUIRED_SCRIPTS[@]}"; do
    fp="${SCRIPT_DIR}/${s}"
    if [[ -f "$fp" ]]; then
        if [[ -x "$fp" ]]; then ok "Found + executable: $s"
        else warn "NOT executable: $s  (fix: chmod +x $fp)"; fi
    else
        fail "Missing script: $fp"
    fi
done

echo ""; echo "--- Samplesheet ---"
if [[ -f "$SAMPLESHEET" ]]; then
    ok "Samplesheet found: $SAMPLESHEET"
else
    fail "Samplesheet not found: $SAMPLESHEET"
fi

if python3 "${SCRIPT_DIR}/validate_samplesheet.py" "$SAMPLESHEET"; then
    ok "Samplesheet valid"
else
    fail "Samplesheet validation failed"
fi

echo ""; echo "--- Tools ---"
for t in bowtie2 samtools bedtools trim_galore fastqc macs2 multiqc python3 Rscript; do
    if command -v "$t" &>/dev/null; then ok "$t in PATH"
    else fail "$t NOT in PATH"; fi
done

if [[ -f "${PICARD_JAR}" ]]; then ok "picard.jar: $PICARD_JAR"
else fail "picard.jar not found: $PICARD_JAR"; fi

if [[ -x "${BEDGRAPH_TO_BIGWIG}" ]]; then ok "bedGraphToBigWig: $BEDGRAPH_TO_BIGWIG"
else fail "bedGraphToBigWig not executable: $BEDGRAPH_TO_BIGWIG"; fi

echo ""; echo "--- Reference files ---"
check_ref() {
    if [[ -f "$2" ]]; then ok "$1: $2"
    else warn "$1 not found (needed at runtime): $2"; fi
}
check_ref "CHROM_SIZES_HUMAN"         "$CHROM_SIZES_HUMAN"
check_ref "CHROM_SIZES_MOUSE"         "$CHROM_SIZES_MOUSE"
check_ref "BLACKLIST_HG38"            "$BLACKLIST_HG38"
check_ref "BLACKLIST_MM39"            "$BLACKLIST_MM39"
check_ref "CHIPQC_ANNOTATION_HG38"    "$CHIPQC_ANNOTATION_HG38"
check_ref "CHIPQC_ANNOTATION_MM39"    "$CHIPQC_ANNOTATION_MM39"
check_ref "CHIPQC_BLACKLIST_HG38_RDS" "$CHIPQC_BLACKLIST_HG38_RDS"
check_ref "CHIPQC_BLACKLIST_MM39_RDS" "$CHIPQC_BLACKLIST_MM39_RDS"
for IDX in "$INDEX_HG38" "$INDEX_MM39"; do
    if ls "${IDX}"*.bt2 &>/dev/null 2>&1 || ls "${IDX}"*.bt2l &>/dev/null 2>&1; then
        ok "Bowtie2 index: $IDX"
    else
        warn "Bowtie2 index not found: $IDX"
    fi
done

echo ""; echo "--- R packages ---"
for pkg in ChIPQC DiffBind BiocParallel GenomicAlignments rtracklayer ggplot2 dplyr; do
    if Rscript -e "if (!requireNamespace('${pkg}', quietly=TRUE)) quit(status=1)" 2>/dev/null; then
        ok "R package: ${pkg}"
    else
        fail "R package NOT installed: ${pkg}"
    fi
done

echo ""; echo "--- Disk space ---"
AVAIL_GB=$(df -BG "${OUTPUT_DIR:-.}" 2>/dev/null | awk 'NR==2{gsub("G",""); print $4}')
if [[ -n "$AVAIL_GB" ]]; then
    if [[ "$AVAIL_GB" -ge 200 ]]; then ok "Disk: ${AVAIL_GB}G available"
    else warn "Low disk: ${AVAIL_GB}G (recommend >=200G)"; fi
fi

echo ""; echo "============================================================"
echo " RESULT: ${PASS} OK | ${WARN} WARNINGS | ${FAIL} FAILURES"
echo "============================================================"
if [[ $FAIL -eq 0 ]]; then
    echo "Pre-flight PASSED -- ready to run."
    exit 0
else
    echo "Pre-flight FAILED -- resolve errors above before running." >&2
    exit 1
fi
