#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Pre-flight smoke test
# Usage: bash scripts/smoke_test.sh config/samplesheet.csv config/config.sh
# Must be run from the fastq2tracks/ root directory.
# =============================================================================
set -euo pipefail

SAMPLESHEET="${1:-config/samplesheet.csv}"
CONFIG="${2:-config/config.sh}"

PASS=0; FAIL=0; WARN=0
ok()   { echo "[OK]   $1"; ((PASS++));  }
fail() { echo "[FAIL] $1" >&2; ((FAIL++)); }
warn() { echo "[WARN] $1"; ((WARN++)); }

echo "============================================================"
echo " fastq2tracks v3.0 — Pre-flight check"
echo " $(date)"
echo "============================================================"

# --- 1. Config file ---
echo ""
echo "--- Config ---"
[[ -f "$CONFIG" ]] && ok "Config file exists: $CONFIG" || { fail "Config not found: $CONFIG"; }
source "$CONFIG" 2>/dev/null || fail "Config failed to source"

# --- 2. Workflow scripts ---
echo ""
echo "--- Scripts ---"
REQUIRED_SCRIPTS=(
    scripts/validate_samplesheet.py
    scripts/trimgalore_batch.sh
    scripts/fastqc_batch.sh
    scripts/bowtie2_align.sh
    scripts/bowtie2_batch.sh
    scripts/picard_dedup.sh
    scripts/picard_dedup_batch.sh
    scripts/blacklist_filter.sh
    scripts/genomecoverage_single.sh
    scripts/genomecoverage_batch.sh
    scripts/merge_replicates.sh
    scripts/macs2_peaks.sh
    scripts/macs2_batch.sh
    scripts/run_chipqc.R
    scripts/prepare_diffbind.R
    scripts/create_ucsc_tracks.sh
    scripts/generate_pipeline_report.sh
)
for s in "${REQUIRED_SCRIPTS[@]}"; do
    if [[ -f "$s" ]]; then
        [[ -x "$s" ]] && ok "Found + executable: $s" || warn "Found but NOT executable: $s (run: chmod +x $s)"
    else
        fail "Missing script: $s"
    fi
done

# --- 3. Sample sheet ---
echo ""
echo "--- Sample sheet ---"
[[ -f "$SAMPLESHEET" ]] && ok "Samplesheet found: $SAMPLESHEET" || fail "Samplesheet not found: $SAMPLESHEET"
python3 scripts/validate_samplesheet.py "$SAMPLESHEET" && ok "Samplesheet valid" || fail "Samplesheet validation failed"

# --- 4. Required tools ---
echo ""
echo "--- Tools ---"
TOOLS=(bowtie2 samtools bedtools trim_galore fastqc macs2 python3 Rscript)
for t in "${TOOLS[@]}"; do
    command -v "$t" &>/dev/null && ok "$t in PATH" || fail "$t NOT found in PATH"
done

# picard
[[ -f "${PICARD_JAR}" ]] && ok "picard.jar: $PICARD_JAR" || fail "picard.jar not found: $PICARD_JAR"

# bedGraphToBigWig
[[ -x "${BEDGRAPH_TO_BIGWIG}" ]] && ok "bedGraphToBigWig: $BEDGRAPH_TO_BIGWIG" || fail "bedGraphToBigWig not found/executable"

# --- 5. Reference files ---
echo ""
echo "--- Reference files ---"
check_ref() { [[ -f "$2" ]] && ok "$1: $2" || warn "$1 not found (will be needed at runtime): $2"; }
check_ref "CHROM_SIZES_HUMAN" "$CHROM_SIZES_HUMAN"
check_ref "CHROM_SIZES_MOUSE" "$CHROM_SIZES_MOUSE"
check_ref "BLACKLIST_HG38"    "$BLACKLIST_HG38"
check_ref "BLACKLIST_MM39"    "$BLACKLIST_MM39"
check_ref "CHIPQC_ANNO_HG38"  "$CHIPQC_ANNOTATION_HG38"
check_ref "CHIPQC_ANNO_MM39"  "$CHIPQC_ANNOTATION_MM39"

# bowtie2 index (check at least one .bt2 file)
for IDX in "$INDEX_HG38" "$INDEX_MM39"; do
    if ls "${IDX}"*.bt2 &>/dev/null 2>&1 || ls "${IDX}"*.bt2l &>/dev/null 2>&1; then
        ok "Bowtie2 index: $IDX"
    else
        warn "Bowtie2 index not found or incomplete: $IDX"
    fi
done

# --- 6. R packages ---
echo ""
echo "--- R packages ---"
Rscript -e "
pkgs <- c('ChIPQC','DiffBind','BiocParallel','GenomicAlignments','rtracklayer','ggplot2','dplyr')
for(p in pkgs){
  if(requireNamespace(p, quietly=TRUE)) cat('[OK]   R package:', p, '\n')
  else cat('[FAIL] R package NOT installed:', p, '\n')
}" 2>/dev/null || warn "Could not check R packages (Rscript failed)"

# --- 7. Disk space ---
echo ""
echo "--- Disk space ---"
OUTDIR="${OUTPUT_DIR:-.}"
AVAIL_GB=$(df -BG "$OUTDIR" 2>/dev/null | awk 'NR==2{gsub("G",""); print $4}')
if [[ -n "$AVAIL_GB" ]]; then
    [[ "$AVAIL_GB" -ge 200 ]] && ok "Disk space: ${AVAIL_GB}G available" || warn "Low disk space: ${AVAIL_GB}G available (recommend >=200G)"
fi

# --- Summary ---
echo ""
echo "============================================================"
echo " RESULT: ${PASS} OK  |  ${WARN} WARNINGS  |  ${FAIL} FAILURES"
echo "============================================================"

[[ $FAIL -eq 0 ]] && echo "Pre-flight PASSED — safe to launch workflow." && exit 0
echo "Pre-flight FAILED — resolve errors above before running." >&2
exit 1
