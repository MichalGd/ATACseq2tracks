#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIFFBIND_R="${REPO_DIR}/scripts/diffbind_analysis.R"
DIFFBIND_SH="${REPO_DIR}/scripts/diffbind_analysis.sh"
DESEQ_R="${REPO_DIR}/scripts/deseq2atac_analysis.R"
DESEQ_SH="${REPO_DIR}/scripts/deseq2atac_analysis.sh"
ANNOTATION_R="${REPO_DIR}/scripts/peak_annotation_helpers.R"
REPORT_PY="${REPO_DIR}/scripts/summarize_differential_accessibility.py"
ENTRYPOINT="${REPO_DIR}/atacseq2tracks.sh"
DA_DOC="${REPO_DIR}/docs/13_differential_accessibility.md"
REAL_PYTHON="$(command -v "${PYTHON_BIN:-python3}")"

grep -q 'DIFFBIND_SUMMITS=100' "${REPO_DIR}/config/config.conf" \
    || { echo 'FAIL DiffBind summit default is not 100' >&2; exit 1; }
grep -q '^RUN_CCRE_ANNOTATION=true$' "${REPO_DIR}/config/config.conf" \
    || { echo 'FAIL cCRE annotation is not enabled by default' >&2; exit 1; }
grep -q 'Set RUN_CCRE_ANNOTATION=false for GTF-only annotation' "$DESEQ_SH" \
    || { echo 'FAIL DESeq2ATAC lacks the explicit GTF-only opt-out' >&2; exit 1; }
grep -q 'Set RUN_CCRE_ANNOTATION=false for GTF-only annotation' "$DIFFBIND_SH" \
    || { echo 'FAIL DiffBind lacks the explicit GTF-only opt-out' >&2; exit 1; }
grep -q 'dba_all, minOverlap = 2, summits = summits' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind all-sample consensus counting is missing' >&2; exit 1; }
grep -q 'peaks = all_consensus, summits = FALSE, filter = 0' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind does not preserve the all-sample universe for its eligible model' >&2; exit 1; }
grep -q 'annotation_canonicalize(all_consensus, genome)' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind consensus is not restricted to canonical chromosomes' >&2; exit 1; }
grep -q 'overlapsAny(all_consensus, blacklist' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind consensus is not explicitly blacklist-filtered' >&2; exit 1; }
grep -q 'IRanges::overlapsAny(all_consensus, blacklist' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind does not use the exported IRanges overlap generic' >&2; exit 1; }
grep -q 'comparison_plan(eligible_order)' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind universal pair planning is missing' >&2; exit 1; }
grep -q 'other peak types remain eligible to run' "$DIFFBIND_SH" \
    || { echo 'FAIL DiffBind broad/narrow failure isolation is missing' >&2; exit 1; }
echo 'OK   DiffBind uses an all-sample consensus and universal eligible-condition pairs'

grep -q 'singleEnd = FALSE, fragments = FALSE' "$DESEQ_R" \
    || { echo 'FAIL paired-end strict fragment counting is missing' >&2; exit 1; }
grep -q 'singleEnd = TRUE, fragments = FALSE' "$DESEQ_R" \
    || { echo 'FAIL single-end read counting is missing' >&2; exit 1; }
grep -q 'isProperPair = TRUE' "$DESEQ_R" \
    || { echo 'FAIL paired-end proper-pair filter is missing' >&2; exit 1; }
grep -q 'model_metadata <- metadata\[metadata\$condition %in% eligible_order' "$DESEQ_R" \
    || { echo 'FAIL singleton conditions are not excluded only at model construction' >&2; exit 1; }
grep -q 'included_in_consensus = TRUE' "$DESEQ_R" \
    || { echo 'FAIL all-sample consensus eligibility contract is missing' >&2; exit 1; }
grep -q 'comparison_plan(eligible_order)' "$DESEQ_R" \
    || { echo 'FAIL DESeq2ATAC universal pair planning is missing' >&2; exit 1; }
grep -q 'pairwise contrasts' "$DESEQ_R" \
    || { echo 'FAIL DESeq2ATAC one-model/many-contrasts contract is missing' >&2; exit 1; }
echo 'OK   DESeq2ATAC counts PE fragments/SE reads and models all eligible pairs'

grep -q 'gene_context' "$ANNOTATION_R" \
    || { echo 'FAIL GTF gene-context annotation is missing' >&2; exit 1; }
grep -q 'nearest_tss_distance_bp' "$ANNOTATION_R" \
    || { echo 'FAIL nearest promoter/TSS annotation is missing' >&2; exit 1; }
grep -q 'ccre_primary_class' "$ANNOTATION_R" \
    || { echo 'FAIL cCRE classification is missing' >&2; exit 1; }
grep -q 'pELS|dELS' "$ANNOTATION_R" \
    || { echo 'FAIL enhancer-like cCRE classification is missing' >&2; exit 1; }
if grep -Eq 'GenomicRanges::(overlapsAny|findOverlaps)' "$DIFFBIND_R" "$DESEQ_R" "$ANNOTATION_R"; then
    echo 'FAIL non-exported overlap generic requested from GenomicRanges' >&2
    exit 1
fi
grep -q 'IRanges::overlapsAny' "$DESEQ_R" \
    || { echo 'FAIL DESeq2ATAC does not use IRanges::overlapsAny' >&2; exit 1; }
grep -q 'IRanges::overlapsAny' "$ANNOTATION_R" \
    || { echo 'FAIL annotation helper does not use IRanges::overlapsAny' >&2; exit 1; }
grep -q 'IRanges::findOverlaps' "$ANNOTATION_R" \
    || { echo 'FAIL annotation helper does not use IRanges::findOverlaps' >&2; exit 1; }
echo 'OK   built-in GTF and cCRE annotation fields are implemented'

grep -q 'differential_accessibility_comparisons.tsv' "$DESEQ_R" \
    || { echo 'FAIL DESeq2ATAC pair-level summary is missing' >&2; exit 1; }
grep -q 'differential_accessibility_comparisons.tsv' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind pair-level summary is missing' >&2; exit 1; }
grep -q 'differential_accessibility_summary.html' "$REPORT_PY" \
    || { echo 'FAIL combined HTML summary is missing' >&2; exit 1; }
grep -q 'reports were generated and automatic cleanup was suppressed' "$ENTRYPOINT" \
    || { echo 'FAIL failure reporting/cleanup suppression is missing' >&2; exit 1; }
grep -q '^ENABLE_AUTOMATIC_CLEANUP=true$' "${REPO_DIR}/config/config.conf" \
    || { echo 'FAIL automatic cleanup is not enabled by default' >&2; exit 1; }
grep -q 'ENABLE_AUTOMATIC_CLEANUP:-true' "$ENTRYPOINT" \
    || { echo 'FAIL legacy-config cleanup fallback is not enabled' >&2; exit 1; }
grep -q '^KEEP_FILTERED_BAMS=true$' "${REPO_DIR}/config/config.conf" \
    || { echo 'FAIL quantitative filtered BAMs are not retained by default' >&2; exit 1; }
grep -q 'KEEP_FILTERED_BAMS:-true' "${REPO_DIR}/scripts/cleanup_intermediates.sh" \
    || { echo 'FAIL legacy-config filtered-BAM retention fallback is unsafe' >&2; exit 1; }
echo 'OK   pair-level TSV/HTML reporting and failure-safe cleanup behavior'

grep -q 'all non-control biological samples' "$DA_DOC" \
    || { echo 'FAIL all-sample consensus policy is not documented' >&2; exit 1; }
grep -q 'pairwise' "$DA_DOC" \
    || { echo 'FAIL universal pairwise behavior is not documented' >&2; exit 1; }
grep -q 'ENCODE4' "$DA_DOC" \
    || { echo 'FAIL human cCRE provenance is not documented' >&2; exit 1; }
grep -q 'liftOver' "$DA_DOC" \
    || { echo 'FAIL mouse cCRE provenance is not documented' >&2; exit 1; }
echo 'OK   universal comparison and regulatory annotation policy is documented'

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$TMP_DIR/bin" "$TMP_DIR/bams" "$TMP_DIR/peaks" "$TMP_DIR/out"
printf 'genome\nhg38\n' > "$TMP_DIR/samplesheet.csv"
printf 'chr1\t1\t2\n' > "$TMP_DIR/blacklist.bed"
cat > "$TMP_DIR/config.conf" <<EOF
BLACKLIST_HG38="$TMP_DIR/blacklist.bed"
R_BIN="$TMP_DIR/bin/Rscript"
MIN_MAPQ=30
RUN_SIMPLE_PEAK_ANNOTATION=false
EOF
cat > "$TMP_DIR/bin/python3" <<'EOF'
#!/usr/bin/env bash
printf 'hg38\n'
EOF
cat > "$TMP_DIR/bin/Rscript" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
out_dir="$5"
peak_type="${13}"
comparison_dir="$out_dir/comparisons/001_B_vs_A"
mkdir -p "$out_dir/plots" "$comparison_dir/plots"
printf 'Status: SUCCESS\nPeak type: %s\nConsensus regions: 3\nPlanned pairwise comparisons: 1\nSuccessful comparisons: 1\nFailed comparisons: 0\n' \
    "$peak_type" > "$out_dir/deseq2atac_summary.txt"
for file in deseq2atac_consensus_peaks.bed deseq2atac_sample_metadata.tsv \
    deseq2atac_all_sample_metadata.tsv deseq2atac_library_summary.tsv \
    deseq2atac_size_factors.tsv deseq2atac_analysis_object.rds \
    deseq2atac_session_info.txt differential_accessibility_condition_eligibility.tsv; do
    printf 'mock\n' > "$out_dir/$file"
done
for file in deseq2atac_consensus_peaks_with_support.tsv.gz \
    deseq2atac_raw_counts.tsv.gz deseq2atac_normalized_counts.tsv.gz; do
    printf 'header\n' | gzip -c > "$out_dir/$file"
done
printf 'module\tpeak_type\tcomparison_id\tnumerator\treference\tnumerator_replicates\treference_replicates\tconsensus_regions\ttested_sites\tsignificant_sites\thigher_in_numerator\thigher_in_reference\talpha\tmin_abs_log2fc\tstatus\tresults_all\tresults_significant\tsummary_file\tmessage\n' > "$out_dir/differential_accessibility_comparisons.tsv"
printf 'DESeq2ATAC\t%s\t001_B_vs_A\tB\tA\t2\t2\t3\t2\t0\t0\t0\t0.05\t0\tSUCCESS\tall\tsig\tsummary\tcompleted_with_zero_significant_sites\n' \
    "$peak_type" >> "$out_dir/differential_accessibility_comparisons.tsv"
for stem in library_sizes_and_size_factors sample_correlation sample_distance \
    pca dispersion_estimates; do
    printf 'mock\n' > "$out_dir/plots/${stem}.png"
    printf 'mock\n' > "$out_dir/plots/${stem}.pdf"
done
EOF
chmod +x "$TMP_DIR/bin/python3" "$TMP_DIR/bin/Rscript"
F2T_CONFIG="$TMP_DIR/config.conf" PATH="$TMP_DIR/bin:$PATH" \
    bash "$DESEQ_SH" "$TMP_DIR/samplesheet.csv" "$TMP_DIR/bams" \
        "$TMP_DIR/peaks" "$TMP_DIR/out"
[[ "$(wc -l < "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv")" -eq 3 ]] \
    || { echo 'FAIL DESeq2ATAC peak-type summary does not contain both analyses' >&2; exit 1; }
grep -q $'^broad\tSUCCESS' "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv" \
    || { echo 'FAIL broad DESeq2ATAC summary row is missing' >&2; exit 1; }
grep -q $'^narrow\tSUCCESS' "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv" \
    || { echo 'FAIL narrow DESeq2ATAC summary row is missing' >&2; exit 1; }
echo 'OK   mocked wrapper executes and validates both peak-type analyses'

mkdir -p "$TMP_DIR/pipeline/deseq2atac/broad" "$TMP_DIR/reports"
cp "$TMP_DIR/out/broad/differential_accessibility_comparisons.tsv" \
    "$TMP_DIR/pipeline/deseq2atac/broad/"
"$REAL_PYTHON" "$REPORT_PY" "$TMP_DIR/pipeline" "$TMP_DIR/reports"
[[ "$(wc -l < "$TMP_DIR/reports/differential_accessibility_summary.tsv")" -eq 2 ]] \
    || { echo 'FAIL combined differential summary row count is wrong' >&2; exit 1; }
grep -q '001_B_vs_A' "$TMP_DIR/reports/differential_accessibility_summary.html" \
    || { echo 'FAIL combined differential HTML lacks comparison' >&2; exit 1; }
echo 'OK   combined pair-level TSV and HTML report generation'

if command -v Rscript >/dev/null 2>&1 && \
   Rscript -e 'p <- c("DESeq2","GenomicAlignments","GenomicRanges","IRanges","Rsamtools","rtracklayer","ggplot2","BiocParallel"); quit(status=ifelse(all(vapply(p, requireNamespace, logical(1), quietly=TRUE)),0,1))' >/dev/null 2>&1; then
    Rscript "${REPO_DIR}/tests/test_iranges_namespace.R"
    Rscript "$DESEQ_R" --self-test
else
    echo 'WARN R/Bioconductor packages unavailable; DESeq2ATAC execution self-test skipped'
fi
