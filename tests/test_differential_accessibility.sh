#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIFFBIND_R="${REPO_DIR}/scripts/diffbind_analysis.R"
DESEQ_R="${REPO_DIR}/scripts/deseq2atac_analysis.R"
DESEQ_SH="${REPO_DIR}/scripts/deseq2atac_analysis.sh"
ENTRYPOINT="${REPO_DIR}/atacseq2tracks.sh"
DA_DOC="${REPO_DIR}/docs/13_differential_accessibility.md"

grep -q 'DIFFBIND_SUMMITS=100' "${REPO_DIR}/config/config.conf" \
    || { echo 'FAIL DiffBind summit default is not 100' >&2; exit 1; }
grep -q 'summits = summits' "$DIFFBIND_R" \
    || { echo 'FAIL configured DiffBind summit width is not passed to dba.count' >&2; exit 1; }
grep -q 'if (length(args) == 3).*as.integer' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind does not retain a configurable summit width' >&2; exit 1; }
echo 'OK   DiffBind defaults to configurable summits=100'

grep -q 'singleEnd = FALSE, fragments = FALSE' "$DESEQ_R" \
    || { echo 'FAIL paired-end strict fragment counting is missing' >&2; exit 1; }
grep -q 'singleEnd = TRUE, fragments = FALSE' "$DESEQ_R" \
    || { echo 'FAIL single-end read counting is missing' >&2; exit 1; }
grep -q 'isProperPair = TRUE' "$DESEQ_R" \
    || { echo 'FAIL paired-end proper-pair filter is missing' >&2; exit 1; }
grep -q 'BAM columns do not match samplesheet' "$DESEQ_R" \
    || { echo 'FAIL sample/count column-order validation is missing' >&2; exit 1; }
echo 'OK   paired-end fragments and single-end reads have distinct count semantics'

grep -q 'canonicalize_ranges' "$DESEQ_R" \
    || { echo 'FAIL canonical chromosome filtering is missing' >&2; exit 1; }
grep -q '!overlapsAny(peaks, blacklist' "$DESEQ_R" \
    || { echo 'FAIL blacklist filtering is missing from consensus construction' >&2; exit 1; }
grep -q 'atom_support >= min_support' "$DESEQ_R" \
    || { echo 'FAIL configurable peak-support threshold is missing' >&2; exit 1; }
grep -q 'any(consensus_support < min_support)' "$DESEQ_R" \
    || { echo 'FAIL consensus support is not validated' >&2; exit 1; }
echo 'OK   canonical, blacklist and minimum-support consensus contract'

grep -q 'peak_type.*broad.*narrow' "$DESEQ_R" \
    || { echo 'FAIL DESeq2ATAC R module does not validate broad/narrow peak type' >&2; exit 1; }
grep -q 'for PEAK_TYPE in broad narrow' "$DESEQ_SH" \
    || { echo 'FAIL DESeq2ATAC wrapper does not run both peak types' >&2; exit 1; }
grep -q 'TYPE_OUT="$OUT_DIR/$PEAK_TYPE"' "$DESEQ_SH" \
    || { echo 'FAIL DESeq2ATAC peak types do not have separate output directories' >&2; exit 1; }
grep -q 'deseq2atac_peak_type_summary.tsv' "$DESEQ_SH" \
    || { echo 'FAIL DESeq2ATAC combined peak-type summary is missing' >&2; exit 1; }
grep -q 'requires macs2_mode=both' "${REPO_DIR}/scripts/smoke_test.sh" \
    || { echo 'FAIL preflight does not require both DESeq2ATAC peak types' >&2; exit 1; }
grep -q 'Exact DESeq2ATAC consensus construction' "$DA_DOC" \
    || { echo 'FAIL exact DESeq2ATAC consensus documentation is missing' >&2; exit 1; }
grep -q 'support-filtered consensus' "$DA_DOC" \
    || { echo 'FAIL DESeq2ATAC support semantics are not documented' >&2; exit 1; }
grep -q 'not recentered to a fixed summit window' "$DA_DOC" \
    || { echo 'FAIL DESeq2ATAC versus DiffBind boundary behavior is not documented' >&2; exit 1; }
echo 'OK   DESeq2ATAC runs independent broad and narrow consensus analyses'

for output in \
    deseq2atac_consensus_peaks.bed \
    deseq2atac_raw_counts.tsv.gz \
    deseq2atac_normalized_counts.tsv.gz \
    deseq2atac_results_all.tsv.gz \
    deseq2atac_results_significant.tsv.gz \
    deseq2atac_summary.txt; do
    grep -q "$output" "$DESEQ_R" \
        || { echo "FAIL expected DESeq2ATAC output is missing from implementation: $output" >&2; exit 1; }
done
for figure in library_sizes_and_size_factors sample_correlation sample_distance \
    pca dispersion_estimates ma volcano; do
    grep -q "$figure" "$DESEQ_R" \
        || { echo "FAIL expected DESeq2ATAC figure is missing: $figure" >&2; exit 1; }
done
grep -q 'gzip -t' "$DESEQ_SH" \
    || { echo 'FAIL compressed outputs are not validated by the wrapper' >&2; exit 1; }
grep -q 'No regions passed FDR' "$DESEQ_R" \
    || { echo 'FAIL explicit zero-significant result handling is missing' >&2; exit 1; }
grep -q 'contrast = c("condition", numerator_condition, reference_condition)' "$DESEQ_R" \
    || { echo 'FAIL explicit DESeq2 contrast direction is missing' >&2; exit 1; }
grep -q 'Positive log2 fold change means' "$DESEQ_R" \
    || { echo 'FAIL contrast direction is not recorded in the summary' >&2; exit 1; }
echo 'OK   DESeq2ATAC output and zero-significant contracts'

grep -q 'if is_done 12a' "$ENTRYPOINT" \
    || { echo 'FAIL DESeq2ATAC checkpoint guard is missing' >&2; exit 1; }
grep -q 'mark_done 12a' "$ENTRYPOINT" \
    || { echo 'FAIL DESeq2ATAC checkpoint completion is missing' >&2; exit 1; }
grep -q 'DESeq2ATAC will still be attempted' "$ENTRYPOINT" \
    || { echo 'FAIL peer modules are not failure-isolated' >&2; exit 1; }
grep -q 'differential_accessibility_summary.tsv' "${REPO_DIR}/scripts/generate_pipeline_report.sh" \
    || { echo 'FAIL combined differential summary is missing from report stage' >&2; exit 1; }
grep -q 'for peak_type in broad narrow' "${REPO_DIR}/scripts/generate_pipeline_report.sh" \
    || { echo 'FAIL report does not include both DESeq2ATAC peak types' >&2; exit 1; }
echo 'OK   independent checkpoint, failure isolation and report integration'

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$TMP_DIR/bin" "$TMP_DIR/bams" "$TMP_DIR/peaks" "$TMP_DIR/out"
printf 'genome\nhg38\n' > "$TMP_DIR/samplesheet.csv"
printf 'chr1\t1\t2\n' > "$TMP_DIR/blacklist.bed"
cat > "$TMP_DIR/config.conf" <<EOF
BLACKLIST_HG38="$TMP_DIR/blacklist.bed"
R_BIN="$TMP_DIR/bin/Rscript"
MIN_MAPQ=30
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
mkdir -p "$out_dir/plots"
printf 'Status: SUCCESS\nPeak type: %s\nConsensus regions: 3\nNonzero tested regions: 2\nSignificant regions: 0\n' \
    "$peak_type" > "$out_dir/deseq2atac_summary.txt"
for file in deseq2atac_consensus_peaks.bed deseq2atac_sample_metadata.tsv \
    deseq2atac_library_summary.tsv deseq2atac_size_factors.tsv \
    deseq2atac_analysis_object.rds deseq2atac_session_info.txt; do
    printf 'mock\n' > "$out_dir/$file"
done
for file in deseq2atac_consensus_peaks_with_support.tsv.gz \
    deseq2atac_raw_counts.tsv.gz deseq2atac_normalized_counts.tsv.gz \
    deseq2atac_results_all.tsv.gz deseq2atac_results_significant.tsv.gz; do
    printf 'header\n' | gzip -c > "$out_dir/$file"
done
for stem in library_sizes_and_size_factors sample_correlation sample_distance \
    pca dispersion_estimates ma volcano; do
    printf 'mock\n' > "$out_dir/plots/${stem}.png"
    printf 'mock\n' > "$out_dir/plots/${stem}.pdf"
done
EOF
chmod +x "$TMP_DIR/bin/python3" "$TMP_DIR/bin/Rscript"
F2T_CONFIG="$TMP_DIR/config.conf" PATH="$TMP_DIR/bin:$PATH" \
    bash "$DESEQ_SH" "$TMP_DIR/samplesheet.csv" "$TMP_DIR/bams" \
        "$TMP_DIR/peaks" "$TMP_DIR/out"
[[ -s "$TMP_DIR/out/broad/deseq2atac_summary.txt" && \
   -s "$TMP_DIR/out/narrow/deseq2atac_summary.txt" ]] \
    || { echo 'FAIL mocked dual DESeq2ATAC outputs are missing' >&2; exit 1; }
[[ "$(wc -l < "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv")" -eq 3 ]] \
    || { echo 'FAIL DESeq2ATAC peak-type summary does not contain both analyses' >&2; exit 1; }
grep -q $'^broad\t' "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv" \
    || { echo 'FAIL broad DESeq2ATAC summary row is missing' >&2; exit 1; }
grep -q $'^narrow\t' "$TMP_DIR/out/deseq2atac_peak_type_summary.tsv" \
    || { echo 'FAIL narrow DESeq2ATAC summary row is missing' >&2; exit 1; }
echo 'OK   mocked wrapper executes and validates both peak-type analyses'

if command -v Rscript >/dev/null 2>&1 && \
   Rscript -e 'p <- c("DESeq2","GenomicAlignments","GenomicRanges","Rsamtools","rtracklayer","ggplot2","BiocParallel"); quit(status=ifelse(all(vapply(p, requireNamespace, logical(1), quietly=TRUE)),0,1))' >/dev/null 2>&1; then
    Rscript "$DESEQ_R" --self-test
else
    echo 'WARN R/Bioconductor packages unavailable; DESeq2ATAC execution self-test skipped'
fi
