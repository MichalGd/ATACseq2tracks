#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIFFBIND_R="${REPO_DIR}/scripts/diffbind_analysis.R"
DIFFBIND_SH="${REPO_DIR}/scripts/diffbind_analysis.sh"
REPORT_SH="${REPO_DIR}/scripts/generate_pipeline_report.sh"

grep -q 'prefilter_diffbind_peaks(ss, out_dir, genome, blacklist_file)' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind does not prefilter every peak set before initial counting' >&2; exit 1; }
grep -q 'sampleSheet = prefiltered_ss_file' "$DIFFBIND_R" \
    || { echo 'FAIL initial DiffBind count does not use the prefiltered sample sheet' >&2; exit 1; }
grep -q 'eligible_ss <- ss_prefiltered' "$DIFFBIND_R" \
    || { echo 'FAIL eligible-sample recount does not use prefiltered peak paths' >&2; exit 1; }
grep -q 'diffbind_peak_prefilter_manifest.tsv' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind peak prefiltering is not audited in a manifest' >&2; exit 1; }
grep -q 'canonical & !blacklisted' "$DIFFBIND_R" \
    || { echo 'FAIL DiffBind prefilter does not combine canonical and blacklist rules' >&2; exit 1; }
grep -q 'diffbind_peak_prefilter_manifest.tsv' "$DIFFBIND_SH" \
    || { echo 'FAIL DiffBind wrapper does not validate its prefilter manifest' >&2; exit 1; }
echo 'OK   DiffBind filters original peaks before initial consensus counting'

grep -q -- '--exclude deeptools' "$REPORT_SH" \
    || { echo 'FAIL incompatible native MultiQC deepTools parser is not excluded' >&2; exit 1; }
grep -q "ignore_images: false" "$REPORT_SH" \
    || { echo 'FAIL MultiQC custom images are not enabled' >&2; exit 1; }
grep -q 'Error converting colou?r' "$REPORT_SH" \
    || { echo 'FAIL MultiQC log audit does not reject colour-conversion errors' >&2; exit 1; }
echo 'OK   MultiQC uses validated static deepTools custom content'

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
MOCK_BIN="${TMP_DIR}/bin"
OUT_DIR="${TMP_DIR}/output"
REPORT_DIR="${OUT_DIR}/reports/pipeline_report_test"
mkdir -p "$MOCK_BIN" "${OUT_DIR}/qc_post_alignment/plots" "$REPORT_DIR"
printf '# test config\n' > "${TMP_DIR}/config.conf"

cat > "${MOCK_BIN}/mock_Rscript" <<'MOCKR'
#!/usr/bin/env bash
set -euo pipefail
output="$3"
mkdir -p "$output"
for required in diffbind_samplesheet_prefiltered.csv diffbind_consensus_peaks.bed \
    differential_accessibility_condition_eligibility.tsv \
    differential_accessibility_comparisons.tsv diffbind_summary.txt diffbind_log.txt; do
    printf 'mock\n' > "${output}/${required}"
done
if [[ "${MOCK_DIFFBIND_OMIT_MANIFEST:-0}" != 1 ]]; then
    printf 'sample_id\tinput_regions\tretained_regions\nmock\t1\t1\n' \
        > "${output}/diffbind_peak_prefilter_manifest.tsv"
fi
MOCKR
chmod +x "${MOCK_BIN}/mock_Rscript"

DIFFBIND_INPUT="${TMP_DIR}/diffbind_input"
DIFFBIND_OUTPUT="${TMP_DIR}/diffbind_output"
mkdir -p "$DIFFBIND_INPUT"
printf 'SampleID,Condition,Peaks\nA,one,/mock/a.bed\n' \
    > "${DIFFBIND_INPUT}/diffbind_samplesheet_hg38_broad.csv"
PATH="${MOCK_BIN}:$PATH" R_BIN=mock_Rscript RUN_SIMPLE_PEAK_ANNOTATION=false \
    F2T_CONFIG="${TMP_DIR}/config.conf" \
    bash "$DIFFBIND_SH" "$DIFFBIND_INPUT" "$DIFFBIND_OUTPUT" >/dev/null
[[ -s "${DIFFBIND_OUTPUT}/diffbind_samplesheet_hg38_broad/diffbind_peak_prefilter_manifest.tsv" ]] \
    || { echo 'FAIL mocked DiffBind wrapper did not validate the prefilter manifest' >&2; exit 1; }

if PATH="${MOCK_BIN}:$PATH" R_BIN=mock_Rscript RUN_SIMPLE_PEAK_ANNOTATION=false \
        MOCK_DIFFBIND_OMIT_MANIFEST=1 F2T_CONFIG="${TMP_DIR}/config.conf" \
        bash "$DIFFBIND_SH" "$DIFFBIND_INPUT" "${TMP_DIR}/diffbind_invalid" >/dev/null 2>&1; then
    echo 'FAIL DiffBind wrapper accepted a missing peak-prefilter manifest' >&2
    exit 1
fi
echo 'OK   DiffBind wrapper rejects incomplete prefilter outputs'

for image in pca_bins pca_peaks correlation_heatmap_pearson \
    correlation_heatmap_spearman correlation_heatmap_pearson_peaks \
    fingerprint heatmap_signal_over_peaks profile_signal_over_peaks; do
    printf 'PNG-%s\n' "$image" > "${OUT_DIR}/qc_post_alignment/plots/${image}.png"
done

cat > "${MOCK_BIN}/python3" <<'MOCKPY'
#!/usr/bin/env bash
set -euo pipefail
summary_dir="$3"
mkdir -p "$summary_dir"
printf 'module\tstatus\nmock\tSUCCESS\n' > "${summary_dir}/differential_accessibility_summary.tsv"
printf '<html><body>mock differential summary</body></html>\n' > "${summary_dir}/differential_accessibility_summary.html"
MOCKPY
chmod +x "${MOCK_BIN}/python3"

cat > "${MOCK_BIN}/multiqc" <<'MOCKMQC'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "$MOCK_MULTIQC_ARGS"
custom_dir="$2"
find "$custom_dir" -maxdepth 1 -type f -name '*_mqc.png' -print | sort > "$MOCK_CUSTOM_IMAGES"
out_dir=""
report_name=""
while (($#)); do
    case "$1" in
        -o) out_dir="$2"; shift 2 ;;
        -n) report_name="$2"; shift 2 ;;
        *) shift ;;
    esac
done
mkdir -p "$out_dir"
printf '<html><body>mock MultiQC</body></html>\n' > "${out_dir}/${report_name}.html"
if [[ "${MOCK_MULTIQC_FAILURE:-0}" == 1 ]]; then
    echo "Oops! The 'deeptools' MultiQC module broke..."
else
    echo 'MultiQC complete'
fi
MOCKMQC
chmod +x "${MOCK_BIN}/multiqc"

export PATH="${MOCK_BIN}:$PATH"
export F2T_CONFIG="${TMP_DIR}/config.conf"
export MOCK_MULTIQC_ARGS="${TMP_DIR}/multiqc_args.txt"
export MOCK_CUSTOM_IMAGES="${TMP_DIR}/custom_images.txt"

bash "$REPORT_SH" "$OUT_DIR" "$REPORT_DIR" html >/dev/null
[[ -s "${REPORT_DIR}/fastq2tracks_unified_$(date +%Y%m%d).html" ]] \
    || { echo 'FAIL mocked MultiQC report was not validated' >&2; exit 1; }
[[ "$(wc -l < "$MOCK_CUSTOM_IMAGES")" -eq 8 ]] \
    || { echo 'FAIL expected eight staged deepTools custom-content images' >&2; exit 1; }
grep -Fxq -- '--exclude' "$MOCK_MULTIQC_ARGS" \
    || { echo 'FAIL MultiQC exclusion flag was not passed' >&2; exit 1; }
grep -Fxq 'deeptools' "$MOCK_MULTIQC_ARGS" \
    || { echo 'FAIL native deepTools module was not excluded' >&2; exit 1; }
echo 'OK   report wrapper embeds all available deepTools plots and validates output'

if MOCK_MULTIQC_FAILURE=1 bash "$REPORT_SH" "$OUT_DIR" "$REPORT_DIR" html >/dev/null 2>&1; then
    echo 'FAIL report wrapper accepted a caught MultiQC module failure' >&2
    exit 1
fi
echo 'OK   caught MultiQC module failures propagate to the workflow'
