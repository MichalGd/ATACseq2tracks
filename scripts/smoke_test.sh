#!/usr/bin/env bash
# ATACseq2tracks v3.2.0 - pre-flight validation
set -euo pipefail

SAMPLESHEET_ARG="${1:?ERROR: pass samplesheet.csv as argument 1}"
CONFIG="${2:?ERROR: pass config.conf as argument 2}"
PASS=0
FAIL=0
WARN=0

ok()   { printf '[OK]   %s\n' "$1"; PASS=$((PASS + 1)); }
fail() { printf '[FAIL] %s\n' "$1" >&2; FAIL=$((FAIL + 1)); }
warn() { printf '[WARN] %s\n' "$1"; WARN=$((WARN + 1)); }

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"
SAMPLESHEET="$SAMPLESHEET_ARG"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION_FILE="${SCRIPT_DIR}/../VERSION"
VERSION="$(tr -d '[:space:]' < "$VERSION_FILE" 2>/dev/null || echo unknown)"

echo "ATACseq2tracks ${VERSION} pre-flight"

for script in validate_samplesheet.py fastqc_batch.sh trimgalore_batch.sh \
    bowtie2_batch.sh picard_dedup_batch.sh blacklist_filter.sh \
    blacklist_filter_batch.sh genomecoverage_single.sh genomecoverage_batch.sh merge_replicates.sh \
    macs2_peaks.sh macs2_batch.sh post_alignment_qc_batch.sh \
    consensus_peak_size_factors.R create_ucsc_tracks.sh ataqv_qc_batch.sh \
    prepare_tss_bed.py extract_ataqv_metrics.py plot_fragment_periodicity.py \
    plot_chrom_coverage.py peak_interpretation.sh prepare_diffbind.R \
    diffbind_analysis.sh diffbind_analysis.R generate_pipeline_report.sh; do
    [[ -f "${SCRIPT_DIR}/${script}" ]] && ok "script: ${script}" || fail "missing script: ${script}"
done

if [[ "${RUN_ATAQV_QC:-true}" == "true" ]]; then
    command -v ataqv >/dev/null 2>&1 && ok "tool: ataqv" || fail "tool not in PATH: ataqv"
    if [[ "${GENERATE_ATAQV_VIEWER:-true}" == "true" ]]; then
        command -v mkarv >/dev/null 2>&1 && ok "tool: mkarv" || fail "tool not in PATH: mkarv"
    fi
    python3 -c 'import matplotlib' >/dev/null 2>&1 \
        && ok "Python package: matplotlib" || fail "Python package not installed: matplotlib"
fi

if [[ "${RUN_PEAK_ANNOTATION:-false}" == "true" ]]; then
    command -v annotatePeaks.pl >/dev/null 2>&1 \
        && ok "tool: annotatePeaks.pl" || fail "tool not in PATH: annotatePeaks.pl"
fi
if [[ "${RUN_MOTIF_ENRICHMENT:-false}" == "true" ]]; then
    command -v findMotifsGenome.pl >/dev/null 2>&1 \
        && ok "tool: findMotifsGenome.pl" || fail "tool not in PATH: findMotifsGenome.pl"
fi

[[ -f "$SAMPLESHEET" ]] || fail "samplesheet not found: $SAMPLESHEET"
[[ -n "${OUTPUT_DIR:-}" ]] || fail "OUTPUT_DIR is not configured"

if [[ -f "$SAMPLESHEET" ]]; then
    python3 "${SCRIPT_DIR}/validate_samplesheet.py" "$SAMPLESHEET" && ok "samplesheet schema" || fail "samplesheet validation"

    mapfile -t USED_GENOMES < <(python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    print("\n".join(sorted({row["genome"].strip() for row in csv.DictReader(handle) if row["genome"].strip()})))
PY
    )
    if (( ${#USED_GENOMES[@]} > 1 )); then
        fail "multiple genome builds in one run (${USED_GENOMES[*]}); v3.2.0 requires one build per run"
    elif (( ${#USED_GENOMES[@]} == 1 )); then
        ok "single genome build: ${USED_GENOMES[0]}"
    fi

    BIOLOGICAL_SAMPLE_COUNT="$(python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    rows = csv.DictReader(handle)
    keys = {(r["sample_id"].strip(), r["replicate"].strip()) for r in rows
            if r["is_control"].strip().lower() not in {"true", "1", "yes"}}
print(len(keys))
PY
    )"
    if [[ "${GENERATE_DESEQ2_CONSENSUS_TRACKS:-true}" == "true" ]] && \
       (( BIOLOGICAL_SAMPLE_COUNT < ${CONSENSUS_MIN_SAMPLES:-2} )); then
        fail "DESeq2 consensus tracks require at least ${CONSENSUS_MIN_SAMPLES:-2} biological samples; found $BIOLOGICAL_SAMPLE_COUNT"
    else
        ok "biological samples: $BIOLOGICAL_SAMPLE_COUNT"
    fi
fi

MACS3_COMMAND="${MACS3_COMMAND:-macs3}"
for tool in bowtie2 samtools bedtools trim_galore fastqc "$MACS3_COMMAND" \
    multiqc python3 "${R_BIN:-Rscript}" bamCoverage multiBamSummary; do
    command -v "$tool" >/dev/null 2>&1 && ok "tool: $tool" || fail "tool not in PATH: $tool"
done

if [[ -n "${PICARD_JAR:-}" && -f "$PICARD_JAR" ]]; then
    ok "Picard JAR: $PICARD_JAR"
else
    fail "PICARD_JAR is missing or invalid: ${PICARD_JAR:-unset}"
fi

if [[ -n "${BEDGRAPH_TO_BIGWIG:-}" && -x "$BEDGRAPH_TO_BIGWIG" ]]; then
    ok "bedGraphToBigWig: $BEDGRAPH_TO_BIGWIG"
else
    warn "BEDGRAPH_TO_BIGWIG is unavailable; deepTools bigWigs still work, but legacy conversion does not"
fi

check_reference() {
    local label="$1" path="$2"
    [[ -f "$path" ]] && ok "$label: $path" || fail "$label not found: $path"
}

for genome in "${USED_GENOMES[@]:-}"; do
    case "$genome" in
        hg38)
            check_reference CHROM_SIZES_HUMAN "${CHROM_SIZES_HUMAN:-}"
            check_reference BLACKLIST_HG38 "${BLACKLIST_HG38:-}"
            if [[ "${RUN_ATAQV_QC:-true}" == "true" ]]; then
                if [[ -n "${TSS_BED_HG38:-}" ]]; then
                    check_reference TSS_BED_HG38 "$TSS_BED_HG38"
                else
                    check_reference GTF_HUMAN "${GTF_HUMAN:-}"
                fi
            fi
            index="${INDEX_HG38:-}"
            ;;
        mm39)
            check_reference CHROM_SIZES_MOUSE "${CHROM_SIZES_MOUSE:-}"
            check_reference BLACKLIST_MM39 "${BLACKLIST_MM39:-}"
            if [[ "${RUN_ATAQV_QC:-true}" == "true" ]]; then
                if [[ -n "${TSS_BED_MM39:-}" ]]; then
                    check_reference TSS_BED_MM39 "$TSS_BED_MM39"
                else
                    check_reference GTF_MOUSE "${GTF_MOUSE:-}"
                fi
            fi
            index="${INDEX_MM39:-}"
            ;;
        "") continue ;;
        *) fail "unsupported genome: $genome"; continue ;;
    esac
    if compgen -G "${index}*.bt2" >/dev/null || compgen -G "${index}*.bt2l" >/dev/null; then
        ok "Bowtie2 index: $index"
    else
        fail "Bowtie2 index not found: $index"
    fi
done

for package in DESeq2 DiffBind ggplot2; do
    if "${R_BIN:-Rscript}" -e "quit(status=ifelse(requireNamespace('${package}', quietly=TRUE), 0, 1))" >/dev/null 2>&1; then
        ok "R package: $package"
    else
        fail "R package not installed: $package"
    fi
done

echo "RESULT: ${PASS} OK | ${WARN} warnings | ${FAIL} failures"
(( FAIL == 0 )) || exit 1
