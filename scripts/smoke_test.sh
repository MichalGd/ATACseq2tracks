#!/usr/bin/env bash
# ATACseq2tracks pre-flight validation
set -euo pipefail

SAMPLESHEET_ARG="${1:?ERROR: pass samplesheet.csv as argument 1}"
CONFIG="${2:?ERROR: pass config.conf as argument 2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_VERSION="$(tr -d '[:space:]' < "${SCRIPT_DIR}/../VERSION" 2>/dev/null || echo unknown)"
INPUT_SANITIZER="${SCRIPT_DIR}/sanitize_text_inputs.py"
PASS=0
FAIL=0
WARN=0

ok()   { printf '[OK]   %s\n' "$1"; PASS=$((PASS + 1)); }
fail() { printf '[FAIL] %s\n' "$1" >&2; FAIL=$((FAIL + 1)); }
warn() { printf '[WARN] %s\n' "$1"; WARN=$((WARN + 1)); }

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
[[ -f "$INPUT_SANITIZER" ]] || { echo "ERROR: input sanitizer not found: $INPUT_SANITIZER" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 is required to validate input text files" >&2; exit 1; }

# This must happen before source: CRLF/UTF BOM artifacts can otherwise make a
# valid Bash configuration fail before preflight begins.
python3 "$INPUT_SANITIZER" "$CONFIG"
# shellcheck disable=SC1090
source "$CONFIG"
if [[ "${TOTAL_CPU_BUDGET:-0}" =~ ^[0-9]+$ ]]; then
    ok "TOTAL_CPU_BUDGET: ${TOTAL_CPU_BUDGET:-0}"
else
    fail "TOTAL_CPU_BUDGET must be a non-negative integer"
fi
case "${RESOURCE_CHECK_MODE:-warn}" in
    warn|error|off) ok "RESOURCE_CHECK_MODE: ${RESOURCE_CHECK_MODE:-warn}" ;;
    *) fail "RESOURCE_CHECK_MODE must be warn, error or off" ;;
esac
for setting_value in \
    "QC_SAMPLE_PARALLEL_JOBS=${QC_SAMPLE_PARALLEL_JOBS:-4}" \
    "ATAQV_PARALLEL_JOBS=${ATAQV_PARALLEL_JOBS:-4}" \
    "TRACK_PARALLEL_JOBS=${TRACK_PARALLEL_JOBS:-2}" \
    "COVERAGE_FILTER_PARALLEL_JOBS=${COVERAGE_FILTER_PARALLEL_JOBS:-2}" \
    "SPIKEIN_PARALLEL_JOBS=${SPIKEIN_PARALLEL_JOBS:-2}" \
    "POOLED_MACS_PARALLEL_JOBS=${POOLED_MACS_PARALLEL_JOBS:-2}" \
    "MERGE_PARALLEL_JOBS=${MERGE_PARALLEL_JOBS:-2}"; do
    setting_name="${setting_value%%=*}"
    value="${setting_value#*=}"
    if [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
        ok "$setting_name: $value"
    else
        fail "$setting_name must be a positive integer"
    fi
done
case "${RUN_CCRE_ANNOTATION:-true}" in
    true|false) ;;
    *) fail "RUN_CCRE_ANNOTATION must be true or false" ;;
esac
case "${ENABLE_AUTOMATIC_CLEANUP:-true}" in
    true) ok "automatic cleanup enabled after full success" ;;
    false) ok "automatic cleanup disabled; all intermediates will be retained" ;;
    *) fail "ENABLE_AUTOMATIC_CLEANUP must be true or false" ;;
esac
for setting_name in GENERATE_CPM_TRACKS GENERATE_DESEQ2_CONSENSUS_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS \
    GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS \
    GENERATE_DROSOPHILA_CONTROL_TRACKS \
    GENERATE_COVERAGE_BIGWIGS GENERATE_COVERAGE_BEDGRAPHS \
    KEEP_NORMALIZATION_POLICY_BAMS KEEP_SPIKEIN_BAMS; do
    setting_value="${!setting_name:-}"
    [[ -n "$setting_value" ]] || {
        case "$setting_name" in
            KEEP_NORMALIZATION_POLICY_BAMS|KEEP_SPIKEIN_BAMS|GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS) setting_value=false ;;
            *) setting_value=true ;;
        esac
    }
    case "$setting_value" in
        true|false) ok "$setting_name: $setting_value" ;;
        *) fail "$setting_name must be true or false" ;;
    esac
done
if [[ "${GENERATE_COVERAGE_BIGWIGS:-true}" != "true" && \
      "${GENERATE_COVERAGE_BEDGRAPHS:-true}" != "true" ]] && \
   [[ "${GENERATE_CPM_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_CONSENSUS_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" || \
      "${GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true}" == "true" || \
      "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
    fail "at least one coverage format must be enabled when a coverage family is enabled"
fi
for setting_name in PERMISSIVE_MIN_MAPQ INTERMEDIATE_MIN_MAPQ; do
    setting_value="${!setting_name:-0}"
    [[ "$setting_value" =~ ^[0-9]+$ ]] \
        && ok "$setting_name: $setting_value" \
        || fail "$setting_name must be a non-negative integer"
done
for setting_name in SPIKEIN_MIN_MAPQ SPIKEIN_MIN_FRAGMENTS_FAIL SPIKEIN_MIN_FRAGMENTS_WARN; do
    setting_value="${!setting_name:-}"
    [[ -n "$setting_value" ]] || {
        case "$setting_name" in
            SPIKEIN_MIN_MAPQ) setting_value=30 ;;
            SPIKEIN_MIN_FRAGMENTS_FAIL) setting_value=1000 ;;
            *) setting_value=10000 ;;
        esac
    }
    [[ "$setting_value" =~ ^[0-9]+$ ]] \
        && ok "$setting_name: $setting_value" \
        || fail "$setting_name must be a non-negative integer"
done
if (( ${SPIKEIN_MIN_FRAGMENTS_WARN:-10000} < ${SPIKEIN_MIN_FRAGMENTS_FAIL:-1000} )); then
    fail "SPIKEIN_MIN_FRAGMENTS_WARN must be >= SPIKEIN_MIN_FRAGMENTS_FAIL"
else
    ok "Drosophila spike-in count thresholds"
fi
if awk -v value="${SPIKEIN_SCALE_TARGET:-1000000}" \
    'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value > 0)}'; then
    ok "SPIKEIN_SCALE_TARGET: ${SPIKEIN_SCALE_TARGET:-1000000}"
else
    fail "SPIKEIN_SCALE_TARGET must be a positive number"
fi
for setting_name in SPIKEIN_WARN_LOW_FRACTION SPIKEIN_WARN_HIGH_FRACTION; do
    setting_value="${!setting_name:-}"
    [[ -n "$setting_value" ]] || {
        [[ "$setting_name" == SPIKEIN_WARN_LOW_FRACTION ]] && setting_value=0.001 || setting_value=0.20
    }
    if awk -v value="$setting_value" \
        'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && value >= 0 && value <= 1)}'; then
        ok "$setting_name: $setting_value"
    else
        fail "$setting_name must be between zero and one"
    fi
done
if awk -v low="${SPIKEIN_WARN_LOW_FRACTION:-0.001}" -v high="${SPIKEIN_WARN_HIGH_FRACTION:-0.20}" \
    'BEGIN{exit !(low < high)}'; then
    ok "Drosophila spike-in fraction thresholds"
else
    fail "SPIKEIN_WARN_LOW_FRACTION must be lower than SPIKEIN_WARN_HIGH_FRACTION"
fi
if [[ "${SPIKEIN_CANONICAL_CONTIGS:-2L,2R,3L,3R,4,X}" =~ ^[A-Za-z0-9_,]+$ ]]; then
    ok "SPIKEIN_CANONICAL_CONTIGS: ${SPIKEIN_CANONICAL_CONTIGS:-2L,2R,3L,3R,4,X}"
else
    fail "SPIKEIN_CANONICAL_CONTIGS contains unsupported characters"
fi
SAMPLESHEET="$SAMPLESHEET_ARG"
[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }
python3 "$INPUT_SANITIZER" "$SAMPLESHEET"
VERSION_FILE="${SCRIPT_DIR}/../VERSION"
VERSION="$(tr -d '[:space:]' < "$VERSION_FILE" 2>/dev/null || echo unknown)"

echo "ATACseq2tracks ${VERSION} pre-flight"

for script in validate_samplesheet.py sanitize_text_inputs.py extract_spikein_samples.py fastqc_batch.sh trimgalore_batch.sh \
    bowtie2_batch.sh picard_dedup_batch.sh blacklist_filter.sh \
    blacklist_filter_batch.sh filter_bam_for_coverage_policy.sh \
    genomecoverage_single.sh genomecoverage_batch.sh merge_replicates.sh \
    macs2_peaks.sh macs2_batch.sh post_alignment_qc_batch.sh \
    consensus_peak_size_factors.R qc_table_helpers.sh track_normalization_helpers.sh parallel_job_helpers.sh \
    generate_filtering_sensitivity_tracks.sh \
    filter_composite_spikein_bam.sh drosophila_spikein_tracks.sh \
    create_ucsc_tracks.sh ataqv_qc_batch.sh \
    prepare_tss_bed.py extract_ataqv_metrics.py plot_fragment_periodicity.py \
    plot_chrom_coverage.py peak_interpretation.sh prepare_diffbind.R \
    diffbind_analysis.sh diffbind_analysis.R deseq2atac_analysis.sh \
    deseq2atac_analysis.R peak_annotation_helpers.R \
    summarize_differential_accessibility.py generate_pipeline_report.sh; do
    [[ -f "${SCRIPT_DIR}/${script}" ]] && ok "script: ${script}" || fail "missing script: ${script}"
done

if [[ "${RUN_ATAQV_QC:-true}" == "true" ]]; then
    command -v ataqv >/dev/null 2>&1 && ok "tool: ataqv" || fail "tool not in PATH: ataqv"
    if [[ "${GENERATE_ATAQV_VIEWER:-true}" == "true" ]]; then
        command -v mkarv >/dev/null 2>&1 && ok "tool: mkarv" || fail "tool not in PATH: mkarv"
    fi
    for package in matplotlib numpy; do
        python3 -c "import ${package}" >/dev/null 2>&1 \
            && ok "Python package: ${package}" || fail "Python package not installed: ${package}"
    done
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
    if [[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
        if python3 "${SCRIPT_DIR}/extract_spikein_samples.py" "$SAMPLESHEET" >/dev/null; then
            ok "dm6 spike-in declarations"
        else
            fail "invalid or absent dm6 spike-in declarations"
        fi
    fi

    mapfile -t USED_GENOMES < <(python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    print("\n".join(sorted({row["genome"].strip() for row in csv.DictReader(handle) if row["genome"].strip()})))
PY
    )
    if (( ${#USED_GENOMES[@]} > 1 )); then
        fail "multiple genome builds in one run (${USED_GENOMES[*]}); v${VERSION} requires one build per run"
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
    if [[ "${GENERATE_DESEQ2_CONSENSUS_TRACKS:-true}" == "true" || \
          "${GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true}" == "true" || \
          "${GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true}" == "true" || \
          "${GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true}" == "true" ]] && \
       (( BIOLOGICAL_SAMPLE_COUNT < ${CONSENSUS_MIN_SAMPLES:-2} )); then
        fail "DESeq2 consensus tracks require at least ${CONSENSUS_MIN_SAMPLES:-2} biological samples; found $BIOLOGICAL_SAMPLE_COUNT"
    else
        ok "biological samples: $BIOLOGICAL_SAMPLE_COUNT"
    fi

    DIFFBIND_SUMMITS_VALUE="${DIFFBIND_SUMMITS:-100}"
    if [[ "$DIFFBIND_SUMMITS_VALUE" =~ ^[0-9]+$ ]]; then
        ok "DiffBind summit half-width: ${DIFFBIND_SUMMITS_VALUE} bp"
    else
        fail "DIFFBIND_SUMMITS must be a non-negative integer"
    fi
    if awk -v value="${DIFFBIND_ALPHA:-0.05}" \
        'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && value > 0 && value < 1)}'; then
        ok "DiffBind FDR alpha: ${DIFFBIND_ALPHA:-0.05}"
    else
        fail "DIFFBIND_ALPHA must be numeric and between zero and one"
    fi
    if awk -v value="${DIFFERENTIAL_MIN_ABS_LOG2FC:-0}" \
        'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && value >= 0)}'; then
        ok "differential minimum absolute log2 fold change: ${DIFFERENTIAL_MIN_ABS_LOG2FC:-0}"
    else
        fail "DIFFERENTIAL_MIN_ABS_LOG2FC must be a non-negative number"
    fi
    if python3 - "$SAMPLESHEET" "${DIFFERENTIAL_CONDITION_ORDER:-}" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    observed = []
    seen_samples = set()
    for row in csv.DictReader(handle):
        if row["is_control"].strip().lower() in {"true", "1", "yes"}:
            continue
        key = (row["sample_id"].strip(), row["replicate"].strip())
        if key in seen_samples:
            continue
        seen_samples.add(key)
        condition = row["condition"].strip()
        if condition not in observed:
            observed.append(condition)
requested = [value.strip() for value in sys.argv[2].split(",") if value.strip()]
if len(requested) != len(set(requested)):
    raise SystemExit("DIFFERENTIAL_CONDITION_ORDER contains duplicate names")
unknown = [value for value in requested if value not in observed]
if unknown:
    raise SystemExit("Unknown condition(s) in DIFFERENTIAL_CONDITION_ORDER: " + ", ".join(unknown))
PY
    then
        ok "universal differential condition order"
    else
        fail "invalid DIFFERENTIAL_CONDITION_ORDER"
    fi

    if [[ "${RUN_DESEQ2ATAC:-true}" == "true" ]]; then
        DESEQ2ATAC_MIN_VALUE="${DESEQ2ATAC_MIN_SAMPLES:-2}"
        if [[ "$DESEQ2ATAC_MIN_VALUE" =~ ^[1-9][0-9]*$ ]] && \
           (( DESEQ2ATAC_MIN_VALUE <= BIOLOGICAL_SAMPLE_COUNT )); then
            ok "DESeq2ATAC minimum broad/narrow peak support: ${DESEQ2ATAC_MIN_VALUE}"
        else
            fail "DESEQ2ATAC_MIN_SAMPLES must be positive and no greater than biological sample count"
        fi
        if awk -v value="${DESEQ2ATAC_ALPHA:-0.05}" \
            'BEGIN{exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && value > 0 && value < 1)}'; then
            ok "DESeq2ATAC FDR alpha: ${DESEQ2ATAC_ALPHA:-0.05}"
        else
            fail "DESEQ2ATAC_ALPHA must be numeric and between zero and one"
        fi
        if python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
bad = []
seen = set()
with open(sys.argv[1], newline="") as handle:
    for row in csv.DictReader(handle):
        if row["is_control"].strip().lower() in {"true", "1", "yes"}:
            continue
        key = (row["sample_id"].strip(), row["replicate"].strip())
        if key in seen:
            continue
        seen.add(key)
        if row["macs2_mode"].strip().lower() != "both":
            bad.append("_bioR".join(key))
if bad:
    print("DESeq2ATAC broad+narrow analyses require macs2_mode=both for: " + ", ".join(bad), file=sys.stderr)
    raise SystemExit(1)
PY
        then
            ok "DESeq2ATAC broad and narrow peaks requested for every biological sample"
        else
            fail "DESeq2ATAC requires macs2_mode=both"
        fi
        if [[ -n "${DESEQ2ATAC_BLOCK_COLUMN:-}" ]]; then
            if python3 - "$SAMPLESHEET" "${DESEQ2ATAC_BLOCK_COLUMN}" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    fields = csv.DictReader(handle).fieldnames or []
raise SystemExit(0 if sys.argv[2] in fields else 1)
PY
            then
                ok "DESeq2ATAC block column: ${DESEQ2ATAC_BLOCK_COLUMN}"
            else
                fail "DESeq2ATAC block column is absent: ${DESEQ2ATAC_BLOCK_COLUMN}"
            fi
        else
            ok "DESeq2ATAC design: ~ condition"
        fi
    fi

    mapfile -t RUN_LAYOUTS < <(python3 - "$SAMPLESHEET" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    rows = csv.DictReader(handle)
    print("\n".join(sorted({r["layout"].strip().upper() for r in rows})))
PY
    )
    if (( ${#RUN_LAYOUTS[@]} > 1 )); then
        fail "mixed PE/SE run (${RUN_LAYOUTS[*]}); use separate samplesheets and output directories"
    elif (( ${#RUN_LAYOUTS[@]} == 1 )); then
        ok "run layout: ${RUN_LAYOUTS[0]}"
    fi
fi

SE_MODE="${SE_SIGNAL_MODE:-read}"
if [[ "${SE_MODE,,}" == "read" ]]; then
    ok "single-end signal mode: read"
else
    fail "unsupported SE_SIGNAL_MODE=${SE_MODE}; v${VERSION} supports read only"
fi

if [[ "${TRACK_STANDARD_CHROMS_ONLY:-true}" == "true" ]]; then
    ok "canonical chromosome universe enabled"
else
    fail "TRACK_STANDARD_CHROMS_ONLY must be true for unified track normalization"
fi

MACS3_COMMAND="${MACS3_COMMAND:-macs3}"
for tool in bowtie2 samtools bedtools trim_galore fastqc "$MACS3_COMMAND" \
    multiqc python3 "${R_BIN:-Rscript}" bamCoverage multiBamSummary gzip java sha256sum; do
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

check_complete_bowtie2_index() {
    local label="$1" prefix="$2" extension component missing
    for extension in bt2 bt2l; do
        missing=0
        for component in 1 2 3 4; do [[ -s "${prefix}.${component}.${extension}" ]] || missing=$((missing + 1)); done
        for component in 1 2; do [[ -s "${prefix}.rev.${component}.${extension}" ]] || missing=$((missing + 1)); done
        if (( missing == 0 )); then
            ok "$label: $prefix"
            return
        fi
    done
    fail "$label is incomplete or missing: $prefix"
}

check_ccre_reference() {
    local label="$1" path="$2" source="$3"
    if [[ -z "$path" ]]; then
        fail "$label is empty while RUN_CCRE_ANNOTATION=true"
        return
    fi
    if [[ ! -s "$path" ]]; then
        fail "$label not found or empty: $path"
        return
    fi
    if [[ "$path" == *.gz ]] && ! gzip -t "$path"; then
        fail "$label is not a readable gzip file: $path"
        return
    fi
    if python3 - "$path" <<'PY'
import gzip, sys
opener = gzip.open if sys.argv[1].lower().endswith(".gz") else open
valid = False
with opener(sys.argv[1], "rt", errors="replace") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        try:
            valid = len(fields) >= 3 and int(fields[1]) >= 0 and int(fields[2]) > int(fields[1])
        except ValueError:
            valid = False
        if valid:
            break
raise SystemExit(0 if valid else 1)
PY
    then
        ok "$label: $path ($source)"
    else
        fail "$label has no valid BED interval: $path"
    fi
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
            if [[ "${RUN_SIMPLE_PEAK_ANNOTATION:-true}" == "true" ]]; then
                check_reference GTF_HUMAN "${GTF_HUMAN:-}"
                if [[ "${RUN_CCRE_ANNOTATION:-true}" == "true" ]]; then
                    check_ccre_reference CCRE_BED_HG38 "${CCRE_BED_HG38:-}" \
                        "${CCRE_SOURCE_HG38:-ENCODE4_GRCh38}"
                else
                    ok "cCRE annotation disabled; using GTF-only annotation"
                fi
            fi
            index="${INDEX_HG38:-}"
            composite_index="${INDEX_HG38_DM6:-}"
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
            if [[ "${RUN_SIMPLE_PEAK_ANNOTATION:-true}" == "true" ]]; then
                check_reference GTF_MOUSE "${GTF_MOUSE:-}"
                if [[ "${RUN_CCRE_ANNOTATION:-true}" == "true" ]]; then
                    check_ccre_reference CCRE_BED_MM39 "${CCRE_BED_MM39:-}" \
                        "${CCRE_SOURCE_MM39:-ENCODE3_mm10_liftOver_mm39}"
                else
                    ok "cCRE annotation disabled; using GTF-only annotation"
                fi
            fi
            index="${INDEX_MM39:-}"
            composite_index="${INDEX_MM39_DM6:-}"
            ;;
        "") continue ;;
        *) fail "unsupported genome: $genome"; continue ;;
    esac
    if compgen -G "${index}*.bt2" >/dev/null || compgen -G "${index}*.bt2l" >/dev/null; then
        ok "Bowtie2 index: $index"
    else
        fail "Bowtie2 index not found: $index"
    fi
    if [[ "${GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS:-false}" == "true" ]]; then
        check_reference CHROM_SIZES_DM6 "${CHROM_SIZES_DM6:-}"
        check_reference BLACKLIST_DM6 "${BLACKLIST_DM6:-}"
        check_complete_bowtie2_index "Composite ${genome}+dm6 Bowtie2 index" "$composite_index"
    fi
done

for package in DESeq2 DiffBind ggplot2 dplyr rtracklayer GenomicAlignments GenomicRanges IRanges Rsamtools BiocParallel; do
    if "${R_BIN:-Rscript}" -e "quit(status=ifelse(requireNamespace('${package}', quietly=TRUE), 0, 1))" >/dev/null 2>&1; then
        ok "R package: $package"
    else
        fail "R package not installed: $package"
    fi
done

if "${R_BIN:-Rscript}" -e '
exports <- getNamespaceExports("IRanges")
quit(status=ifelse(all(c("overlapsAny", "findOverlaps") %in% exports), 0, 1))
' >/dev/null 2>&1; then
    ok "IRanges overlap generics exported"
else
    fail "IRanges must export overlapsAny and findOverlaps"
fi

echo "RESULT: ${PASS} OK | ${WARN} warnings | ${FAIL} failures"
(( FAIL == 0 )) || exit 1
