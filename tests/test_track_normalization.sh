#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_DIR}/scripts/track_normalization_helpers.sh"
PYTHON_BIN="${PYTHON_BIN:-python3}"

TEST_TMP="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-tracks.XXXXXX")"
trap 'rm -rf "$TEST_TMP"' EXIT

assert_equal() {
    local expected="$1" actual="$2" label="$3"
    if [[ "$actual" != "$expected" ]]; then
        printf 'FAIL %s: expected <%s>, got <%s>\n' "$label" "$expected" "$actual" >&2
        exit 1
    fi
    printf 'OK   %s\n' "$label"
}

declare -a args=()
deeptools_signal_args PE args
assert_equal '--extendReads --samFlagInclude 66' "${args[*]}" 'PE uses one properly paired fragment record'
deeptools_signal_args SE args
assert_equal '0' "${#args[@]}" 'SE retains read-based signal'
assert_equal 'fragment' "$(signal_unit_for_layout PE)" 'PE signal unit'
assert_equal 'read' "$(signal_unit_for_layout SE)" 'SE signal unit'
assert_equal 'read' "$(configured_se_signal_mode)" 'SE read mode is explicit'
SE_SIGNAL_MODE='tn5'
if deeptools_signal_args SE args >/dev/null 2>&1; then
    echo 'FAIL unsupported Tn5 SE mode was accepted' >&2
    exit 1
fi
SE_SIGNAL_MODE='read'
echo 'OK   unsupported SE signal modes are rejected'

PEAKS="${TEST_TMP}/input.bed"
CANONICAL="${TEST_TMP}/canonical.bed"
printf 'chr1\t10\t20\nchr19\t20\t30\nchr20\t30\t40\nchrX\t40\t50\nchrY\t50\t60\nchrM\t60\t70\nchrUn_A\t70\t80\n' > "$PEAKS"
canonicalize_peak_file "$PEAKS" "$CANONICAL" mm39
assert_equal $'chr1\t10\t20\nchr19\t20\t30\nchrX\t40\t50\nchrY\t50\t60' "$(cat "$CANONICAL")" \
    'mm39 canonical peak universe excludes chr20, mitochondrial and noncanonical contigs'
grep -q 'canonical_contigs_from_bam' "${REPO_DIR}/scripts/genomecoverage_single.sh" \
    || { echo 'FAIL CPM track does not apply the shared canonical universe' >&2; exit 1; }
[[ "$(grep -c 'canonical_contigs_from_bam' "${REPO_DIR}/scripts/post_alignment_qc_batch.sh")" -ge 2 ]] \
    || { echo 'FAIL post-alignment bigWig/bedGraph paths do not share the canonical universe' >&2; exit 1; }
grep -q 'canonicalize_peak_file' "${REPO_DIR}/scripts/post_alignment_qc_batch.sh" \
    || { echo 'FAIL consensus peaks do not apply the shared canonical universe' >&2; exit 1; }
echo 'OK   shared canonical chromosome universe contract'

size_factor='0.8'
cohort_geomean='1000'
consensus_scale="$(awk "BEGIN{print 1/$size_factor}")"
robust_scale="$(awk "BEGIN{print 1000000/($size_factor*$cohort_geomean)}")"
cohort_constant="$(awk "BEGIN{print 1000000/$cohort_geomean}")"
assert_equal "$robust_scale" "$(awk "BEGIN{print $consensus_scale*$cohort_constant}")" \
    'robust CPM is a cohort-wide rescaling of DESeq2 consensus signal'

GENOME_SCRIPT="${REPO_DIR}/scripts/genomecoverage_single.sh"
POST_SCRIPT="${REPO_DIR}/scripts/post_alignment_qc_batch.sh"
grep -q '_CPM.bw' "$GENOME_SCRIPT" || { echo 'FAIL CPM bigWig name is missing' >&2; exit 1; }
if grep -q '_CPM.bedGraph' "$GENOME_SCRIPT"; then
    echo 'FAIL CPM bedGraph must not be generated' >&2
    exit 1
fi
grep -q '_DESeq2Consensus.bw' "$POST_SCRIPT" || { echo 'FAIL DESeq2 consensus bigWig name is missing' >&2; exit 1; }
grep -q '_DESeq2Consensus.bedGraph' "$POST_SCRIPT" || { echo 'FAIL DESeq2 consensus bedGraph name is missing' >&2; exit 1; }
grep -q '_DESeq2RobustCPM.bw' "$POST_SCRIPT" || { echo 'FAIL robust CPM bigWig name is missing' >&2; exit 1; }
grep -q '_DESeq2RobustCPM.bedGraph' "$POST_SCRIPT" || { echo 'FAIL robust CPM bedGraph name is missing' >&2; exit 1; }
grep -q 'fpm(dds, robust = TRUE)' "${REPO_DIR}/scripts/consensus_peak_size_factors.R" \
    || { echo 'FAIL DESeq2 robust FPM equivalence check is missing' >&2; exit 1; }
echo 'OK   requested track output contract'

# Execute the CPM generator with small command mocks. This verifies that PE
# coverage and the CPM denominator use the same single-record fragment
# definition, the expected output is non-empty, and no CPM bedGraph appears.
MOCK_BIN="${TEST_TMP}/mock-bin"
MOCK_OUT="${TEST_TMP}/mock-output"
MOCK_LOG="${TEST_TMP}/mock-commands.log"
mkdir -p "$MOCK_BIN" "$MOCK_OUT"
export MOCK_LOG

cat > "${MOCK_BIN}/samtools" <<'MOCK_SAMTOOLS'
#!/usr/bin/env bash
set -euo pipefail
printf 'samtools %s\n' "$*" >> "$MOCK_LOG"
case "${1:-}" in
    quickcheck) exit 0 ;;
    idxstats)
        printf 'chr1\t1000\t4\t0\nchrM\t100\t2\t0\nchrUn_A\t50\t1\t0\n*\t0\t0\t0\n'
        ;;
    view)
        for arg in "$@"; do
            [[ "$arg" == '-c' ]] && { printf '2\n'; exit 0; }
        done
        while [[ $# -gt 0 ]]; do
            if [[ "$1" == '-o' ]]; then
                printf 'mock-bam\n' > "$2"
                exit 0
            fi
            shift
        done
        ;;
    index)
        bam="${@: -1}"
        printf 'mock-index\n' > "${bam}.bai"
        ;;
esac
MOCK_SAMTOOLS

cat > "${MOCK_BIN}/bamCoverage" <<'MOCK_BAMCOVERAGE'
#!/usr/bin/env bash
set -euo pipefail
printf 'bamCoverage %s\n' "$*" >> "$MOCK_LOG"
output=''
while [[ $# -gt 0 ]]; do
    if [[ "$1" == '--outFileName' ]]; then output="$2"; shift 2; else shift; fi
done
[[ -n "$output" ]] || exit 1
printf 'mock-bigwig\n' > "$output"
MOCK_BAMCOVERAGE
chmod +x "${MOCK_BIN}/samtools" "${MOCK_BIN}/bamCoverage"

MOCK_CONFIG="${TEST_TMP}/config.conf"
cat > "$MOCK_CONFIG" <<'MOCK_CONFIG_EOF'
TRACK_STANDARD_CHROMS_ONLY=true
TRACK_BIN_SIZE=10
THREADS_BIGWIG=1
BAMCOVERAGE_COMMON_ARGS=""
SE_SIGNAL_MODE="read"
MOCK_CONFIG_EOF

MOCK_BAM="${TEST_TMP}/sample_dedup_blFilt.bam"
printf 'mock-input\n' > "$MOCK_BAM"
PATH="${MOCK_BIN}:$PATH" F2T_CONFIG="$MOCK_CONFIG" \
    bash "$GENOME_SCRIPT" "$MOCK_BAM" mm39 "$MOCK_OUT" PE >/dev/null

[[ -s "${MOCK_OUT}/sample_dedup_blFilt_CPM.bw" ]] \
    || { echo 'FAIL mocked CPM bigWig is empty or missing' >&2; exit 1; }
if find "$MOCK_OUT" -maxdepth 1 -name '*_CPM.bedGraph' -print -quit | grep -q .; then
    echo 'FAIL mocked CPM generation produced a bedGraph' >&2
    exit 1
fi
grep -q 'samtools view -c -f 66 -F 4' "$MOCK_LOG" \
    || { echo 'FAIL PE CPM denominator did not count first mates once' >&2; exit 1; }
grep -q 'bamCoverage .*--normalizeUsing CPM .*--extendReads --samFlagInclude 66' "$MOCK_LOG" \
    || { echo 'FAIL PE CPM coverage did not use fragment arguments' >&2; exit 1; }
grep -q $'fragment\tCPM\t2\t' "${MOCK_OUT}/sample_dedup_blFilt_CPM.metadata.tsv" \
    || { echo 'FAIL CPM metadata lacks the fragment normalization count' >&2; exit 1; }
echo 'OK   mocked PE CPM generation and non-empty output'

MOCK_SE_BAM="${TEST_TMP}/sample_se_dedup_blFilt.bam"
printf 'mock-se-input\n' > "$MOCK_SE_BAM"
PATH="${MOCK_BIN}:$PATH" F2T_CONFIG="$MOCK_CONFIG" \
    bash "$GENOME_SCRIPT" "$MOCK_SE_BAM" mm39 "$MOCK_OUT" SE >/dev/null
[[ -s "${MOCK_OUT}/sample_se_dedup_blFilt_CPM.bw" ]] \
    || { echo 'FAIL mocked SE CPM bigWig is empty or missing' >&2; exit 1; }
grep -q 'samtools view -c -F 4 .*sample_se.*canonical.bam' "$MOCK_LOG" \
    || { echo 'FAIL SE CPM denominator did not count retained reads once' >&2; exit 1; }
SE_BAMCOVERAGE_LOG="$(grep 'bamCoverage .*sample_se_dedup_blFilt_CPM.bw' "$MOCK_LOG")"
[[ "$SE_BAMCOVERAGE_LOG" != *'--extendReads'* && "$SE_BAMCOVERAGE_LOG" != *'--samFlagInclude'* ]] \
    || { echo 'FAIL SE read coverage was artificially extended or pair-filtered' >&2; exit 1; }
grep -q $'read\tCPM\t2\t' "${MOCK_OUT}/sample_se_dedup_blFilt_CPM.metadata.tsv" \
    || { echo 'FAIL SE CPM metadata lacks the read normalization count' >&2; exit 1; }
echo 'OK   mocked SE read-based CPM generation and non-empty output'

BAMCOVERAGE_COMMON_ARGS=''
for spec in \
    'sample_DESeq2Consensus.bw bigwig 1.25' \
    'sample_DESeq2Consensus.bedGraph bedgraph 1.25' \
    'sample_DESeq2RobustCPM.bw bigwig 1250' \
    'sample_DESeq2RobustCPM.bedGraph bedgraph 1250'; do
    read -r filename format scale <<< "$spec"
    PATH="${MOCK_BIN}:$PATH" write_scaled_coverage_track \
        "$MOCK_BAM" "${MOCK_OUT}/${filename}" "$format" "$scale" PE 10 1
    [[ -s "${MOCK_OUT}/${filename}" ]] \
        || { echo "FAIL requested scaled output is empty: ${filename}" >&2; exit 1; }
done
echo 'OK   non-empty DESeq2 bigWig and bedGraph output contract'

ENTRYPOINT="${REPO_DIR}/atacseq2tracks.sh"
for stage in 7 8 10 13 14; do
    grep -q "if is_done ${stage}" "$ENTRYPOINT" \
        || { echo "FAIL checkpoint guard missing for Step ${stage}" >&2; exit 1; }
    grep -q "mark_done ${stage}" "$ENTRYPOINT" \
        || { echo "FAIL checkpoint completion missing for Step ${stage}" >&2; exit 1; }
done
grep -q '_dedup_blFilt_CPM.bw.*SKIP' "${REPO_DIR}/scripts/genomecoverage_batch.sh" \
    || { echo 'FAIL existing non-empty CPM output skip behavior is missing' >&2; exit 1; }
echo 'OK   checkpoint/resume and existing-output skip contract'

"$PYTHON_BIN" "${REPO_DIR}/scripts/validate_samplesheet.py" \
    "${REPO_DIR}/config/samplesheet_example_atac_se.csv" >/dev/null \
    || { echo 'FAIL SE ATAC example samplesheet is invalid' >&2; exit 1; }

MIXED_SAMPLESHEET="${TEST_TMP}/mixed_layouts.csv"
cat > "$MIXED_SAMPLESHEET" <<'MIXED_EOF'
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
sample_pe,/data/pe_R1.fastq.gz,/data/pe_R2.fastq.gz,PE,mm39,ATAC-seq,accessibility,WT,none,cell,1,1,FALSE,,narrow,/ref/mm39.blacklist.bed,sample_pe
sample_se,/data/se.fastq.gz,,SE,mm39,ATAC-seq,accessibility,WT,none,cell,2,1,FALSE,,narrow,/ref/mm39.blacklist.bed,sample_se
MIXED_EOF
if "$PYTHON_BIN" "${REPO_DIR}/scripts/validate_samplesheet.py" "$MIXED_SAMPLESHEET" >/dev/null 2>&1; then
    echo 'FAIL mixed PE/SE samplesheet was accepted' >&2
    exit 1
fi
echo 'OK   SE example is valid and mixed-layout samplesheets are rejected'
