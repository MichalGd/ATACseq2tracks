#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${REPO_DIR}/config/config.conf"
ENTRYPOINT="${REPO_DIR}/atacseq2tracks.sh"
FILTER="${REPO_DIR}/scripts/filter_bam_for_coverage_policy.sh"
POLICY_MODULE="${REPO_DIR}/scripts/generate_filtering_sensitivity_tracks.sh"
POST="${REPO_DIR}/scripts/post_alignment_qc_batch.sh"
HELPERS="${REPO_DIR}/scripts/track_normalization_helpers.sh"

assert_grep() {
    local pattern="$1" file="$2" label="$3"
    grep -Eq -- "$pattern" "$file" || { echo "FAIL $label" >&2; exit 1; }
    echo "OK   $label"
}

[[ "$(tr -d '[:space:]' < "${REPO_DIR}/VERSION")" == "4.2.0" ]] \
    || { echo 'FAIL VERSION is not 4.2.0' >&2; exit 1; }
echo 'OK   v4.2.0 preserves v4.1.0 coverage-policy contract'

for setting in \
    GENERATE_CPM_TRACKS \
    GENERATE_DESEQ2_CONSENSUS_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS \
    GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS \
    GENERATE_COVERAGE_BIGWIGS \
    GENERATE_COVERAGE_BEDGRAPHS; do
    assert_grep "^${setting}=true$" "$CONFIG" "$setting defaults to enabled"
done
assert_grep '^PERMISSIVE_MIN_MAPQ=0$' "$CONFIG" 'permissive MAPQ is explicit'
assert_grep '^INTERMEDIATE_MIN_MAPQ=0$' "$CONFIG" 'intermediate MAPQ is explicit'
assert_grep '^KEEP_NORMALIZATION_POLICY_BAMS=false$' "$CONFIG" 'policy BAM cleanup defaults on'

assert_grep 'GENERATE_CPM_TRACKS:-true' "$ENTRYPOINT" 'CPM family has a runtime switch'
assert_grep 'GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS:-true' "$ENTRYPOINT" 'permissive family has a runtime switch'
assert_grep 'GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS:-true' "$ENTRYPOINT" 'intermediate family has a runtime switch'
assert_grep 'GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS:-true' "$POST" 'stringent family has a runtime switch'
assert_grep 'GENERATE_COVERAGE_BIGWIGS:-true' "$POLICY_MODULE" 'bigWig format has a runtime switch'
assert_grep 'GENERATE_COVERAGE_BEDGRAPHS:-true' "$POLICY_MODULE" 'bedGraph format has a runtime switch'

assert_grep '-f 2 -F 2828' "$FILTER" 'PE policy filter keeps primary proper pairs and excludes secondary/supplementary records'
assert_grep '-F 2820' "$FILTER" 'SE policy filter excludes secondary/supplementary records'
if grep -Eq -- '-F (3852|3844)' "$FILTER"; then
    echo 'FAIL policy filter unexpectedly excludes duplicate flags directly' >&2
    exit 1
fi
echo 'OK   duplicate retention/removal is defined by the source BAM branch'

assert_grep 'source_dir="\$\{OUTPUT_DIR\}/bams"' "$POLICY_MODULE" 'permissive branch uses pre-dedup BAMs'
assert_grep 'source_dir="\$\{OUTPUT_DIR\}/dedupBams"' "$POLICY_MODULE" 'intermediate branch uses Picard-deduplicated BAMs'
assert_grep 'CONSENSUS_PEAK' "$POLICY_MODULE" 'all policies consume the fixed consensus peak BED'
assert_grep 'consensus_peak_sha256' "$POLICY_MODULE" 'metadata records the consensus peak checksum'
assert_grep 'consensus_peak_size_factors.R' "$POLICY_MODULE" 'each additional policy estimates its own DESeq2 factors'

grep -Fq '_DESeq2RobustCPM_${title}' "$POLICY_MODULE" \
    || { echo 'FAIL permissive/intermediate robust filename template is missing' >&2; exit 1; }
echo 'OK   permissive/intermediate robust filename template is present'
assert_grep '_DESeq2RobustCPM_Stringent.bw' "$POST" 'explicit stringent bigWig name is present'
assert_grep '_DESeq2RobustCPM_Stringent.bedGraph' "$POST" 'explicit stringent bedGraph name is present'
assert_grep '_DESeq2RobustCPM.bw' "$POST" 'legacy robust bigWig alias is retained'

TEST_TMP="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-policy-test.XXXXXX")"
trap 'rm -rf "$TEST_TMP"' EXIT
mkdir -p "${TEST_TMP}/bin"
cat > "${TEST_TMP}/bin/samtools" <<'MOCK_SAMTOOLS'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" == view ]] || exit 1
shift
if [[ " $* " == *' -c '* ]]; then
    if [[ " $* " == *' -q 30 '* ]]; then printf '4\n'
    elif [[ " $* " == *' -q 10 '* ]]; then printf '6\n'
    elif [[ " $* " == *' -q 1 '* ]]; then printf '8\n'
    else printf '10\n'
    fi
else
    printf 'r1\t66\tchr1\t1\t0\t50M\t=\t101\t150\tA\tI\tXS:i:0\n'
    printf 'r2\t66\tchr1\t2\t30\t50M\t=\t102\t150\tA\tI\n'
    printf 'r3\t66\tchr1\t3\t10\t50M\t=\t103\t150\tA\tI\tXS:i:-1\n'
fi
MOCK_SAMTOOLS
chmod +x "${TEST_TMP}/bin/samtools"
# shellcheck disable=SC1090
source "$HELPERS"
observed="$(PATH="${TEST_TMP}/bin:$PATH" signal_mapq_bin_counts mock.bam PE)"
[[ "$observed" == $'2\t2\t2\t4\t2' ]] \
    || { echo "FAIL MAPQ-bin fragment diagnostics: $observed" >&2; exit 1; }
echo 'OK   MAPQ bins and XS diagnostics count PE fragments once'

assert_grep 'KEEP_NORMALIZATION_POLICY_BAMS:-false' "$ENTRYPOINT" 'temporary policy BAM cleanup is configurable'
assert_grep 'coverage_filtering_policy_bams.*-type f' "$ENTRYPOINT" 'policy BAM cleanup is narrowly targeted'
assert_grep 'if is_done 10b' "$ENTRYPOINT" 'filtering sensitivity stage is independently checkpointed'

# Exercise the policy orchestrator with small command mocks. Existing mock
# policy BAMs deliberately take the resume path; filtering flags are tested
# statically above, while this verifies matrices, factors, tracks and metadata.
MOCK_ROOT="${TEST_TMP}/policy-module"
MOCK_BIN="${MOCK_ROOT}/bin"
MOCK_OUTPUT="${MOCK_ROOT}/output"
mkdir -p "$MOCK_BIN" "$MOCK_OUTPUT/bams" "$MOCK_OUTPUT/dedupBams" \
    "$MOCK_OUTPUT/coverage_filtering_policy_bams/permissive" \
    "$MOCK_OUTPUT/coverage_filtering_policy_bams/intermediate"

for key in sampleA_bioR1 sampleB_bioR1; do
    printf 'bam\n' > "${MOCK_OUTPUT}/bams/${key}.bam"
    printf 'bam\n' > "${MOCK_OUTPUT}/dedupBams/${key}_dedup.bam"
    printf 'bam\n' > "${MOCK_OUTPUT}/coverage_filtering_policy_bams/permissive/${key}_permissive.bam"
    printf 'bam\n' > "${MOCK_OUTPUT}/coverage_filtering_policy_bams/intermediate/${key}_intermediate.bam"
done
printf 'chr1\t10\t30\n' > "${MOCK_ROOT}/consensus.bed"
printf 'chr1\t0\t5\n' > "${MOCK_ROOT}/blacklist.bed"
cat > "${MOCK_ROOT}/samplesheet.csv" <<EOF
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
sampleA,/data/A_R1.fastq.gz,/data/A_R2.fastq.gz,PE,mm39,ATAC-seq,accessibility,A,none,cell,1,1,FALSE,,both,${MOCK_ROOT}/blacklist.bed,sampleA
sampleB,/data/B_R1.fastq.gz,/data/B_R2.fastq.gz,PE,mm39,ATAC-seq,accessibility,B,none,cell,1,1,FALSE,,both,${MOCK_ROOT}/blacklist.bed,sampleB
EOF

cat > "${MOCK_BIN}/samtools" <<'MOCK_POLICY_SAMTOOLS'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
    quickcheck) exit 0 ;;
    view)
        shift
        if [[ " $* " == *' -c '* ]]; then
            if [[ " $* " == *' -q 30 '* ]]; then printf '4\n'
            elif [[ " $* " == *' -q 10 '* ]]; then printf '6\n'
            elif [[ " $* " == *' -q 1 '* ]]; then printf '8\n'
            else printf '10\n'
            fi
        else
            printf 'r1\t66\tchr1\t1\t0\t50M\t=\t101\t150\tA\tI\tXS:i:0\n'
            printf 'r2\t66\tchr1\t2\t30\t50M\t=\t102\t150\tA\tI\n'
        fi
        ;;
    *) exit 0 ;;
esac
MOCK_POLICY_SAMTOOLS

cat > "${MOCK_BIN}/multiBamSummary" <<'MOCK_MULTIBAM'
#!/usr/bin/env bash
set -euo pipefail
npz=''; tab=''
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) npz="$2"; shift 2 ;;
        --outRawCounts) tab="$2"; shift 2 ;;
        *) shift ;;
    esac
done
printf 'npz\n' > "$npz"
printf '#chr\tstart\tend\tsampleA_bioR1\tsampleB_bioR1\nchr1\t10\t30\t10\t20\n' > "$tab"
MOCK_MULTIBAM

cat > "${MOCK_BIN}/Rscript" <<'MOCK_RSCRIPT'
#!/usr/bin/env bash
set -euo pipefail
table_dir="$4"; raw="$5"; norm="$6"
mkdir -p "$table_dir"
printf 'sample_id\tkey\tgenome\tsize_factor\ttrack_scale_factor\tdeseq2_consensus_scale\tconsensus_count_sum\tcohort_geometric_mean_column_sum\trobust_effective_library_size\tdeseq2_robust_cpm_scale\tbasis\n' > "${table_dir}/consensus_sizeFactors.tsv"
printf 'sampleA\tsampleA_bioR1\tmm39\t0.5\t2\t2\t10\t14.1421356\t7.0710678\t141421.356\tDESeq2_poscounts_on_consensus_peaks\n' >> "${table_dir}/consensus_sizeFactors.tsv"
printf 'sampleB\tsampleB_bioR1\tmm39\t1.5\t0.6666667\t0.6666667\t20\t14.1421356\t21.2132034\t47140.452\tDESeq2_poscounts_on_consensus_peaks\n' >> "${table_dir}/consensus_sizeFactors.tsv"
printf 'region\tsampleA_bioR1\tsampleB_bioR1\nchr1:10-30\t10\t20\n' > "$raw"
cp "$raw" "$norm"
MOCK_RSCRIPT

cat > "${MOCK_BIN}/bamCoverage" <<'MOCK_POLICY_BAMCOVERAGE'
#!/usr/bin/env bash
set -euo pipefail
output=''
while [[ $# -gt 0 ]]; do
    if [[ "$1" == --outFileName ]]; then output="$2"; shift 2; else shift; fi
done
printf 'track\n' > "$output"
MOCK_POLICY_BAMCOVERAGE
chmod +x "${MOCK_BIN}/"*

cat > "${MOCK_ROOT}/config.conf" <<'MOCK_POLICY_CONFIG'
TRACK_PARALLEL_JOBS=2
COVERAGE_FILTER_PARALLEL_JOBS=2
THREADS_DEEPTOOLS=1
THREADS_BIGWIG=1
TRACK_BIN_SIZE=10
TRACK_STANDARD_CHROMS_ONLY=true
SE_SIGNAL_MODE=read
R_BIN=Rscript
GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS=true
GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS=true
GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS=false
GENERATE_COVERAGE_BIGWIGS=true
GENERATE_COVERAGE_BEDGRAPHS=true
BAMCOVERAGE_COMMON_ARGS=""
MOCK_POLICY_CONFIG

PATH="${MOCK_BIN}:$PATH" F2T_CONFIG="${MOCK_ROOT}/config.conf" \
    bash "$POLICY_MODULE" "${MOCK_ROOT}/samplesheet.csv" "$MOCK_OUTPUT" \
        "${MOCK_ROOT}/consensus.bed" >/dev/null

for policy in permissive intermediate; do
    title="${policy^}"
    for key in sampleA_bioR1 sampleB_bioR1; do
        [[ -s "${MOCK_OUTPUT}/bigwig_deseq2_robust_cpm/${policy}/${key}_DESeq2RobustCPM_${title}.bw" ]] \
            || { echo "FAIL missing mocked ${policy} bigWig for ${key}" >&2; exit 1; }
        [[ -s "${MOCK_OUTPUT}/bigwig_deseq2_robust_cpm/${policy}/${key}_DESeq2RobustCPM_${title}.bedGraph" ]] \
            || { echo "FAIL missing mocked ${policy} bedGraph for ${key}" >&2; exit 1; }
    done
    [[ -s "${MOCK_OUTPUT}/coverage_filtering_sensitivity/${policy}/tables/consensus_sizeFactors.tsv" ]] \
        || { echo "FAIL missing mocked ${policy} factor table" >&2; exit 1; }
done
[[ "$(wc -l < "${MOCK_OUTPUT}/coverage_filtering_sensitivity/track_normalization_metadata.tsv")" -eq 5 ]] \
    || { echo 'FAIL combined policy metadata does not contain four samples plus header' >&2; exit 1; }
[[ "$(tail -n +2 "${MOCK_OUTPUT}/coverage_filtering_sensitivity/track_normalization_metadata.tsv" | cut -f18 | sort -u | wc -l)" -eq 1 ]] \
    || { echo 'FAIL filtering policies did not record one fixed consensus checksum' >&2; exit 1; }
echo 'OK   mocked policy module writes policy-specific factors, tracks and combined metadata'

echo 'OK   v4.1.0 configurable coverage-policy regression checks'
