#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_DIR}/scripts/qc_table_helpers.sh"

TEST_TMP="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-hotfix.XXXXXX")"
trap 'rm -rf "$TEST_TMP"' EXIT

assert_equal() {
    local expected="$1" actual="$2" label="$3"
    if [[ "$actual" != "$expected" ]]; then
        printf 'FAIL %s: expected <%s>, got <%s>\n' "$label" "$expected" "$actual" >&2
        exit 1
    fi
    printf 'OK   %s\n' "$label"
}

PICARD_METRICS="${TEST_TMP}/dup_metrics.txt"
printf '%s\n' '## METRICS CLASS picard.sam.DuplicationMetrics' > "$PICARD_METRICS"
printf 'LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n' >> "$PICARD_METRICS"
printf '\t0\t100\t0\t0\t0\t20\t1\t0.2375\t2300468\n' >> "$PICARD_METRICS"

duplication_pct="$(picard_duplication_pct "$PICARD_METRICS")"
assert_equal '23.75' "$duplication_pct" 'Picard percentage is selected by header name'

SIZE_FACTORS="${TEST_TMP}/consensus_sizeFactors.tsv"
printf 'sample_id\tkey\tgenome\tsize_factor\ttrack_scale_factor\tconsensus_count_sum\tcohort_geometric_mean_column_sum\tdeseq2_robust_cpm_scale\tbasis\n' > "$SIZE_FACTORS"
printf 'sample-a\tsample-a_bioR1\tmm39\t0.8125\t1.2307692308\t800\t1000\t1230.7692308\tconsensus_peak_counts\n' >> "$SIZE_FACTORS"
printf 'sample-b\tsample-b_bioR1\tmm39\t1.25\t0.8\t1250\t1000\t800\tconsensus_peak_counts\n' >> "$SIZE_FACTORS"

size_factor="$(consensus_size_factor_for_key "$SIZE_FACTORS" 'sample-a_bioR1')"
assert_equal '0.8125' "$size_factor" 'DESeq2 size factor is selected by header name'

robust_scale="$(consensus_table_value_for_key "$SIZE_FACTORS" 'sample-a_bioR1' deseq2_robust_cpm_scale)"
assert_equal '1230.7692308' "$robust_scale" 'robust CPM scale is selected by header name'

if consensus_size_factor_for_key "$SIZE_FACTORS" 'missing_bioR1' >/dev/null 2>&1; then
    echo 'FAIL missing DESeq2 key was accepted' >&2
    exit 1
fi
echo 'OK   missing DESeq2 key is rejected'

printf 'sample_id\tkey\tgenome\tsize_factor\ninvalid\tinvalid_bioR1\tmm39\t0\n' > "$SIZE_FACTORS"
if consensus_size_factor_for_key "$SIZE_FACTORS" 'invalid_bioR1' >/dev/null 2>&1; then
    echo 'FAIL non-positive DESeq2 size factor was accepted' >&2
    exit 1
fi
echo 'OK   non-positive DESeq2 size factor is rejected'
