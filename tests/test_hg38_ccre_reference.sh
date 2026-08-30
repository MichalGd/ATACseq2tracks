#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UTILITY="${REPO_DIR}/utilities/prepare_encode4_hg38_ccre.sh"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

grep -q 'Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz' "$UTILITY"
grep -q 'EXPECTED_RECORDS=2348854' "$UTILITY"
grep -q 'GRCh38/hg38' "$UTILITY"

printf 'chr1\t100\t200\tEH38D1\tEH38E1\tpELS,CTCF-bound\nchrX\t300\t450\tEH38D2\tEH38E2\tPLS\n' \
    | gzip -c > "$TMP_DIR/valid.bed.gz"
bash "$UTILITY" --validate-only "$TMP_DIR/valid.bed.gz" >/dev/null

printf 'chr1\t100\t200\tEH38D1\tEH38E1\tunknown-class\n' \
    | gzip -c > "$TMP_DIR/invalid-class.bed.gz"
if bash "$UTILITY" --validate-only "$TMP_DIR/invalid-class.bed.gz" >/dev/null 2>&1; then
    echo 'FAIL unsupported cCRE class was accepted' >&2
    exit 1
fi

printf 'chrM\t100\t200\tEH38D1\tEH38E1\tpELS\n' \
    | gzip -c > "$TMP_DIR/noncanonical.bed.gz"
if bash "$UTILITY" --validate-only "$TMP_DIR/noncanonical.bed.gz" >/dev/null 2>&1; then
    echo 'FAIL noncanonical cCRE interval was accepted' >&2
    exit 1
fi

echo 'OK   ENCODE4 hg38 cCRE reference utility validation'
