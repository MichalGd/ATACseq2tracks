#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SANITIZER="${REPO_DIR}/scripts/sanitize_text_inputs.py"
TEST_TMP="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-inputs.XXXXXX")"
trap 'rm -rf "$TEST_TMP"' EXIT

CONFIG="${TEST_TMP}/config.conf"
SAMPLESHEET="${TEST_TMP}/samplesheet.csv"
CLEAN="${TEST_TMP}/clean.conf"
INVALID="${TEST_TMP}/invalid.conf"
PYTHON_BIN="${PYTHON_BIN:-python3}"

"$PYTHON_BIN" - "$CONFIG" "$SAMPLESHEET" "$CLEAN" "$INVALID" <<'PY'
import codecs
import sys
from pathlib import Path

config, samplesheet, clean, invalid = map(Path, sys.argv[1:])
config.write_bytes(
    codecs.BOM_UTF8
    + b'SAMPLESHEET="/tmp/samplesheet.csv"\r\n'
    + b'OUTPUT_DIR="/tmp/output"\r\n'
    + b'THREADS_PARALLEL_JOBS=2\r\n'
)
samplesheet.write_bytes(
    (
        "sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,"
        "treatment,cell_type,replicate,tech_replicate,is_control,control_id,"
        "macs2_mode,blacklist,output_prefix\r\n"
        "sample-a,/tmp/a_R1.fastq.gz,/tmp/a_R2.fastq.gz,PE,mm39,ATAC-seq,"
        "accessibility,case,none,cell,1,1,FALSE,,both,/tmp/blacklist.bed,"
        "sample-a\r\n"
    ).encode("utf-16")
)
clean.write_bytes(b'OUTPUT_DIR="/tmp/clean"\n')
invalid.write_bytes(b'OUTPUT_DIR="/tmp/\xff"\r\n')
PY

cp "$CONFIG" "${CONFIG}.expected-original"
cp "$SAMPLESHEET" "${SAMPLESHEET}.expected-original"
cp "$INVALID" "${INVALID}.expected-original"

"$PYTHON_BIN" "$SANITIZER" "$CONFIG" "$SAMPLESHEET" "$CLEAN"

mapfile -t config_backups < <(compgen -G "${CONFIG}.windows-artifact-backup.*" || true)
mapfile -t samplesheet_backups < <(compgen -G "${SAMPLESHEET}.windows-artifact-backup.*" || true)
mapfile -t clean_backups < <(compgen -G "${CLEAN}.windows-artifact-backup.*" || true)

[[ ${#config_backups[@]} -eq 1 ]] || { echo "FAIL expected one config backup" >&2; exit 1; }
[[ ${#samplesheet_backups[@]} -eq 1 ]] || { echo "FAIL expected one samplesheet backup" >&2; exit 1; }
[[ ${#clean_backups[@]} -eq 0 ]] || { echo "FAIL clean file received a backup" >&2; exit 1; }
cmp -s "${CONFIG}.expected-original" "${config_backups[0]}" \
    || { echo "FAIL config backup is not byte-for-byte identical" >&2; exit 1; }
cmp -s "${SAMPLESHEET}.expected-original" "${samplesheet_backups[0]}" \
    || { echo "FAIL samplesheet backup is not byte-for-byte identical" >&2; exit 1; }

"$PYTHON_BIN" - "$CONFIG" "$SAMPLESHEET" <<'PY'
import sys
from pathlib import Path

for filename in sys.argv[1:]:
    payload = Path(filename).read_bytes()
    assert not payload.startswith((b"\xef\xbb\xbf", b"\xff\xfe", b"\xfe\xff"))
    assert b"\r" not in payload
    payload.decode("utf-8")
PY

# A second pass must be idempotent and must not create more backups.
"$PYTHON_BIN" "$SANITIZER" "$CONFIG" "$SAMPLESHEET" >/dev/null
[[ $(compgen -G "${CONFIG}.windows-artifact-backup.*" | wc -l) -eq 1 ]]
[[ $(compgen -G "${SAMPLESHEET}.windows-artifact-backup.*" | wc -l) -eq 1 ]]

# Unknown encodings are rejected without changing or backing up the input.
if "$PYTHON_BIN" "$SANITIZER" "$INVALID" >/dev/null 2>&1; then
    echo "FAIL invalid encoding was accepted" >&2
    exit 1
fi
cmp -s "$INVALID" "${INVALID}.expected-original" \
    || { echo "FAIL invalid input was modified" >&2; exit 1; }
if compgen -G "${INVALID}.windows-artifact-backup.*" >/dev/null; then
    echo "FAIL invalid input received an unnecessary backup" >&2
    exit 1
fi

echo "OK   Windows artifacts are backed up and normalized"
echo "OK   clean inputs are unchanged and normalization is idempotent"
echo "OK   unknown encodings are rejected without modification"
