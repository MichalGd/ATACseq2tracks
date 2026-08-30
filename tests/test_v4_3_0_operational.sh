#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-v430.XXXXXX")"
trap 'rm -rf -- "$TMP_ROOT"' EXIT

touch "$TMP_ROOT/a_R1.fastq.gz" "$TMP_ROOT/a_R2.fastq.gz" \
      "$TMP_ROOT/b_R1.fastq.gz" "$TMP_ROOT/b_R2.fastq.gz" "$TMP_ROOT/blacklist.bed"
printf 'x' > "$TMP_ROOT/a_R1.fastq.gz"
printf 'x' > "$TMP_ROOT/a_R2.fastq.gz"
printf 'x' > "$TMP_ROOT/b_R1.fastq.gz"
printf 'x' > "$TMP_ROOT/b_R2.fastq.gz"
cat > "$TMP_ROOT/samplesheet.csv" <<EOF
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
sampleA,$TMP_ROOT/a_R1.fastq.gz,$TMP_ROOT/a_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,conditionA,none,cells,1,1,FALSE,,both,$TMP_ROOT/blacklist.bed,sampleA
sampleA,$TMP_ROOT/b_R1.fastq.gz,$TMP_ROOT/b_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,conditionA,none,cells,1,2,FALSE,,both,$TMP_ROOT/blacklist.bed,sampleA
EOF
cat > "$TMP_ROOT/config.conf" <<EOF
SAMPLESHEET="$TMP_ROOT/samplesheet.csv"
OUTPUT_DIR="$TMP_ROOT/output"
EOF

"$PYTHON_BIN" "$REPO_DIR/scripts/resolve_config.py" \
    --config "$TMP_ROOT/config.conf" --template "$REPO_DIR/config/config.conf" \
    --shell-output "$TMP_ROOT/resolved.conf" --tsv-output "$TMP_ROOT/resolved.tsv"
grep -Fq "SAMPLESHEET=" "$TMP_ROOT/resolved.conf"
echo "OK   literal configuration is resolved without shell execution"

cat > "$TMP_ROOT/unsafe.conf" <<'EOF'
SAMPLESHEET="$(touch /tmp/atacseq2tracks-unsafe)"
OUTPUT_DIR=/tmp/output
EOF
! "$PYTHON_BIN" "$REPO_DIR/scripts/resolve_config.py" \
    --config "$TMP_ROOT/unsafe.conf" --template "$REPO_DIR/config/config.conf" \
    --shell-output "$TMP_ROOT/unsafe-resolved.conf" --tsv-output "$TMP_ROOT/unsafe.tsv" >/dev/null 2>&1
echo "OK   shell expressions are rejected"

cat > "$TMP_ROOT/duplicate.conf" <<EOF
SAMPLESHEET="$TMP_ROOT/samplesheet.csv"
SAMPLESHEET="$TMP_ROOT/samplesheet.csv"
OUTPUT_DIR="$TMP_ROOT/output"
EOF
! "$PYTHON_BIN" "$REPO_DIR/scripts/resolve_config.py" \
    --config "$TMP_ROOT/duplicate.conf" --template "$REPO_DIR/config/config.conf" \
    --shell-output "$TMP_ROOT/duplicate-resolved.conf" --tsv-output "$TMP_ROOT/duplicate.tsv" >/dev/null 2>&1
echo "OK   duplicate configuration keys are rejected"

bash "$REPO_DIR/atacseq2tracks.sh" --config "$TMP_ROOT/config.conf" --plan >/dev/null
[[ "$(awk -F '\t' 'NR==2 {print $3}' "$TMP_ROOT/output/metadata/technical_merge_audit.tsv")" == "2" ]]
[[ "$(($(wc -l < "$TMP_ROOT/output/metadata/biological_libraries.tsv") - 1))" == "1" ]]
grep -Fq 'before_trim' "$TMP_ROOT/output/metadata/technical_merge_audit.tsv"
echo "OK   technical lanes map to one biological library before trimming"

bash "$REPO_DIR/atacseq2tracks.sh" --help | grep -Fq -- '--preflight-only'
bash "$REPO_DIR/bin/atacseq2tracks" --version | grep -Fxq '4.3.2'
ln -s "$REPO_DIR/bin/atacseq2tracks" "$TMP_ROOT/atacseq2tracks"
"$TMP_ROOT/atacseq2tracks" --version | grep -Fxq '4.3.2'
echo "OK   activation-free launcher resolves an external command symlink"
bash "$REPO_DIR/utilities/regenerate_reports.sh" --help | grep -Fq 'Regenerates reports only'
echo "OK   activation-free launcher and operational interfaces are exposed"

grep -Fq 'cleanup_manifest.tsv' "$REPO_DIR/scripts/cleanup_intermediates.sh"
grep -Fq -- '--from-stage NAME' "$REPO_DIR/atacseq2tracks.sh"
echo "OK   cleanup audit and named recovery are implemented"
