#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/atacseq2tracks-spikein.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT
REAL_PYTHON="$(command -v python3 2>/dev/null || true)"
HAVE_REAL_PYTHON=false
if [[ -n "$REAL_PYTHON" ]] && "$REAL_PYTHON" -c 'import csv' >/dev/null 2>&1; then
    HAVE_REAL_PYTHON=true
fi

# The pure scaling helper must implement C * declared_ratio / retained_dm6.
# shellcheck disable=SC1091
source "${REPO_DIR}/scripts/track_normalization_helpers.sh"
[[ "$(spikein_scale_factor 1000000 1 2000)" == "500" ]]
[[ "$(spikein_scale_factor 1000000 0.5 2000)" == "250" ]]
! spikein_scale_factor 1000000 1 0 >/dev/null 2>&1
echo "OK   dm6 spike-in scaling formula"

grep -Fq 'host.header_filtered.bam' "${REPO_DIR}/scripts/filter_composite_spikein_bam.sh"
grep -Fq 'sub(/^SN:dm6__/,"SN:",$i)' "${REPO_DIR}/scripts/filter_composite_spikein_bam.sh"
echo "OK   species-split BAM headers are restricted to their canonical universe"

mkdir -p "$TMP_ROOT/bin" "$TMP_ROOT/index" "$TMP_ROOT/trimmed" \
    "$TMP_ROOT/output/spikein/composite_bams" \
    "$TMP_ROOT/output/spikein/dedup_bams" \
    "$TMP_ROOT/output/spikein/stringent_host_bams" \
    "$TMP_ROOT/output/spikein/stringent_dm6_bams"

for component in 1 2 3 4; do printf x > "$TMP_ROOT/index/hg38_dm6.${component}.bt2"; done
for component in 1 2; do printf x > "$TMP_ROOT/index/hg38_dm6.rev.${component}.bt2"; done
printf 'chr1\t0\t1\n' > "$TMP_ROOT/host.blacklist.bed"
printf 'chr2L\t0\t1\n' > "$TMP_ROOT/dm6.blacklist.bed"
printf 'chr2L\t23513712\n' > "$TMP_ROOT/dm6.chrom.sizes"

cat > "$TMP_ROOT/samplesheet.csv" <<EOF
sample_id,replicate,layout,genome,is_control,blacklist,spikein_genome,spikein_to_host_ratio,spikein_stage
sampleA,1,PE,hg38,FALSE,$TMP_ROOT/host.blacklist.bed,dm6,1,pre_tagmentation_nuclei
EOF

cat > "$TMP_ROOT/config.conf" <<EOF
INDEX_HG38_DM6="$TMP_ROOT/index/hg38_dm6"
INDEX_MM39_DM6="$TMP_ROOT/index/mm39_dm6"
BLACKLIST_DM6="$TMP_ROOT/dm6.blacklist.bed"
CHROM_SIZES_DM6="$TMP_ROOT/dm6.chrom.sizes"
GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS=true
GENERATE_DROSOPHILA_CONTROL_TRACKS=true
GENERATE_COVERAGE_BIGWIGS=true
GENERATE_COVERAGE_BEDGRAPHS=true
SPIKEIN_MIN_MAPQ=30
SPIKEIN_SCALE_TARGET=1000000
SPIKEIN_MIN_FRAGMENTS_FAIL=1000
SPIKEIN_MIN_FRAGMENTS_WARN=10000
SPIKEIN_WARN_LOW_FRACTION=0.001
SPIKEIN_WARN_HIGH_FRACTION=0.20
SPIKEIN_PARALLEL_JOBS=1
THREADS_ALIGN=1
THREADS_SAMTOOLS=1
THREADS_BIGWIG=1
TRACK_BIN_SIZE=10
SE_SIGNAL_MODE=read
BAMCOVERAGE_COMMON_ARGS=""
EOF

key=sampleA_bioR1
for bam in \
    "$TMP_ROOT/output/spikein/composite_bams/${key}_host_dm6.bam" \
    "$TMP_ROOT/output/spikein/dedup_bams/${key}_host_dm6_dedup.bam" \
    "$TMP_ROOT/output/spikein/stringent_host_bams/${key}_host_stringent.bam" \
    "$TMP_ROOT/output/spikein/stringent_dm6_bams/${key}_dm6_stringent.bam"; do
    printf mock > "$bam"
done

cat > "$TMP_ROOT/bin/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  quickcheck) exit 0 ;;
  --version) echo 'samtools 1.mock' ;;
  view)
    if [[ " $* " == *" -c "* ]]; then
      bam="${!#}"
      if [[ "$bam" == *dm6_stringent.bam ]]; then echo 2000
      elif [[ "$bam" == *host_stringent.bam ]]; then echo 8000
      elif [[ " $* " == *" -q 30 "* ]]; then echo 9000
      else echo 10000
      fi
    else
      exit 1
    fi
    ;;
  *) echo "unexpected mocked samtools call: $*" >&2; exit 1 ;;
esac
EOF
cat > "$TMP_ROOT/bin/bamCoverage" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
output=''
scale=''
while (( $# )); do
  case "$1" in
    --outFileName) output="$2"; shift 2 ;;
    --scaleFactor) scale="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ "$scale" == "500" || "$scale" == "1" ]] \
    || { echo "wrong mocked scale: $scale" >&2; exit 1; }
printf 'mock-track\n' > "$output"
EOF
cat > "$TMP_ROOT/bin/bowtie2" <<'EOF'
#!/usr/bin/env bash
[[ "${1:-}" == "--version" ]] && { echo 'bowtie2 2.mock'; exit 0; }
exit 1
EOF
cat > "$TMP_ROOT/bin/python3" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf 'key\tsample_id\treplicate\tlayout\thost_genome\thost_blacklist\tspikein_genome\tspikein_to_host_ratio\tspikein_stage\tis_control\n'
printf 'sampleA_bioR1\tsampleA\t1\tPE\thg38\t%s\tdm6\t1\tpre_tagmentation_nuclei\tfalse\n' "$TEST_HOST_BLACKLIST"
EOF
chmod +x "$TMP_ROOT/bin/samtools" "$TMP_ROOT/bin/bamCoverage" \
    "$TMP_ROOT/bin/bowtie2" "$TMP_ROOT/bin/python3"

PATH="$TMP_ROOT/bin:$PATH" F2T_CONFIG="$TMP_ROOT/config.conf" \
    TEST_HOST_BLACKLIST="$TMP_ROOT/host.blacklist.bed" \
    bash "${REPO_DIR}/scripts/drosophila_spikein_tracks.sh" \
    "$TMP_ROOT/samplesheet.csv" "$TMP_ROOT/trimmed" "$TMP_ROOT/output"

BW="$TMP_ROOT/output/bigwig_spikein/stringent/${key}_SpikeInDM6_Stringent.bw"
BG="$TMP_ROOT/output/bigwig_spikein/stringent/${key}_SpikeInDM6_Stringent.bedGraph"
TABLE="$TMP_ROOT/output/spikein/tables/spikein_normalization.tsv"
DM6_RAW_BW="$TMP_ROOT/output/bigwig_spikein/dm6_control/${key}_dm6_StringentRaw.bw"
DM6_RAW_BG="$TMP_ROOT/output/bigwig_spikein/dm6_control/${key}_dm6_StringentRaw.bedGraph"
DM6_CPM_BW="$TMP_ROOT/output/bigwig_spikein/dm6_control/${key}_dm6_StringentCPM.bw"
DM6_CPM_BG="$TMP_ROOT/output/bigwig_spikein/dm6_control/${key}_dm6_StringentCPM.bedGraph"
[[ -s "$BW" && -s "$BG" && -s "$TABLE" ]]
[[ -s "$DM6_RAW_BW" && -s "$DM6_RAW_BG" && -s "$DM6_CPM_BW" && -s "$DM6_CPM_BG" ]]
awk -F '\t' 'NR==2 {exit !($12==2000 && $18==1000000 && $19==500 && $20==1 && $21==500)}' "$TABLE"
grep -Fq $'sampleA_bioR1\tlow_dm6_count\t2000\t10000' \
    "$TMP_ROOT/output/spikein/tables/spikein_warnings.tsv"
echo "OK   mocked stringent spike-in tracks, metadata and warning"
echo "OK   raw and CPM dm6 UCSC control-track outputs"

grep -Fq 'ucsc_tracks_dm6.txt' "${REPO_DIR}/atacseq2tracks.sh"
grep -Fq '*_dm6_StringentRaw' "${REPO_DIR}/scripts/create_ucsc_tracks.sh"
grep -Fq '*_dm6_StringentCPM' "${REPO_DIR}/scripts/create_ucsc_tracks.sh"
echo "OK   dm6-specific UCSC track-definition contract"

# Legacy sheets remain valid when spike-in mode is not requested, while
# partial declarations are rejected.
if [[ "$HAVE_REAL_PYTHON" == "true" ]]; then
    cat > "$TMP_ROOT/partial.csv" <<EOF
sample_id,replicate,layout,genome,is_control,blacklist,spikein_genome
sampleA,1,PE,hg38,FALSE,$TMP_ROOT/host.blacklist.bed,dm6
EOF
    ! "$REAL_PYTHON" "${REPO_DIR}/scripts/validate_samplesheet.py" \
        "$TMP_ROOT/partial.csv" >/dev/null 2>&1
    echo "OK   incomplete spike-in metadata is rejected"
else
    echo "WARN Python runtime unavailable; validator execution deferred to Linux smoke test"
fi

echo "OK   v4.2.0 Drosophila spike-in regression checks"
