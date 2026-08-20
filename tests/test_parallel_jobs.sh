#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/scripts/parallel_job_helpers.sh"

fail() { echo "FAIL $*" >&2; exit 1; }
ok() { echo "OK   $*"; }

parallel_require_positive_integer TEST 4 || fail "positive integer rejected"
if parallel_require_positive_integer TEST 0 2>/dev/null; then
    fail "zero job limit accepted"
fi
ok "parallel limits require positive integers"

parallel_pool_init 2
parallel_pool_submit pass bash -c 'exit 0'
parallel_pool_submit fail_job bash -c 'exit 7'
if parallel_pool_wait_all; then
    fail "failed child process did not fail the pool"
fi
[[ "$PARALLEL_POOL_FAILURES" -eq 1 ]] || fail "unexpected failure count"
[[ "$(parallel_failed_labels_csv)" == "fail_job" ]] || fail "failed label was not retained"
ok "child failures propagate after all submitted jobs are collected"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
ordered_worker() {
    local index="$1" delay="$2"
    sleep "$delay"
    printf '%s\n' "$index" > "${tmp}/${index}.row"
}
parallel_pool_init 3
parallel_pool_submit sample_0 ordered_worker 0 0.20
parallel_pool_submit sample_1 ordered_worker 1 0.05
parallel_pool_submit sample_2 ordered_worker 2 0.10
parallel_pool_wait_all || fail "ordered workers failed"
for index in 0 1 2; do cat "${tmp}/${index}.row"; done > "${tmp}/merged"
[[ "$(tr '\n' ' ' < "${tmp}/merged")" == "0 1 2 " ]] || fail "deterministic parent merge failed"
ok "worker completion order does not alter merged sample order"

parallel_write_timing_row "${tmp}/timing" test sample 100 109 SUCCESS 4 8
[[ "$(cut -f5 "${tmp}/timing")" == 9 ]] || fail "elapsed timing is incorrect"
ok "per-job timing rows include elapsed seconds"

grep -q 'QC_SAMPLE_PARALLEL_JOBS:-4' "${ROOT}/scripts/post_alignment_qc_batch.sh" || fail "QC default missing"
grep -q 'TRACK_PARALLEL_JOBS:-2' "${ROOT}/scripts/post_alignment_qc_batch.sh" || fail "track default missing"
grep -q 'ATAQV_PARALLEL_JOBS:-4' "${ROOT}/scripts/ataqv_qc_batch.sh" || fail "ataqv default missing"
grep -q 'POOLED_MACS_PARALLEL_JOBS:-2' "${ROOT}/scripts/macs2_batch.sh" || fail "pooled MACS default missing"
grep -q 'MERGE_PARALLEL_JOBS:-2' "${ROOT}/scripts/merge_replicates.sh" || fail "merge default missing"
grep -q 'COVERAGE_FILTER_PARALLEL_JOBS:-' "${ROOT}/scripts/generate_filtering_sensitivity_tracks.sh" || fail "coverage-filter default missing"
ok "all six configurable defaults are wired into their batch stages"

mock_bin="${tmp}/bin"
mkdir -p "$mock_bin" "${tmp}/bams" "${tmp}/merge_out" "${tmp}/peaks"
printf '#!/usr/bin/env bash\nset -euo pipefail\ncmd="$1"; shift\ncase "$cmd" in\nmerge) out=""; while (( $# )); do case "$1" in -@) shift 2;; -f) out="$2"; shift 2;; *) shift;; esac; done; printf "BAM\\n" > "$out";;\nindex) bam="${@: -1}"; printf "BAI\\n" > "${bam}.bai";;\nquickcheck) exit 0;;\nidxstats) printf "chr1\\t1000000\\t100\\t0\\n*\\t0\\t0\\t0\\n";;\nview) out=""; count=false; while (( $# )); do case "$1" in -c) count=true; shift;; -o) out="$2"; shift 2;; -@|-f|-F) shift 2;; -b) shift;; *) shift;; esac; done; if $count; then echo 100; elif [[ -n "$out" ]]; then printf "BAM\\n" > "$out"; fi;;\nesac\n' > "${mock_bin}/samtools"
printf '#!/usr/bin/env bash\nset -euo pipefail\nout=""\nwhile (( $# )); do case "$1" in --outFileName|-o) out="$2"; shift 2;; *) shift;; esac; done\nprintf "track\\n" > "$out"\n' > "${mock_bin}/bamCoverage"
printf '#!/usr/bin/env bash\nset -euo pipefail\noutdir=""; name=""; broad=false\nwhile (( $# )); do case "$1" in --outdir) outdir="$2"; shift 2;; -n) name="$2"; shift 2;; --broad) broad=true; shift;; *) shift;; esac; done\nmkdir -p "$outdir"\nif $broad; then suffix=broadPeak; else suffix=narrowPeak; fi\nprintf "chr1\\t0\\t100\\n" > "${outdir}/${name}_peaks.${suffix}"\n' > "${mock_bin}/macs3"
chmod +x "${mock_bin}/samtools" "${mock_bin}/bamCoverage" "${mock_bin}/macs3"

samplesheet="${tmp}/samplesheet.csv"
printf '%s\n' \
  'sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix' \
  'sampleA1,a_R1.fastq.gz,a_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,A,none,cell,1,1,FALSE,,both,,sampleA1' \
  'sampleA2,b_R1.fastq.gz,b_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,A,none,cell,2,1,FALSE,,both,,sampleA2' \
  'sampleB1,c_R1.fastq.gz,c_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,B,none,cell,1,1,FALSE,,both,,sampleB1' \
  'sampleB2,d_R1.fastq.gz,d_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,B,none,cell,2,1,FALSE,,both,,sampleB2' \
  > "$samplesheet"
for key in sampleA1_bioR1 sampleA2_bioR2 sampleB1_bioR1 sampleB2_bioR2; do
    printf 'BAM\n' > "${tmp}/bams/${key}_dedup_blFilt.bam"
done
config="${tmp}/config.conf"
printf 'MERGE_PARALLEL_JOBS=2\nPOOLED_MACS_PARALLEL_JOBS=2\nTHREADS_PARALLEL_JOBS=2\nTHREADS_SAMTOOLS=1\nTHREADS_BIGWIG=1\nTRACK_STANDARD_CHROMS_ONLY=true\nSE_SIGNAL_MODE=read\nMACS3_COMMAND="%s/macs3"\n' \
    "$mock_bin" > "$config"
export F2T_CONFIG="$config"
PATH="${mock_bin}:$PATH" bash "${ROOT}/scripts/merge_replicates.sh" \
    "$samplesheet" "${tmp}/bams" "${tmp}/merge_out" hg38 >/dev/null
[[ "$(find "${tmp}/merge_out" -maxdepth 1 -name '*_CPM.bw' | wc -l | tr -d ' ')" == 2 ]] \
    || fail "parallel merge did not create both group tracks"
[[ "$(tail -n +2 "${tmp}/merge_out/merge_job_timing.tsv" | wc -l | tr -d ' ')" == 2 ]] \
    || fail "merge timing table does not contain one row per group"
ok "mocked replicate groups run through the bounded merge pool"

PATH="${mock_bin}:$PATH" bash "${ROOT}/scripts/macs2_batch.sh" \
    "$samplesheet" "${tmp}/bams" "${tmp}/peaks" >/dev/null
[[ "$(find "${tmp}/peaks/per_replicate" -name '*_peaks.narrowPeak' | wc -l | tr -d ' ')" == 4 ]] \
    || fail "per-replicate MACS jobs are incomplete"
[[ "$(find "${tmp}/peaks/pooled" -name '*_peaks.broadPeak' | wc -l | tr -d ' ')" == 2 ]] \
    || fail "pooled MACS jobs are incomplete"
[[ "$(tail -n +2 "${tmp}/peaks/macs3_job_timing.tsv" | wc -l | tr -d ' ')" == 6 ]] \
    || fail "MACS timing table does not contain four replicate and two pooled rows"
ok "mocked replicate and pooled MACS3 jobs complete through bounded pools"
