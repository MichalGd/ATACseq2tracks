# 06 — Pipeline steps

[← Running](05_running.md) | [Next: Outputs →](07_outputs.md)

---

## Step 0 — Pre-flight checks

**Script:** `scripts/smoke_test.sh`

Runs before any analysis. Checks:

- Config file exists and sources without error
- All required scripts are present and executable
- Samplesheet passes schema validation (`validate_samplesheet.py`)
- All required tools are in `PATH`
- `picard.jar` and `bedGraphToBigWig` are present
- All reference files (indices, chrom sizes, blacklists, RDS objects) exist
- Required R packages are installed
- At least 50 GB disk space available in `OUTPUT_DIR`

A failure at step 0 stops the pipeline before any files are written.

---

## Step 1 — FastQC (raw reads)

**Script:** `scripts/fastqc_batch.sh`

Runs FastQC on all raw FASTQ files listed in the samplesheet, then generates a MultiQC summary.

- Parallelised: up to `THREADS_PARALLEL_JOBS` samples simultaneously
- Output: `fastQC/fastQC_unTrimmed/` and `multiQC/multiQC_unTrimmed/`

---

## Step 2 — Adapter trimming

**Script:** `scripts/trimgalore_batch.sh`

Trims adapters and low-quality bases using Trim Galore.

- **SE samples:** produces `<stem>_trimmed.fq.gz`
- **PE samples:** produces `<stem>_val_1.fq.gz` and `<stem>_val_2.fq.gz`
- **Technical replicates:** multiple runs with the same `sample_id` + `replicate` are merged before trimming
- Threads: `THREADS_TRIMGALORE` cores per sample
- Output: `trimmedFastq/`

> Requires Trim Galore v2.2.0+ (Oxidized Edition). Older versions use `-o` instead of `--output_dir` and will fail.

---

## Step 3 — FastQC (trimmed reads)

**Script:** `scripts/fastqc_batch.sh`

Same as step 1 but on trimmed FASTQs.

- Output: `fastQC/fastQC_trimmed/` and `multiQC/multiQC_trimmed/`

---

## Step 4 — Alignment

**Script:** `scripts/bowtie2_batch.sh` (dispatches `bowtie2_align.sh`)

Aligns trimmed reads to the reference genome specified per sample in the samplesheet.

- SE: `bowtie2 -x <index> -U <fastq>`
- PE: `bowtie2 -x <index> -1 <R1> -2 <R2>`
- Output sorted, indexed BAM: `bams/<sample_id>.bam`
- Threads: `THREADS_ALIGN` per sample
- MultiQC summary: `multiQC/multiQC_alignments/`

---

## Step 5 — Deduplication

**Script:** `scripts/picard_dedup_batch.sh` (dispatches `picard_dedup.sh`)

Marks and removes PCR duplicates using Picard MarkDuplicates.

Two sub-steps:
1. **Add read group (RG) tags** — required by Picard; uses `sample_id` as SM tag
2. **Mark duplicates** — removes optical and PCR duplicates

Output naming: `dedupBams/<sample_id>_bioR<replicate>_dedup.bam`

> The `_bioR<N>` suffix encodes the biological replicate number from the samplesheet. All downstream scripts look up BAMs using this naming pattern.

- MultiQC summary: `multiQC/multiQC_deduplication/`

---

## Step 6 — Blacklist filtering

**Script:** `scripts/blacklist_filter_batch.sh` (dispatches `blacklist_filter.sh`)

Removes reads overlapping blacklisted genomic regions using `bedtools intersect -v`.

- Blacklist BED path read from samplesheet column 16 (per-sample)
- Input: `dedupBams/<sample_id>_bioR<N>_dedup.bam`
- Output: `filteredBams/<sample_id>_bioR<N>_dedup_blFilt.bam`

The filtered BAMs are the starting point for all downstream steps (7–11).

---

## Step 7 — Genome coverage tracks

**Script:** `scripts/genomecoverage_batch.sh` (dispatches `genomecoverage_single.sh`)

Generates normalised coverage tracks from filtered BAMs.

- bedGraph from `bedtools genomecov`
- Library-size normalisation (RPKM)
- Converted to bigwig with `bedGraphToBigWig`
- Output: `bigwig/<sample_id>_bioR<N>_dedup_blFilt.bw`
- bedGraphs retained in `bedGraph/` and `NormBedGraph/`

---

## Step 8 — Replicate merging and merged tracks

**Script:** `scripts/merge_replicates.sh`

Groups biological replicates by `factor__condition__treatment__cell_type__genome` and generates a merged bigwig track for each group.

- Uses `samtools merge` on filtered BAMs from all replicates in a group
- Generates a single normalised bigwig for the merged BAM
- Output: `bigwig_merged/<group_key>.bw`

---

## Step 9 — Peak calling

**Script:** `scripts/macs2_batch.sh` (dispatches `macs2_peaks.sh`)

Calls peaks for all IP samples using MACS2. Runs in two modes:

**Per-replicate:** one peak set per IP sample
**Pooled:** one peak set per condition group (merged BAMs)

Both always produce **narrow** (`narrowPeak`) **and broad** (`broadPeak`) output, regardless of `macs2_mode`.

- Control BAM is resolved via `control_id` → samplesheet replicate lookup
- Output: `peaks/per_replicate/<sample_id>/narrow/` and `.../broad/`
- Output: `peaks/pooled/<group>/narrow/` and `.../broad/`

---

## Step 10 — ChIPQC

**Script:** `scripts/run_chipqc.R`

Runs the Bioconductor `ChIPQC` package to assess ChIP-seq quality.

- Called twice per genome: once for narrow peaks, once for broad peaks
- Uses pre-built annotation RDS (`CHIPQC_ANNOTATION_HG38/MM39`) and blacklist RDS
- Parallelised via `BiocParallel::MulticoreParam` with `THREADS_CHIPQC` workers
- Outputs per run: HTML report, `chipqc_samplesheet.csv`, `chipqc_metrics_summary.csv`, `chipqc_frip.csv`
- Output directory: `chipqc/ChIPQC_<genome>_<narrow|broad>/`

> Key metric to check: **FRiP** (Fraction of Reads in Peaks). Values > 0.01 are generally acceptable; > 0.05 is good for most marks.

---

## Step 11 — DiffBind samplesheet preparation

**Script:** `scripts/prepare_diffbind.R`

Generates DiffBind-compatible samplesheets from the filtered BAMs and peak files.

- Two output CSVs per genome: one for narrow peaks, one for broad
- Control BAM paths resolved per-sample using the samplesheet `control_id` and `replicate`
- Output: `diffbind/diffbind_samplesheet_<genome>_<narrow|broad>.csv`

See [Downstream: DiffBind](08_diffbind.md) for usage.

---

## Step 12 — UCSC track hub

**Script:** `scripts/create_ucsc_tracks.sh`

Generates a UCSC Genome Browser trackDb text file pointing to the bigwig files.

- Scans `bigwig/` and `bigwig_merged/`
- Output: `reports/ucsc_trackdb.txt`

> You must serve the bigwig files on a web server and configure the base URL in `config.conf` before loading the trackdb into UCSC.

---

## Step 13 — Pipeline report

**Script:** `scripts/generate_pipeline_report.sh`

Generates a summary HTML report of the full pipeline run.

- Collects output counts, MultiQC paths, ChIPQC paths, and DiffBind CSV paths
- Output: `reports/pipeline_report_<YYYYMMDD>.html`

---

[← Running](05_running.md) | [Next: Outputs →](07_outputs.md)
