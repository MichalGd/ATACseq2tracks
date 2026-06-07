# 06 — Pipeline steps

[← Running](05_running.md) | [Next: Outputs →](07_outputs.md)

---

## Step 0 — Pre-flight checks

**Script:** `scripts/smoke_test.sh`

Runs before any analysis. Checks:

- Config file exists and sources without error
- All required scripts are present and executable
- Samplesheet passes schema validation (`validate_samplesheet.py`)
- All required tools are in `PATH` (including `deeptools`, `python3`)
- `picard.jar` and `bedGraphToBigWig` are present
- All reference files (indices, chrom sizes, blacklists) exist
- Required Python packages are importable (`matplotlib`, `pandas`, `numpy`)
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
- Duplication metrics (`*_dup_metrics.txt`) are **reused by Step 10** — no re-deduplication is needed.

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

## Step 10 — Post-alignment QC (deepTools)

**Scripts:** `scripts/post_alignment_qc_batch.sh` + `scripts/plot_chrom_coverage.py`

Runs a comprehensive post-alignment quality control module based on deepTools and standard
command-line tools. This replaces the legacy ChIPQC module.

The module runs in four sequential phases:

**Phase 1 — Per-sample metrics**
- Total aligned reads (`samtools view -c`)
- Mitochondrial read fraction (`samtools idxstats`)
- Duplicate rate (read from pre-computed Picard metrics)
- Peak counts (narrow and broad per replicate)
- FRiP — Fraction of Reads in Peaks (`bedtools intersect`); skipped gracefully if zero peaks

**Phase 2 — Chromosome-level karyogram**
- Per-sample chromosome coverage barcode plots replicating the ChIPQC "ChIP Peaks over Chromosomes" panel
- One row per chromosome in cytogenetic order; X-axis = chromosomal position in bp
- Signal from `bamCoverage --binSize 100000 --normalizeUsing RPKM` (bedtools fallback available)
- Per-sample PNGs + multi-sample panel for direct comparison

**Phase 3 — Genome-wide deepTools QC**
- `plotFingerprint` (Lorenz curve) — IP enrichment vs input
- `multiBamSummary bins` + `plotCorrelation` — Spearman correlation heatmap (10 kb bins)
- `plotPCA` — genome-wide PCA

**Phase 4 — Consensus peaks and peak-centric QC**
- Reproducible peak consensus derived from narrow peaks supported by at least two replicates when possible
- `multiBamSummary BED-file` over consensus peaks + `plotCorrelation` + `plotPCA`
- `computeMatrix reference-point` + `plotHeatmap` + `plotProfile` over consensus peaks
- DESeq2 size factor estimation from consensus peak counts
- Peak-normalised bigwig generation using consensus-block scaling
- FRiP over consensus peak set for all samples

> **Zero-peak safety:** samples with zero peaks are retained in all summary tables and
> flagged `NO_PEAKS`. They are skipped only for peak-centric steps, not for fingerprint,
> bins-level correlation, or PCA.

Config parameter required:
```bash
THREADS_DEEPTOOLS=8   # set to available CPUs / 2 as a reasonable default
```

Output directory: `qc_post_alignment/`
See [Post-alignment QC (deepTools)](12_post_alignment_qc.md) for full documentation.

---

## Step 11 — DiffBind samplesheet preparation

**Script:** `scripts/prepare_diffbind.R`

Generates DiffBind-compatible samplesheets from the filtered BAMs and peak files.

- Two output CSVs per genome: one for narrow peaks, one for broad
- Control BAM paths resolved per-sample using the samplesheet `control_id` and `replicate`
- Preserves optional `batch` metadata if provided in the samplesheet
- Output: `diffbind/diffbind_samplesheet_<genome>_<narrow|broad>.csv`

See [Downstream: DiffBind](08_diffbind.md) for usage.

---

## Step 12 — DiffBind differential analysis

**Scripts:** `scripts/diffbind_analysis.sh`, `scripts/diffbind_analysis.R`

Runs differential accessibility analysis on the prepared DiffBind samplesheets.

- Reads the prefabricated DiffBind CSVs from `diffbind/`
- Counts reads over peaks using `DiffBind::dba.count()`
- Normalises using DiffBind defaults
- Builds contrasts by `Condition`
- Runs `DiffBind::dba.analyze()` and exports results
- Generates PCA, heatmap, MA plot, and volcano plot
- Output: `diffbind_results/<sample_sheet>/`

See [13 — Differential accessibility](13_differential_accessibility.md).

---

## Step 13 — UCSC track hub

**Script:** `scripts/create_ucsc_tracks.sh`

Generates a UCSC Genome Browser trackDb text file pointing to the bigwig files.

- Scans `bigwig/` and `bigwig_merged/`
- Output: `reports/ucsc_trackdb.txt`

> You must serve the bigwig files on a web server and configure the base URL in `config.conf` before loading the trackdb into UCSC.

---

## Step 13 — Pipeline report

**Script:** `scripts/generate_pipeline_report.sh`

Generates a summary HTML report of the full pipeline run.

- Collects output counts, MultiQC paths, QC summary table, and DiffBind CSV paths
- Output: `reports/pipeline_report_<YYYYMMDD>.html`

---

[← Running](05_running.md) | [Next: Outputs →](07_outputs.md)
