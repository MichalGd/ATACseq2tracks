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

**Script:** `scripts/bowtie2_batch.sh`

Aligns trimmed reads to the reference genome specified per sample in the samplesheet.

- SE: `bowtie2 -x <index> -U <fastq>`
- PE: `bowtie2 -x <index> -1 <R1> -2 <R2>`
- Output sorted, indexed BAM: `bams/<sample_id>.bam`
- Threads: `THREADS_ALIGN` per sample
- MultiQC summary: `multiQC/multiQC_alignments/`

---

## Step 5 — Deduplication

**Script:** `scripts/picard_dedup_batch.sh`

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

## Step 6s — Optional Drosophila spike-in calibration

**Scripts:** `scripts/drosophila_spikein_tracks.sh` and
`scripts/filter_composite_spikein_bam.sh`

When `GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS=true`, trimmed reads are
competitively aligned to the configured hg38+dm6 or mm39+dm6 composite index.
The composite BAM is deduplicated, split by species and stringently restricted
to primary canonical MAPQ≥30 alignments with species-matched blacklist removal.
PE fragments are counted once; SE reads are counted once.

Raw stringent host coverage is scaled by
`SPIKEIN_SCALE_TARGET × spikein_to_host_ratio / retained_dm6_count`. This branch
does not use DESeq2 factors. It writes `_SpikeInDM6_Stringent` bigWig/bedGraph
tracks and complete audit tables. See the
[v4.2.0 update note](v4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md).

---

## Step 7 — Genome coverage tracks

**Script:** `scripts/genomecoverage_batch.sh` (dispatches `genomecoverage_single.sh`)

Generates normalised coverage tracks from filtered BAMs.

- Paired-end signal is one properly paired fragment, selected through its first mate and extended across the insert
- Single-end signal uses `SE_SIGNAL_MODE="read"`: one filtered alignment record per observation, without extension or an inferred fragment
- Canonical autosomes plus X/Y only; mitochondrial and noncanonical contigs are excluded
- deepTools CPM with exact scaling and the same fragment/read definition in its denominator
- Outputs: `bigwig/<sample_id>_bioR<N>_dedup_blFilt_CPM.{bw,bedGraph}`
- `GENERATE_CPM_TRACKS=false` disables this coverage family; the two global
  format switches can independently disable bigWig or bedGraph generation

---

## Step 8 — Replicate merging and merged tracks

**Script:** `scripts/merge_replicates.sh`

Groups biological replicates by `factor__condition__treatment__cell_type__genome` and generates a merged bigwig track for each group.

- Uses `samtools merge` on filtered BAMs from all replicates in a group
- Generates CPM bigWig and bedGraph files for the merged BAM with the same layout-aware signal definition
- Outputs: `bigwig_merged/<group_key>_CPM.{bw,bedGraph}`

---

## Step 9 — Peak calling

**Script:** `scripts/macs2_batch.sh` (dispatches `macs2_peaks.sh`)

Calls peaks for all non-control samples using MACS3. Paired-end libraries use fragment-aware `BAMPE`; single-end ATAC libraries use the configured shift and extension.

The SE MACS3 shift/extension affects peak calling only. Normalization tracks and consensus counts remain read-based and do not use MACS3 pseudo-fragment extension.

**Per-replicate:** one peak set per IP sample
**Pooled:** one peak set per condition group (merged BAMs)

The samplesheet `macs2_mode` selects narrow, broad, both, or none. The legacy column and script names are retained for compatibility, but the executable is consistently `macs3`.

- Control BAM is resolved via `control_id` → samplesheet replicate lookup
- Output: `peaks/per_replicate/<sample_id>/narrow/` and `.../broad/`
- Output: `peaks/pooled/<group>/narrow/` and `.../broad/`

---

## Step 10 — Post-alignment QC (deepTools)

**Scripts:** `scripts/post_alignment_qc_batch.sh`, `scripts/plot_chrom_coverage.py`, and `scripts/ataqv_qc_batch.sh`

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
- Signal from `bamCoverage --binSize 100000 --normalizeUsing RPKM` with the shared canonical and fragment/read semantics
- Per-sample PNGs + multi-sample panel for direct comparison

**Phase 3 — Genome-wide deepTools QC**
- `plotFingerprint` (Lorenz curve) — IP enrichment vs input
- `multiBamSummary bins` + `plotCorrelation` — Spearman correlation heatmap (10 kb bins)
- `plotPCA` — genome-wide PCA

**Phase 4 — Consensus peaks and peak-centric QC**
- Consensus peaks derived from the selected peak type and supported by at least `CONSENSUS_MIN_SAMPLES` distinct biological sample keys (default two; no silent one-sample fallback)
- `multiBamSummary BED-file` over consensus peaks + `plotCorrelation` + `plotPCA`
- `computeMatrix reference-point` + `plotHeatmap` + `plotProfile` over consensus peaks
- Canonical-chromosome filtering of peak inputs before consensus construction
- Fragment/read-aware consensus-peak counting using the same semantics as coverage tracks
- DESeq2 size factor estimation from consensus peak counts
- DESeq2-consensus bigWig/bedGraph generation using inverse size factors, without an additional CPM divisor
- DESeq2 robust-CPM bigWig/bedGraph generation using `1e6 / (size_factor * geometric_mean(column_sums))`
- Per-sample normalization metadata, including the CPM count, column sum, cohort constant and both scales
- FRiP over consensus peak set for all samples

**ATAC-specific QC**
- ENCODE-style TSS enrichment and short-to-mononucleosomal ratio from `ataqv`
- Compressed ataqv JSON and optional local viewer
- Paired-end fragment-size bins, nucleosome-periodicity metrics, and PNG/PDF plots
- Explicit not-applicable status for single-end fragment periodicity

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

## Step 10b — Read-filtering sensitivity coverage

**Scripts:** `scripts/generate_filtering_sensitivity_tracks.sh` and
`scripts/filter_bam_for_coverage_policy.sh`

This independently checkpointed stage uses the fixed Step 10 consensus BED.
The permissive branch filters pre-dedup Bowtie2 BAMs at MAPQ 0; the intermediate
branch filters Picard-deduplicated BAMs at MAPQ 0. Both retain only primary,
non-supplementary records, apply identical canonical-chromosome and blacklist
rules, and count one proper fragment for PE data or one read for SE data.

Each enabled policy obtains a separate consensus count matrix, DESeq2
`poscounts` size factors and robust cohort constant. The stringent family uses
the existing deduplicated `MIN_MAPQ=30` final BAM and Step 10 factor table. The
permissive and intermediate tracks are sensitivity outputs, not replacements
for the stringent/current production track.

Outputs include policy tracks under `bigwig_deseq2_robust_cpm/<policy>/`,
matrices and factors under `coverage_filtering_sensitivity/<policy>/`, and a
combined filtering/scale-factor metadata table. Temporary policy BAMs are
removed after full workflow success unless
`KEEP_NORMALIZATION_POLICY_BAMS=true`.

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
- Builds a `minOverlap=2` all-sample consensus with
  `DiffBind::dba.count(summits=DIFFBIND_SUMMITS)`; the ATAC default is a
  100-bp half-width (approximately 201-bp windows)
- Restricts the resulting universe to canonical nuclear chromosomes and removes
  intervals overlapping the configured genome-matched blacklist
- Keeps every sample in consensus construction, then recounts model-eligible
  samples on that fixed universe without a second summit recentering
- Normalises peak counts through DiffBind and explicitly uses the DESeq2 analysis method
- Fits one multi-condition model and exports every pair of conditions having at
  least two biological samples; singleton conditions are excluded only here
- Runs `DiffBind::dba.analyze(method = DBA_DESEQ2)` and exports DESeq2 results
- Generates PCA, heatmap, MA plot, and volcano plot
- Output: `diffbind_results/<sample_sheet>/`

### Peer Step 12a - DESeq2ATAC

**Scripts:** `scripts/deseq2atac_analysis.sh`, `scripts/deseq2atac_analysis.R`

Runs two independent raw-count DESeq2 analyses over replicate-supported MACS3
peak regions: one broad-peak analysis and one narrow-peak analysis.

- Requires both peak types for every biological sample (`macs2_mode=both`)
- Constructs the broad and narrow consensuses independently
- For each peak type, disjoins sample peaks into atomic segments, retains atoms
  with support >=2 by default (`DESEQ2ATAC_MIN_SAMPLES`), then reduces retained
  atoms to sorted non-overlapping regions
- Applies the shared canonical chromosome and blacklist universe before support
  is calculated
- Counts each PE proper pair once or each SE alignment once
- Uses raw integer counts and DESeq2 `poscounts` size factors
- Keeps all samples in consensus construction/raw descriptive counts, fits one
  model to conditions with at least two biological samples, and exports all pairs
- Supports an explicitly configured, full-rank block term
- Adds GTF gene context/nearest TSS and, by default, required genome-matched
  cCRE annotation; `RUN_CCRE_ANNOTATION=false` selects GTF-only annotation
- Exports compressed matrices/results and publication-quality PNG/PDF diagnostics
- Treats zero FDR hits as a successful, explicitly documented result
- Uses `.checkpoints/step12a.done`, independently of DiffBind Step 12
- Output: `deseq2atac/broad/`, `deseq2atac/narrow/`, and a peak-type summary

Both differential modules and both peak types are attempted before a module
failure is propagated. Pair-level TSV/HTML reports are generated and cleanup is
suppressed before a failed run exits.

See [13 — Differential accessibility](13_differential_accessibility.md).
See [16 — Peak annotation](16_peak_annotation.md) for the annotation columns,
category precedence, cCRE sources and interpretation limits.

---

## Step 13 — UCSC track definitions

**Script:** `scripts/create_ucsc_tracks.sh`

Generates UCSC Genome Browser custom-track definitions for every bigWig family.

- Processes `bigwig/`, `bigwig_deseq2_consensus/`, the legacy robust directory,
  each `bigwig_deseq2_robust_cpm/<policy>/` directory, and `bigwig_merged/`
  whenever that directory contains bigWigs
- Output: `ucsc_tracks.txt` inside each processed bigWig directory

> You must serve the bigwig files on a web server and configure the base URL in `config.conf` before loading the trackdb into UCSC.

---

## Step 14 — Pipeline report

**Script:** `scripts/generate_pipeline_report.sh`

Generates a summary HTML report of the full pipeline run.

- Collects output counts, MultiQC paths, QC summaries, and all pair-level
  differential results/statuses
- Writes `differential_accessibility_summary.tsv` and `.html`
- Output: `reports/pipeline_report_<YYYYMMDD>.html`

---

## Post-success cleanup

Cleanup is not a checkpointed analysis step. With the default
`ENABLE_AUTOMATIC_CLEANUP=true`, it runs only after Steps 13 and 14 complete and
all enabled differential modules have returned successfully. The default `KEEP_*`
policy deletes pre-dedup BAMs, trimmed FASTQs, pre-blacklist deduplicated BAMs,
temporary permissive/intermediate policy BAMs and raw bedGraphs, while retaining filtered quantitative BAMs and every final
track, peak, QC, differential result and report. A failed differential module
still permits reporting, then causes a non-zero exit before cleanup. Set
`ENABLE_AUTOMATIC_CLEANUP=false` to retain all intermediates.

---

[← Running](05_running.md) | [Next: Outputs →](07_outputs.md)
