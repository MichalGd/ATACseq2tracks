# 12 — Post-alignment QC (deepTools)

[← Blacklist filtering](11_blacklist_filtering.md) | [← Pipeline steps](06_pipeline_steps.md)

---

## Overview

In v4.0.0, Step 10 combines the assay-general deepTools module with ATAC-specific `ataqv` QC. For non-control ATAC-seq samples, `ataqv` calculates ENCODE-style TSS enrichment and the short-to-mononucleosomal ratio. Paired-end samples additionally receive compact full-scan fragment histograms, nucleosome-free/mono/di/tri fractions, NFR-to-mono ratio, local peak spacing, and 300-dpi PNG plus vector PDF periodicity plots. Single-end fragment periodicity is reported as not applicable.

Step 10 of the ATACseq2tracks pipeline runs a post-alignment QC module based on
[deepTools](https://deeptools.readthedocs.io/) and standard command-line tools.
It replaces the legacy ChIPQC (Bioconductor) module and produces robust,
publication-ready QC outputs that work across all supported assay types, including
samples with very few or zero peaks.

**Scripts:**
- `scripts/post_alignment_qc_batch.sh` — main orchestrator
- `scripts/plot_chrom_coverage.py` — chromosome-level karyogram plots

**Inputs:**
- Blacklist-filtered BAMs from `filteredBams/`
- MACS2 peak files from `peaks/`
- Normalised bigWig files from `bigwig/`
- Picard deduplication metrics from `dedupBams/`
- Samplesheet (IP samples only; controls excluded from multi-sample plots)

**Output directory:** `qc_post_alignment/`

ATAC-specific scripts are `ataqv_qc_batch.sh`, `prepare_tss_bed.py`, `extract_ataqv_metrics.py`, and `plot_fragment_periodicity.py`. Their outputs are written under `qc_post_alignment/atac_qc/`. Detailed ataqv JSON is compressed; fragment histograms contain one row per length bin rather than one row per fragment.

## Bounded parallel execution

Per-sample metrics, chromosome-coverage QC and consensus FRiP use
`QC_SAMPLE_PARALLEL_JOBS=4` by default. DESeq2 consensus/robust track generation
uses `TRACK_PARALLEL_JOBS=2`. ATAC-specific ataqv plus periodicity processing
uses `ATAQV_PARALLEL_JOBS=4`; `mkarv` still runs once after all sample JSON files
have completed. Each worker writes sample-specific temporary rows, and the
parent process merges tables in samplesheet order. A failed required child job
makes its stage fail after the other submitted jobs have been collected.

Timing records are written to:

- `qc_post_alignment/tables/parallel_job_timing.tsv`;
- `qc_post_alignment/atac_qc/tables/ataqv_job_timing.tsv`.

The thread settings remain per worker. In particular, maximum ataqv CPU demand
is approximately `ATAQV_PARALLEL_JOBS × THREADS_ATAQV`, and track generation is
approximately `TRACK_PARALLEL_JOBS × THREADS_BIGWIG`. Lower the job count first
if shared storage becomes saturated.

---

## Supported assay types

| Assay | Notes |
|---|---|
| ChIP-seq | Full support. All metrics applicable. |
| ChIPmentation | Full support. Expect higher duplication rates. |
| CUT&RUN | Full support. Typically narrow peaks; FRiP often >0.3. |
| CUT&TAG | Full support. Similar to CUT&RUN. |
| ATAC-seq | Full support. Use merged broad+narrow consensus peaks. |

---

## Processing phases

The module runs in four sequential phases.

---

### Phase 1 — Per-sample metrics

Collects a per-sample QC summary from existing pipeline outputs. No re-alignment
or re-deduplication is performed.

| Metric | Source | Description |
|---|---|---|
| Total aligned reads | `samtools view -c` | Total reads in filtered BAM |
| Mitochondrial reads / % | `samtools idxstats` chrM | High chrM suggests poor IP or cell lysis |
| Duplicate rate | Picard `*_dup_metrics.txt` | Read from pre-computed Step 5 metrics |
| Number of narrow peaks | MACS2 `narrowPeak` | Line count of per-replicate narrowPeak file |
| Number of broad peaks | MACS2 `broadPeak` | Line count of per-replicate broadPeak file |
| FRiP (narrow) | `bedtools intersect` | Fraction of Reads in Peaks — narrow peak set |
| FRiP (broad) | `bedtools intersect` | Fraction of Reads in Peaks — broad peak set |

> **Zero-peak safety:** FRiP is only calculated when >=1 peak exists. Samples with
> zero peaks are retained in the summary table and flagged `NO_PEAKS` — they are
> never silently dropped.

Output: `qc_post_alignment/tables/qc_summary.tsv`

---

### Phase 2 — Chromosome-level karyogram plots

**Script:** `scripts/plot_chrom_coverage.py`

Generates per-sample coverage plots showing signal distribution across all
chromosomes — replicating the ChIPQC "ChIP Peaks over Chromosomes" panel style:

- One row per canonical nuclear chromosome (autosomes, chrX and chrY) in cytogenetic order
- X-axis: chromosomal position in bp, shared scale across all chromosomes
- Signal bars: 100 kb bin RPKM from `bamCoverage`, using fragment-based PE or read-based SE signal
- Grey background bar: chromosome length reference
- Chromosome labels on the right margin

Two outputs are produced:
- `<KEY>_karyogram.png` — one plot per sample
- `karyogram_all_samples.png` — multi-panel grid with all IP samples together

**Tools used:**
- `deeptools bamCoverage --binSize 100000 --normalizeUsing RPKM`
- `matplotlib` / `pandas` for plotting

A read-based bedtools fallback is deliberately not used, because it would violate the paired-end fragment signal definition.

Output directory: `qc_post_alignment/plots/chromosome_coverage/`
Per-chromosome table: `qc_post_alignment/tables/per_chromosome_reads.tsv`

---

### Phase 3 — Genome-wide deepTools QC

Uses all IP sample BAMs together for multi-sample comparisons.

#### Fingerprint plot

```bash
deeptools plotFingerprint \
    --bamfiles <all IP BAMs> \
    --plotFile fingerprint.png \
    --outRawCounts fingerprint_metrics.tsv
```

A fingerprint (Lorenz curve) shows how reads are distributed across the genome.
A flat diagonal = uniform coverage (typical of input/control). A steep curve
toward the top-right = reads concentrated at enriched regions (good IP).

Output: `qc_post_alignment/plots/fingerprint.png`
Metrics: `qc_post_alignment/tables/fingerprint_metrics.tsv`

#### Sample correlation (Spearman, genome-wide bins)

```bash
deeptools multiBamSummary bins \
    --bamfiles <all IP BAMs> \
    --binSize 10000 \
    --outFileName bins_summary.npz

deeptools plotCorrelation \
    --corData bins_summary.npz \
    --corMethod spearman \
    --whatToPlot heatmap \
    --plotFile correlation_heatmap_bins.png
```

Genome-wide correlation across 10 kb bins. Biological replicates of the same
condition should show Spearman r >= 0.90. Different conditions or factors
should cluster separately.

Output: `qc_post_alignment/plots/correlation_heatmap_bins.png`

#### PCA (genome-wide bins)

```bash
deeptools plotPCA \
    --corData bins_summary.npz \
    --plotFile pca_bins.png \
    --outFileNameData pca_bins.tsv
```

PCA of genome-wide signal. Replicates of the same condition cluster together.
Outlier samples are immediately visible.

Output: `qc_post_alignment/plots/pca_bins.png`
Data: `qc_post_alignment/matrices/pca_bins.tsv`

---

### Phase 4 — Consensus peaks + peak-centric QC

Builds a reproducible consensus peak set from narrow peaks and performs peak-centric
multi-sample comparisons.

#### Consensus peak set construction

When replicate narrow peak calls are available, the pipeline builds a reproducible
consensus peak set from regions observed in at least two per-replicate narrow peaks.
If no reproducible narrow consensus can be computed, it falls back to the merged
union of narrow and broad peaks.

```bash
cat peaks/per_replicate/*/narrow/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > merged_narrow.bed
cat peaks/per_replicate/*/broad/*.broadPeak  | sort -k1,1 -k2,2n | bedtools merge > merged_broad.bed
bedtools multiinter -i peaks/per_replicate/*/narrow/*.narrowPeak \
    | awk '$4 >= 2 {print $1"\t"$2"\t"$3}' \
    | bedtools merge -i stdin > consensus_peaks_supported.bed
# Consensus = reproducible narrow peaks when replicates exist
cp consensus_peaks_supported.bed consensus_peaks.bed || true
if [[ ! -s consensus_peaks.bed ]]; then
  cat merged_narrow.bed merged_broad.bed | sort -k1,1 -k2,2n | bedtools merge > consensus_peaks.bed
fi
```

Samples with zero peaks still contribute signal via their BAM — they simply
do not add regions to the consensus set.

Output: `qc_post_alignment/peak_sets/consensus_peaks.bed`

Additional consensus outputs:
- `qc_post_alignment/matrices/consensus_peak_counts.tsv` — raw consensus peak count matrix
- `qc_post_alignment/matrices/consensus_peak_normCounts.tsv` — DESeq2 normalized consensus peak counts
- `qc_post_alignment/tables/consensus_sizeFactors.tsv` — DESeq2 size factors, consensus column sums, cohort geometric mean and robust CPM scales
- `qc_post_alignment/tables/track_normalization_metadata.tsv` — per-sample fragment/read count used for CPM plus every applied track scale

All peak inputs are first restricted to the same canonical autosome/X/Y universe used for tracks. Paired-end peak counts use one properly paired first-mate record per fragment and extend it across the insert; single-end counts remain read-based.

Two DESeq2-derived track families are generated. `*_DESeq2Consensus` multiplies raw coverage by `1 / s_j`. `*_DESeq2RobustCPM` multiplies it by `1e6 / (s_j * G)`, where `G = exp(mean(log(colSums(K))))`. This matches `DESeq2::fpm(robust=TRUE)` on the consensus count matrix. Because `1e6/G` is common to the cohort, robust CPM is only a cohort-wide unit conversion of the consensus track, not an additional sample-specific normalization.

#### Signal correlation and PCA over consensus peaks

```bash
deeptools multiBamSummary BED-file \
    --BED consensus_peaks.bed \
    --bamfiles <all IP BAMs> \
    --outFileName peaks_summary.npz

deeptools plotCorrelation ... --plotFile correlation_heatmap_peaks.png
deeptools plotPCA         ... --plotFile pca_peaks.png
```

Peak-centric correlation is more stringent than genome-wide bins because
background noise regions are excluded. Replicates should show r >= 0.90
over consensus peaks.

Outputs:
- `qc_post_alignment/plots/correlation_heatmap_peaks.png`
- `qc_post_alignment/plots/pca_peaks.png`

#### Signal heatmap and average profile over consensus peaks

```bash
deeptools computeMatrix reference-point \
    --referencePoint center \
    --regionsFileName consensus_peaks.bed \
    --scoreFileName <all IP bigWigs> \
    --upstream 2000 --downstream 2000 \
    --outFileName matrix_peaks.gz

deeptools plotHeatmap --matrixFile matrix_peaks.gz --plotFile heatmap_peaks.png
deeptools plotProfile --matrixFile matrix_peaks.gz --plotFile profile_peaks.png
```

The heatmap shows per-peak signal intensity for all samples, sorted by mean
signal. The profile shows average enrichment ±2 kb around peak centres.

Outputs:
- `qc_post_alignment/plots/heatmap_peaks.png`
- `qc_post_alignment/plots/profile_peaks.png`

#### FRiP over consensus peaks

In addition to per-replicate FRiP (Phase 1), each sample is also scored against
the shared consensus peak set. This gives a comparable metric across all samples
regardless of individual peak call quality.

Output: `qc_post_alignment/tables/frip_consensus.tsv`

---

## Output directory structure

```
qc_post_alignment/
|
+-- tables/
|   +-- qc_summary.tsv               # Per-sample: reads, dup%, mito%, peaks, FRiP
|   +-- qc_warnings.tsv              # Flagged samples with warning codes
|   +-- frip_consensus.tsv           # FRiP over merged consensus peak set
|   +-- consensus_sizeFactors.tsv    # DESeq2 consensus peak size factors
|   +-- track_normalization_metadata.tsv # CPM counts and DESeq2 scales
|   +-- per_chromosome_reads.tsv     # Read counts per chromosome per sample
|   +-- fingerprint_metrics.tsv      # deepTools fingerprint raw values
|
+-- plots/
|   +-- chromosome_coverage/
|   |   +-- <KEY>_karyogram.png          # Per-sample chromosome coverage
|   |   +-- karyogram_all_samples.png    # All samples in one panel
|   |   +-- <KEY>_100kb.bedGraph         # 100 kb bin coverage (intermediate)
|   |   +-- <KEY>_per_chrom.tsv          # Per-chromosome totals (intermediate)
|   +-- fingerprint.png                  # deepTools Lorenz curve
|   +-- correlation_heatmap_bins.png     # Genome-wide Spearman correlation
|   +-- correlation_heatmap_peaks.png    # Peak-centric Spearman correlation
|   +-- pca_bins.png                     # PCA on genome-wide bins
|   +-- pca_peaks.png                    # PCA on consensus peaks
|   +-- heatmap_peaks.png                # Signal heatmap over peaks
|   +-- profile_peaks.png                # Average profile over peaks
|
+-- matrices/
|   +-- bins_summary.npz                 # multiBamSummary bins matrix
|   +-- peaks_summary.npz                # multiBamSummary peaks matrix
|   +-- consensus_peak_counts.tsv        # Raw counts over consensus peak regions
|   +-- consensus_peak_normCounts.tsv    # DESeq2 normalized consensus counts
|   +-- signal_matrix.gz                 # computeMatrix output
|   +-- pca_bins.tsv                     # PCA coordinates (text)
|
+-- peak_sets/
|   +-- merged_narrow.bed                # Union of all narrow peaks
|   +-- merged_broad.bed                 # Union of all broad peaks
|   +-- consensus_peaks.bed              # Final merged consensus set
|
+-- logs/
    +-- post_alignment_qc_<date>.log
```

---

## Quality warning flags

Flags are written to `qc_warnings.tsv`. Thresholds are **assay-dependent** —
treat them as starting points, not hard cutoffs.

| Flag | Condition | Typical interpretation |
|---|---|---|
| `LOW_READS` | < 5 M aligned reads | Library too small; consider re-sequencing |
| `HIGH_DUPLICATION` | Dup rate > 80% | Over-amplification or very low input |
| `HIGH_MITO` | chrM reads > 10% | Poor cell lysis or mitochondrial contamination |
| `NO_PEAKS` | 0 peaks called | IP failure or extreme low-input; check fingerprint |
| `FEW_PEAKS` | < 200 peaks | Very low enrichment |
| `LOW_FRIP_NARROW` | FRiP narrow < 1% | Poor signal-to-noise |
| `LOW_FRIP_BROAD` | FRiP broad < 1% | As above for broad peaks |

> CUT&RUN and CUT&TAG typically show FRiP 0.1–0.6. ChIP-seq FRiP is usually
> 0.01–0.20. ATAC-seq FRiP depends on nucleosomal enrichment settings. Always
> interpret flags in the context of the specific assay and antibody target.

---

## Config parameter

Add to `config.conf`:

```bash
THREADS_DEEPTOOLS=16   # threads for deepTools steps; set to available CPUs / 2
```

---

## Comparison with legacy ChIPQC

| Feature | ChIPQC (removed) | deepTools QC (current) |
|---|---|---|
| FRiP | yes | yes (per-replicate + consensus) |
| Duplication rate | yes | yes (from Picard metrics) |
| Chromosome coverage karyogram | yes | yes (`plot_chrom_coverage.py`) |
| Fingerprint (Lorenz curve) | no | yes |
| PCA | no | yes (bins + peaks) |
| Sample correlation heatmap | no | yes (bins + peaks) |
| Signal heatmap over peaks | no | yes |
| Average profile over peaks | no | yes |
| Consensus/merged peak set | no | yes |
| Mitochondrial reads % | no | yes |
| Zero-peak sample handling | crashes | graceful — flagged, not dropped |
| Requires R / Bioconductor | yes | no |
| Requires RDS annotation files | yes | no |
| MultiQC integration | partial | yes (validated static deepTools plots embedded as custom content) |

---

## Rerunning Step 10 only

```bash
rm /path/to/project/.checkpoints/step10.done

nohup bash /path/to/ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/project/config/config.conf \
    >> /path/to/project/rerun_step10.log 2>&1 &
```

---

[← Blacklist filtering](11_blacklist_filtering.md) | [← Pipeline steps](06_pipeline_steps.md)
