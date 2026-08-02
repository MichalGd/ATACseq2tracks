# 07 — Outputs

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)

---

## Output directory tree

> v3.2.0 naming: RPM tracks end in `_RPM.bw`; consensus-scaled tracks are stored in `bigwig_deseq2_consensus/` and end in `_DESeq2Consensus.bw`. ATAC-specific QC is under `qc_post_alignment/atac_qc/`, including `tables/ataqv_selected_metrics.tsv`, `tables/nucleosome_periodicity_metrics.tsv`, compressed `ataqv_metrics/*.ataqv.json.gz`, static fragment-periodicity PNG/PDF plots, and the optional `ataqv_viewer/`.

```
<OUTPUT_DIR>/
│
├── .checkpoints/                      # Step completion flags
│   ├── step1.done
│   ├── step2.done
│   └── ...
│
├── fastQC/
│   ├── fastQC_unTrimmed/              # Per-sample FastQC HTML + zip (raw)
│   └── fastQC_trimmed/               # Per-sample FastQC HTML + zip (trimmed)
│
├── multiQC/
│   ├── multiQC_unTrimmed/             # MultiQC HTML report for raw reads
│   ├── multiQC_trimmed/               # MultiQC HTML report for trimmed reads
│   ├── multiQC_alignments/            # Bowtie2 alignment summary
│   └── multiQC_deduplication/         # Picard duplication metrics
│
├── trimmedFastq/                      # Trimmed FASTQs (deleted if KEEP_TRIMMED_FASTQ=false)
│
├── bams/                              # Sorted, indexed BAMs from Bowtie2
│                                      # Deleted if KEEP_INTERMEDIATE_BAMS=false
│
├── dedupBams/                         # BAMs after Picard deduplication
│   └── <sample_id>_bioR<N>_dedup.bam
│
├── filteredBams/                      # BAMs after blacklist filtering
│   └── <sample_id>_bioR<N>_dedup_blFilt.bam
│
├── bedGraph/                          # Raw coverage bedGraphs
├── NormBedGraph/                      # Library-size normalised bedGraphs
│
├── bigwig/                            # Per-replicate bigwig tracks
│   └── <sample_id>_bioR<N>_dedup_blFilt_RPM.bw
│
├── bigwig_deseq2_consensus/           # DESeq2 consensus-peak-scaled bigWig tracks
│   └── <sample_id>_bioR<N>_dedup_blFilt_DESeq2Consensus.bw
│
├── bigwig_merged/                     # Condition-group merged bigwig tracks
│   └── <factor>__<condition>__<treatment>__<cell_type>__<genome>.bw
│
├── peaks/
│   ├── per_replicate/
│   │   └── <sample_id>/
│   │       ├── narrow/
│   │       │   ├── <sample_id>_peaks.narrowPeak
│   │       │   ├── <sample_id>_summits.bed
│   │       │   └── <sample_id>_peaks.xls
│   │       └── broad/
│   │           ├── <sample_id>_peaks.broadPeak
│   │           └── <sample_id>_peaks.xls
│   └── pooled/
│       └── <factor>__<condition>__<treatment>__<cell_type>__<genome>/
│           ├── narrow/
│           └── broad/
│
├── qc_post_alignment/                 # deepTools post-alignment QC (Step 10)
│   ├── tables/
│   │   ├── qc_summary.tsv             # Per-sample: reads, dup%, mito%, peaks, FRiP
│   │   ├── qc_warnings.tsv            # Flagged samples with warning codes
│   │   ├── frip_consensus.tsv         # FRiP over merged consensus peak set
│   │   ├── consensus_sizeFactors.tsv  # DESeq2 consensus peak size factors
│   │   ├── per_chromosome_reads.tsv   # Read counts per chromosome per sample
│   │   └── fingerprint_metrics.tsv    # deepTools fingerprint raw values
│   ├── plots/
│   │   ├── chromosome_coverage/
│   │   │   ├── <KEY>_karyogram.png    # Per-sample chromosome coverage (ChIPQC style)
│   │   │   └── karyogram_all_samples.png
│   │   ├── fingerprint.png            # Lorenz curve (deepTools plotFingerprint)
│   │   ├── correlation_heatmap_bins.png
│   │   ├── correlation_heatmap_peaks.png
│   │   ├── pca_bins.png
│   │   ├── pca_peaks.png
│   │   ├── heatmap_peaks.png
│   │   └── profile_peaks.png
│   ├── matrices/
│   │   ├── bins_summary.npz
│   │   ├── peaks_summary.npz
│   │   ├── consensus_peak_counts.tsv
│   │   ├── consensus_peak_normCounts.tsv
│   │   └── signal_matrix.gz
│   ├── peak_sets/
│   │   ├── merged_narrow.bed
│   │   ├── merged_broad.bed
│   │   └── consensus_peaks.bed
│   └── logs/
│       └── post_alignment_qc_<date>.log
│
├── diffbind/
│   ├── diffbind_samplesheet_hg38_narrow.csv
│   ├── diffbind_samplesheet_hg38_broad.csv
│   └── (mm39 equivalents if applicable)
│
├── logs/                              # Batch log files per step
│
└── reports/
    ├── ucsc_trackdb.txt
    └── pipeline_report_<YYYYMMDD>.html
```

---

## Key output files explained

### bigwig tracks

Bigwig files (`.bw`) contain per-base coverage normalised by library size. Load them directly into:
- [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)
- [IGV](https://igv.org/)
- [deepTools](https://deeptools.readthedocs.io/) for heatmaps and profile plots

Per-replicate tracks are in `bigwig/`. Merged (averaged across replicates) tracks are in `bigwig_merged/`.

### MACS2 peak files

| File | Description |
|---|---|
| `*_peaks.narrowPeak` | BED6+4 format: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak |
| `*_peaks.broadPeak` | BED6+3 format: same without summit position |
| `*_summits.bed` | Single-nucleotide summit positions (narrow only) — useful for motif analysis |

> Use `-log10(qValue)` (column 9) as the primary filtering criterion. A common threshold is `qValue > 2` (i.e. FDR < 0.01).

### Post-alignment QC (`qc_post_alignment/`)

The main summary is `tables/qc_summary.tsv`. Key columns:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier from samplesheet |
| `aligned_reads` | Total reads in filtered BAM |
| `mito_pct` | Mitochondrial read fraction (%) |
| `dup_rate` | Duplication rate from Picard metrics |
| `n_narrow_peaks` | Number of MACS2 narrow peaks |
| `n_broad_peaks` | Number of MACS2 broad peaks |
| `frip_narrow` | FRiP over per-replicate narrow peaks |
| `frip_broad` | FRiP over per-replicate broad peaks |
| `frip_consensus` | FRiP over merged consensus peak set |
| `warnings` | Comma-separated warning flags (e.g. `NO_PEAKS,LOW_READS`) |

FRiP thresholds (assay-dependent — see [Post-alignment QC](12_post_alignment_qc.md)):

| FRiP value | Interpretation |
|---|---|
| < 0.01 | Poor enrichment — check antibody and library quality |
| 0.01–0.05 | Acceptable for some marks (H3K27me3, broad marks) |
| > 0.05 | Good |
| > 0.20 | Excellent (typical for H3K27ac, H3K4me3, CTCF) |

> CUT&RUN and CUT&TAG typically show FRiP 0.1–0.6. ATAC-seq FRiP depends on nucleosomal enrichment settings.

### DiffBind samplesheets

Ready-to-use input for the Bioconductor `DiffBind` package.
See [Downstream: DiffBind](08_diffbind.md).

---

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)
