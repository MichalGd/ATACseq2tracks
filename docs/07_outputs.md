# 07 — Outputs

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)

---

## Output directory tree

> Track naming: CPM tracks use `_CPM.{bw,bedGraph}`. Consensus-scaled files use `_DESeq2Consensus.{bw,bedGraph}`. Robust filtering-policy tracks are stored under `bigwig_deseq2_robust_cpm/<permissive|intermediate|stringent>/` with the policy in the suffix. The legacy `_DESeq2RobustCPM.{bw,bedGraph}` stringent alias is retained for compatibility. Optional dm6-calibrated tracks use `_SpikeInDM6_Stringent.{bw,bedGraph}`.

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
├── bigwig/                            # Per-replicate fragment/read CPM tracks
│   ├── <sample_id>_bioR<N>_dedup_blFilt_CPM.bw
│   ├── <sample_id>_bioR<N>_dedup_blFilt_CPM.bedGraph
│   └── ucsc_tracks.txt
│
├── bigwig_deseq2_consensus/           # Reciprocal-size-factor coverage
│   ├── <sample_id>_bioR<N>_dedup_blFilt_DESeq2Consensus.bw
│   ├── <sample_id>_bioR<N>_dedup_blFilt_DESeq2Consensus.bedGraph
│   └── ucsc_tracks.txt
│
├── bigwig_deseq2_robust_cpm/          # DESeq2 robust CPM/FPM-style coverage
│   ├── <sample_id>_bioR<N>_dedup_blFilt_DESeq2RobustCPM.bw
│   ├── <sample_id>_bioR<N>_dedup_blFilt_DESeq2RobustCPM.bedGraph
│   ├── permissive/*_DESeq2RobustCPM_Permissive.{bw,bedGraph}
│   ├── intermediate/*_DESeq2RobustCPM_Intermediate.{bw,bedGraph}
│   ├── stringent/*_DESeq2RobustCPM_Stringent.{bw,bedGraph}
│   └── ucsc_tracks.txt
│
├── coverage_filtering_sensitivity/
│   ├── track_normalization_metadata.tsv
│   ├── permissive/{matrices,tables}/
│   └── intermediate/{matrices,tables}/
│
├── bigwig_spikein/stringent/           # Optional dm6-calibrated stringent host tracks
│   ├── <sample_id>_bioR<N>_SpikeInDM6_Stringent.bw
│   └── <sample_id>_bioR<N>_SpikeInDM6_Stringent.bedGraph
│
├── spikein/
│   ├── tables/spikein_normalization.tsv
│   ├── tables/spikein_warnings.tsv
│   ├── tables/spikein_provenance.tsv
│   └── logs/
│
├── bigwig_merged/                     # Condition-group merged CPM tracks
│   ├── <factor>_<condition>_<treatment>_<cell_type>_merged_CPM.{bw,bedGraph}
│   └── ucsc_tracks.txt
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
│   │   ├── track_normalization_metadata.tsv # Counts, constants and scales used for tracks
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
├── diffbind_results/
│   ├── diffbind_samplesheet_<genome>_narrow/  # Shared universe/model + comparisons/
│   └── diffbind_samplesheet_<genome>_broad/   # Shared universe/model + comparisons/
│
├── deseq2atac/
│   ├── deseq2atac_peak_type_summary.tsv
│   ├── broad/                         # Shared broad universe/model + comparisons/
│   └── narrow/                        # Shared narrow universe/model + comparisons/
│
├── logs/                              # Batch log files per step
│
└── reports/
    ├── differential_accessibility_summary.tsv
    ├── differential_accessibility_summary.html
    └── pipeline_report_<YYYYMMDD>.html
```

---

## Key output files explained

### bigwig tracks

BigWig files (`.bw`) contain binned coverage under the normalization named in each suffix. Load them directly into:
- [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)
- [IGV](https://igv.org/)
- [deepTools](https://deeptools.readthedocs.io/) for heatmaps and profile plots

For paired-end samples, coverage represents fragments counted once and extended over their inserts; single-end coverage represents reads. All five families use canonical autosomes plus X/Y. All robust policies use the same fixed consensus BED but obtain separate count matrices, DESeq2 size factors and cohort constants. Within one policy, robust CPM is a cohort-wide rescaling of its DESeq2-consensus signal. The permissive and intermediate tracks are sensitivity outputs; none of the DESeq2 families is an absolute cross-study calibration.

The optional dm6 family is external-reference calibration rather than DESeq2
normalization. It scales raw stringent host coverage by
`SPIKEIN_SCALE_TARGET × spikein_to_host_ratio / retained_dm6_observations`.
Interpret it only with `spikein_normalization.tsv` and `spikein_warnings.tsv`.
It can support global-shift comparisons when the reference was added correctly,
but it does not remove study or protocol batch effects.

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

### DESeq2ATAC outputs

`deseq2atac/broad/` and `deseq2atac/narrow/` contain separate BED4 consensus
sets, support tables, raw and normalized count matrices, sample/library metadata,
DESeq2 size factors, serialized objects, session information and shared PNG/PDF
diagnostics. Complete and FDR-filtered result tables and contrast-specific plots
are under `comparisons/<comparison_id>/`. The top-level
`deseq2atac_peak_type_summary.tsv` records broad/narrow status and pair counts.

For each peak type, per-sample peaks are canonical/blacklist filtered and reduced,
then disjoined across samples. Atomic segments supported by at least
`DESEQ2ATAC_MIN_SAMPLES` biological samples are retained and adjacent retained
segments are reduced into the final non-overlapping consensus. This is not a
one-sample union. Narrow DESeq2ATAC regions retain narrowPeak-derived boundaries;
unlike DiffBind, they are not recentered to 201-bp summit windows.

The complete result table preserves independent-filtered `padj=NA` entries. The
significant table remains a valid header-only compressed table when no region
passes `DESEQ2ATAC_ALPHA`; this is a successful result, not a workflow failure.

All non-control biological samples contribute to each consensus and all-sample
raw count matrix. Only conditions with at least two biological samples enter the
model. Every pair among eligible conditions is exported; a singleton condition
is not compared but still receives all other workflow outputs.

When built-in annotation is enabled, each peak-type directory also contains a
compressed consensus annotation table; its gene-context, nearest-TSS and
cCRE fields are joined to pair-level results by default. With
`RUN_CCRE_ANNOTATION=false`, the tables retain GTF fields without cCRE class
assignments.

See [Peak annotation](16_peak_annotation.md) for every column, gene-context
precedence, primary-cCRE selection, source provenance and interpretation limits.

`reports/differential_accessibility_summary.{tsv,html}` records status, tested,
significant, higher-in-numerator and higher-in-reference counts for every
method, peak type and condition pair. Counts are not deduplicated loci and should
not be added across methods or peak types.

The unified MultiQC report excludes MultiQC's native deepTools parser because
MultiQC 1.35 can reject `plotPCA` auxiliary tables and overwrite repeated sample
names. The authoritative deepTools tables remain under `qc_post_alignment/`;
their already-rendered PCA, correlation, fingerprint and peak-signal PNGs are
embedded in MultiQC as self-contained custom-content sections. The adjacent
`*.multiqc.log` is checked for caught module, validation and colour failures.

---

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)
