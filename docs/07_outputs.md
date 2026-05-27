# 07 — Outputs

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)

---

## Output directory tree

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
│   └── <sample_id>_bioR<N>_dedup_blFilt.bw
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
├── chipqc/
│   ├── ChIPQC_hg38_narrow/
│   │   ├── ChIPQCReport.html          # Main QC report
│   │   ├── chipqc_samplesheet.csv     # Samples that passed BAM check
│   │   ├── chipqc_metrics_summary.csv # QC metrics table
│   │   └── chipqc_frip.csv            # FRiP scores per sample
│   ├── ChIPQC_hg38_broad/
│   └── ChIPQC_mm39_narrow/            # Only if mm39 samples present
│   └── ChIPQC_mm39_broad/
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

### ChIPQC FRiP scores (`chipqc_frip.csv`)

FRiP = Fraction of Reads in Peaks.

| FRiP value | Interpretation |
|---|---|
| < 0.01 | Poor enrichment — check antibody and library quality |
| 0.01–0.05 | Acceptable for some marks (H3K27me3, broad marks) |
| > 0.05 | Good |
| > 0.20 | Excellent (typical for H3K27ac, H3K4me3, CTCF) |

### DiffBind samplesheets

Ready-to-use input for the Bioconductor `DiffBind` package.
See [Downstream: DiffBind](08_diffbind.md).

---

[← Pipeline steps](06_pipeline_steps.md) | [Next: DiffBind →](08_diffbind.md)
