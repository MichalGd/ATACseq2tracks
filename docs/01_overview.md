# 01 — Overview

[← Index](README.md) | [Next: Quick start →](02_quickstart.md)

---

## What is fastq2tracks?

**fastq2tracks** is a modular, checkpoint-resumable pipeline that converts raw FASTQ files from chromatin-profiling assays into analysis-ready outputs:

- Normalised **bigwig tracks** for genome browsers (individual replicates and merged)
- **MACS2 peaks** in both narrow and broad format for every IP sample
- **ChIPQC** quality reports with FRiP scores and coverage metrics
- **DiffBind-ready samplesheets** for differential binding analysis

A single **CSV samplesheet** describes every sample. A single **config file** sets all paths and parameters. The pipeline runs as a single bash command and can be safely interrupted and resumed — each step writes a checkpoint file, and completed steps are skipped on rerun.

---

## Supported assay types

| Assay | Notes |
|---|---|
| ChIP-seq | Primary target. SE and PE. hg38 and mm39. |
| ATAC-seq | Supported. Use `macs2_mode=both`; adjust fragment size in config. |
| CUT&RUN | Supported. Typically narrow peaks only. |
| CUT&TAG | Supported. Similar to CUT&RUN. |

---

## Workflow overview

```mermaid
flowchart TD
    A([Raw FASTQ files
SE or PE]) --> S0[Step 0 — Pre-flight
smoke_test.sh]
    S0 --> S1[Step 1 — FastQC raw
fastqc_batch.sh]
    S1 --> S2[Step 2 — Adapter trimming
trimgalore_batch.sh]
    S2 --> S3[Step 3 — FastQC trimmed
fastqc_batch.sh]
    S3 --> S4[Step 4 — Alignment
bowtie2_batch.sh]
    S4 --> S5[Step 5 — Deduplication
picard_dedup_batch.sh]
    S5 --> S6[Step 6 — Blacklist filtering
blacklist_filter_batch.sh]
    S6 --> S7[Step 7 — Coverage tracks
genomecoverage_batch.sh]
    S6 --> S8[Step 8 — Merged tracks
merge_replicates.sh]
    S6 --> S9[Step 9 — Peak calling
macs2_batch.sh]
    S7 --> OUT1[(bigwig/
bedGraph/)]
    S8 --> OUT2[(bigwig_merged/)]
    S9 --> OUT3[(peaks/
narrow + broad)]
    OUT3 --> S10[Step 10 — ChIPQC
run_chipqc.R]
    S6 --> S10
    OUT3 --> S11[Step 11 — DiffBind prep
prepare_diffbind.R]
    S6 --> S11
    OUT1 --> S12[Step 12 — UCSC tracks
create_ucsc_tracks.sh]
    OUT2 --> S12
    S10 --> S13[Step 13 — Report
generate_pipeline_report.sh]
    S11 --> S13
    S12 --> S13
    S13 --> FINAL([HTML report
DiffBind CSVs
UCSC trackdb.txt])

    style A fill:#d4edda,stroke:#28a745
    style FINAL fill:#cce5ff,stroke:#004085
    style S0 fill:#fff3cd,stroke:#856404
```

---

## Design principles

| Principle | Implementation |
|---|---|
| **Single entry point** | `fastq2tracks.sh --config config.conf` runs everything |
| **No hardcoded paths** | All paths and parameters in `config/config.conf` |
| **Safe resume** | `.checkpoints/stepN.done` flags; failed/removed steps rerun |
| **Mixed genomes** | hg38 and mm39 rows can coexist in one samplesheet |
| **SE and PE** | Read layout per sample via samplesheet column `layout` |
| **Replicate tracking** | Biological and technical replicates explicit in samplesheet |
| **Control linkage** | `control_id` column links each IP to its input |

---

## fastq2tracks vs rnaseq2tracksP

Both pipelines share the same samplesheet-driven, config-file-parameterised, checkpoint-resumable design and can be run on samples from the same experiment.

| Feature | **fastq2tracks** | **rnaseq2tracksP** |
|---|---|---|
| Assay | ChIP-seq, ATAC-seq, CUT&RUN, CUT&TAG | RNA-seq |
| Aligner | Bowtie2 (gapped-free) | STAR / HISAT2 (splice-aware) |
| Duplicate removal | Picard MarkDuplicates | Optional |
| Peak calling | MACS2 narrow + broad | — |
| QC module | ChIPQC (FRiP, TSS enrichment) | RSeQC / MultiQC |
| Track normalisation | Library-size / spike-in | TPM / CPM |
| Downstream prep | DiffBind samplesheets | DESeq2 / edgeR count matrices |
| Genome support | hg38, mm39 | hg38, mm39 |
| Config system | `config.conf` (bash source) | `config.conf` (bash source) |

---

[← Index](README.md) | [Next: Quick start →](02_quickstart.md)
