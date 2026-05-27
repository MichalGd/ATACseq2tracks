# fastq2tracks v3.0.4

**ChIP-seq / NGS track-generation, peak-calling, and QC workflow**

[![Version](https://img.shields.io/badge/version-3.0.4-blue)](CHANGELOG.md)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Genome](https://img.shields.io/badge/genome-hg38%20%7C%20mm39-orange)]()

fastq2tracks processes raw FASTQ files from ChIP-seq and related NGS assays (ATAC-seq, CUT&RUN, CUT&TAG) through a full analysis pipeline — adapter trimming, alignment, deduplication, blacklist filtering, bigwig track generation, MACS2 peak calling (narrow **and** broad for every IP sample), ChIPQC quality assessment, and DiffBind-ready output preparation — all driven from a single CSV samplesheet.

---

## Table of Contents

- [What's new in v3.0](#whats-new-in-v30)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Samplesheet format](#samplesheet-format)
- [Configuration file](#configuration-file)
- [Pipeline steps](#pipeline-steps)
- [Output structure](#output-structure)
- [Rerunning individual steps](#rerunning-individual-steps)
- [Server resource recommendations](#server-resource-recommendations)
- [Downstream: DiffBind](#downstream-diffbind)
- [Known issues and workarounds](#known-issues-and-workarounds)
- [Citation](#citation)

---

## What's new in v3.0

| Feature | v2.1 | v3.0 |
|---|---|---|
| Read layout | PE only | SE **and** PE per sample |
| Blacklist filtering | No | Yes — hg38 + mm39 |
| Replicates | Not tracked | Biological + technical replicates |
| Control linkage | Manual | Explicit `control_id` in samplesheet |
| MACS2 | No | Narrow **and** broad for every IP sample |
| ChIPQC | No | BiocParallel, pre-built hg38/mm39 annotations |
| Merged replicate tracks | No | `bigwig_merged/` for each condition group |
| DiffBind prep | No | Samplesheets for narrow + broad per genome |
| Smoke test | No | Pre-flight check before every run |
| Conda environment | No | `environment.yml` provided |
| Samplesheet validation | No | `validate_samplesheet.py` |
| Checkpoint system | No | Per-step `.done` files, safe resume |
| Script paths | Hardcoded absolute | Relative to `fastq2tracks/scripts/` |

---

## Repository layout

```
fastq2tracks/
├── fastq2tracks.sh            # Master script — run this
├── environment.yml            # Conda environment for collaborators
├── CHANGELOG.md
├── LICENSE
├── README.md                  # This file
├── config/
│   ├── config.conf            # All paths, tool settings, thread counts
│   ├── samplesheet_template.csv
│   └── samplesheet_example.csv
└── scripts/
    ├── smoke_test.sh          # Pre-flight checks
    ├── validate_samplesheet.py
    ├── trimgalore_batch.sh    # Adapter trimming (SE+PE, tech-rep merging)
    ├── fastqc_batch.sh        # FastQC on raw and trimmed FASTQs
    ├── bowtie2_align.sh       # Single-sample aligner wrapper
    ├── bowtie2_batch.sh       # Samplesheet-driven alignment (hg38 + mm39)
    ├── picard_dedup.sh        # Single-sample Picard MarkDuplicates wrapper
    ├── picard_dedup_batch.sh  # Batch deduplication
    ├── blacklist_filter.sh    # Single-sample blacklist filtering
    ├── blacklist_filter_batch.sh
    ├── genomecoverage_single.sh  # Single-sample bedGraph + bigwig
    ├── genomecoverage_batch.sh   # Batch track generation
    ├── merge_replicates.sh    # Merge replicates → normalized tracks
    ├── macs2_peaks.sh         # Single-sample MACS2 (narrow + broad)
    ├── macs2_batch.sh         # Batch peak calling
    ├── run_chipqc.R           # ChIPQC with BiocParallel
    ├── prepare_diffbind.R     # DiffBind samplesheet generation
    ├── create_ucsc_tracks.sh  # UCSC track hub helper
    └── generate_pipeline_report.sh
```

> **Always run the workflow from the parent directory** containing `fastq2tracks/`, not from inside the `fastq2tracks/` folder itself.

---

## Requirements

### Software (managed via conda — see [Installation](#installation))

| Tool | Version tested | Role |
|---|---|---|
| bowtie2 | ≥ 2.5 | Alignment |
| samtools | ≥ 1.18 | BAM processing |
| bedtools | ≥ 2.31 | Interval operations |
| trim_galore | 2.2.0 (Oxidized Edition) | Adapter trimming |
| fastqc | ≥ 0.12 | Read QC |
| macs2 | ≥ 2.2.9 | Peak calling |
| multiqc | ≥ 1.20 | Aggregated QC reports |
| picard | ≥ 3.1 (`.jar`) | Deduplication |
| bedGraphToBigWig | UCSC kent | bigwig conversion |
| R + Bioconductor | R ≥ 4.3 | ChIPQC, DiffBind prep |

### R packages

```r
BiocManager::install(c(
  "ChIPQC", "DiffBind", "BiocParallel",
  "GenomicAlignments", "rtracklayer"
))
install.packages(c("ggplot2", "dplyr"))
```

### Reference files (local paths — configure in `config/config.conf`)

| File | Config variable |
|---|---|
| hg38 bowtie2 index | `INDEX_HG38` |
| mm39 bowtie2 index | `INDEX_MM39` |
| hg38 chrom sizes | `CHROM_SIZES_HUMAN` |
| mm39 chrom sizes | `CHROM_SIZES_MOUSE` |
| hg38 blacklist BED (ENCODE ENCFF356LFX) | `BLACKLIST_HG38` |
| mm39 blacklist BED (Boyle lab) | `BLACKLIST_MM39` |
| hg38 ChIPQC annotation RDS | `CHIPQC_ANNOTATION_HG38` |
| mm39 ChIPQC annotation RDS | `CHIPQC_ANNOTATION_MM39` |
| hg38 ChIPQC blacklist RDS | `CHIPQC_BLACKLIST_HG38_RDS` |
| mm39 ChIPQC blacklist RDS | `CHIPQC_BLACKLIST_MM39_RDS` |

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/MichalGd/fastq2tracks.git
cd fastq2tracks

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate fastq2tracks

# 3. Make all scripts executable
chmod +x fastq2tracks.sh scripts/*.sh scripts/*.py scripts/*.R

# 4. Copy and edit the config for your project
mkdir -p /path/to/my_project/config
cp config/config.conf              /path/to/my_project/config/config.conf
cp config/samplesheet_template.csv /path/to/my_project/config/samplesheet.csv
# Edit both files — set reference paths, output directory, thread counts
```

---

## Quick start

```bash
# 1. Activate environment
conda activate fastq2tracks

# 2. Validate your samplesheet (optional but recommended)
python3 fastq2tracks/scripts/validate_samplesheet.py /path/to/config/samplesheet.csv

# 3. Run the full pipeline
cd /path/to/containing_folder   # parent of fastq2tracks/
nohup bash fastq2tracks/fastq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    > /path/to/my_project/run.log 2>&1 &

# Monitor progress
tail -f /path/to/my_project/run.log
```

---

## Samplesheet format

The samplesheet is a CSV with one row per library. Use `tech_replicate` to handle multiple sequencing runs of the same library — they are merged at the trimming step before alignment.

```
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,
cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,
blacklist,chipqc_annotation,output_prefix
```

| Column | Description | Values |
|---|---|---|
| `sample_id` | Unique sample identifier | string |
| `fastq_1` | Absolute path to R1 FASTQ | path |
| `fastq_2` | Absolute path to R2 FASTQ (PE only) | path or empty |
| `layout` | Read layout | `SE` or `PE` |
| `genome` | Reference genome | `hg38` or `mm39` |
| `assay` | Assay type | `ChIP-seq`, `ATAC-seq`, etc. |
| `factor` | Antibody target or histone mark | e.g. `H3K27ac`, `CTCF` |
| `condition` | Biological condition | e.g. `treated`, `rest` |
| `treatment` | Fixation protocol or treatment | e.g. `DSG_FA`, `FA` |
| `cell_type` | Cell line or tissue | e.g. `NHEK`, `LymphocyteT` |
| `replicate` | Biological replicate number | integer |
| `tech_replicate` | Technical replicate number | integer (1 if none) |
| `is_control` | Is this an input/IgG control? | `TRUE` or `FALSE` |
| `control_id` | `sample_id` of the matched control | string or empty for controls |
| `macs2_mode` | Peak calling mode | `both`, `narrow`, `broad`, `none` |
| `blacklist` | Path to blacklist BED file | path |
| `chipqc_annotation` | Path to ChIPQC annotation RDS | path |
| `output_prefix` | Prefix for output files | string |

### Example rows

```csv
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,chipqc_annotation,output_prefix
NHEK_Input_bioR1,/data/NHEK_Input_R1.fastq.gz,,SE,hg38,ChIP-seq,Input,untreated,none,NHEK,1,1,TRUE,,none,/ref/blacklist_hg38.bed,/ref/anno_hg38_chipqc.rds,NHEK_Input_bioR1
NHEK_H3K27ac_day0_bioR1,/data/NHEK_H3K27ac_d0_R1.fastq.gz,,SE,hg38,ChIP-seq,H3K27ac,day0,none,NHEK,1,1,FALSE,NHEK_Input_bioR1,both,/ref/blacklist_hg38.bed,/ref/anno_hg38_chipqc.rds,NHEK_H3K27ac_day0_bioR1
Tco_rest_1_SATB1_bioR1,/data/Tco_SATB1_r1_R1.fastq.gz,/data/Tco_SATB1_r1_R2.fastq.gz,PE,hg38,ChIP-seq,SATB1,rest,DSG_FA,LymphocyteT,1,1,FALSE,Tco_rest_1_Input_bioR1,both,/ref/blacklist_hg38.bed,/ref/anno_hg38_chipqc.rds,Tco_rest_1_SATB1_bioR1
```

> **Controls:** Set `is_control=TRUE` and `macs2_mode=none`. Leave `control_id` empty for control rows.

> **Technical replicates:** If a library was sequenced in multiple runs, add one row per run with the same `sample_id` and `replicate` but different `tech_replicate` values (1, 2, …). They are merged automatically before trimming.

> **MACS2 modes:** `both` always produces narrow **and** broad peaks. Use `none` for control/input samples only.

---

## Configuration file

Copy `config/config.conf` and edit for your server. Minimal required variables:

```bash
# ── Project ───────────────────────────────────────────────────────────────────
SAMPLESHEET="/path/to/config/samplesheet.csv"
OUTPUT_DIR="/path/to/analysis/"

# ── Reference genomes ─────────────────────────────────────────────────────────
INDEX_HG38="/path/to/bowtie2/hg38"
INDEX_MM39="/path/to/bowtie2/mm39"
CHROM_SIZES_HUMAN="/path/to/hs38n.chrom.sizes"
CHROM_SIZES_MOUSE="/path/to/mm39n.chrom.sizes"

# ── Blacklists ────────────────────────────────────────────────────────────────
BLACKLIST_HG38="/path/to/blacklist_hg38_ENCFF356LFX.bed"
BLACKLIST_MM39="/path/to/blacklist_mm39_Boyle.bed"

# ── ChIPQC pre-built annotation RDS objects ───────────────────────────────────
CHIPQC_ANNOTATION_HG38="/path/to/anno_hg38_chipqc.rds"
CHIPQC_ANNOTATION_MM39="/path/to/anno_mm39_chipqc.rds"
CHIPQC_BLACKLIST_HG38_RDS="/path/to/blacklist_hg38.rds"
CHIPQC_BLACKLIST_MM39_RDS="/path/to/blacklist_mm39.rds"

# ── Tools ─────────────────────────────────────────────────────────────────────
PICARD_JAR="/path/to/picard.jar"
R_BIN="Rscript"

# ── Threads ───────────────────────────────────────────────────────────────────
THREADS_PARALLEL_JOBS=8
THREADS_ALIGN=16
THREADS_SAMTOOLS=16
THREADS_TRIMGALORE=8
THREADS_FASTQC=10
THREADS_BIGWIG=16
THREADS_CHIPQC=20

# ── Cleanup ───────────────────────────────────────────────────────────────────
KEEP_INTERMEDIATE_BAMS=false
KEEP_TRIMMED_FASTQ=false
```

---

## Pipeline steps

| Step | Script | Input → Output |
|---|---|---|
| 0 Pre-flight | `smoke_test.sh` | config + samplesheet → pass/fail |
| 1 FastQC raw | `fastqc_batch.sh` | raw FASTQs → `fastQC/fastQC_unTrimmed/` |
| 2 Trimming | `trimgalore_batch.sh` | raw FASTQs → `trimmedFastq/` |
| 3 FastQC trimmed | `fastqc_batch.sh` | trimmed FASTQs → `fastQC/fastQC_trimmed/` |
| 4 Alignment | `bowtie2_batch.sh` | trimmed FASTQs → `bams/*.bam` |
| 5 Deduplication | `picard_dedup_batch.sh` | `bams/` → `dedupBams/*_dedup.bam` |
| 6 Blacklist filter | `blacklist_filter_batch.sh` | `dedupBams/` → `filteredBams/*_filtered.bam` |
| 7 Coverage tracks | `genomecoverage_batch.sh` | `filteredBams/` → `bigwig/*.bw` + `bedGraph/` |
| 8 Merged tracks | `merge_replicates.sh` | `filteredBams/` → `bigwig_merged/*.bw` |
| 9 Peak calling | `macs2_batch.sh` | `filteredBams/` → `peaks/*/narrow/` + `.../broad/` |
| 10 ChIPQC | `run_chipqc.R` | filtered BAMs + peaks → `chipqc/` |
| 11 DiffBind prep | `prepare_diffbind.R` | filtered BAMs + peaks → `diffbind/*.csv` |
| 12 UCSC tracks | `create_ucsc_tracks.sh` | `bigwig/` → `reports/ucsc_trackdb.txt` |
| 13 Report | `generate_pipeline_report.sh` | all outputs → `reports/pipeline_report_<date>.html` |

---

## Output structure

```
<OUTPUT_DIR>/
├── .checkpoints/               # Step completion flags
├── fastQC/fastQC_unTrimmed/
├── fastQC/fastQC_trimmed/
├── multiQC/                    # MultiQC HTML reports per step
├── trimmedFastq/               # *_trimmed.fq.gz (SE) | *_val_1/2.fq.gz (PE)
├── bams/                       # <KEY>.bam + .bai
├── dedupBams/                  # <KEY>_dedup.bam
├── filteredBams/               # <KEY>_filtered.bam (blacklist removed)
├── bedGraph/
├── NormBedGraph/
├── bigwig/                     # Per-replicate tracks
├── bigwig_merged/              # Condition-merged tracks
├── peaks/
│   ├── per_replicate/<sample_id>/narrow/   # *_peaks.narrowPeak
│   ├── per_replicate/<sample_id>/broad/    # *_peaks.broadPeak
│   └── pooled/<group>/narrow/ and /broad/
├── chipqc/
│   ├── ChIPQC_report_hg38_narrow/
│   ├── ChIPQC_report_hg38_broad/
│   └── (mm39 equivalents)
├── diffbind/
│   ├── diffbind_samplesheet_hg38_narrow.csv
│   ├── diffbind_samplesheet_hg38_broad.csv
│   └── (mm39 equivalents)
├── logs/
└── reports/
```

---

## Rerunning individual steps

```bash
# Rerun one step — e.g. step 9 (MACS2)
rm /path/to/analysis/.checkpoints/step9.done
nohup bash fastq2tracks/fastq2tracks.sh \
    --config /path/to/config/config.conf >> run_rerun.log 2>&1 &

# Rerun everything from scratch
rm -rf /path/to/analysis/.checkpoints/
```

---

## Server resource recommendations

| Step | Config variable | Recommended value |
|---|---|---|
| Parallel samples | `THREADS_PARALLEL_JOBS` | 8 |
| Bowtie2 | `THREADS_ALIGN` | 16 |
| Samtools sort | `THREADS_SAMTOOLS` | 16 |
| TrimGalore | `THREADS_TRIMGALORE` | 8 |
| FastQC | `THREADS_FASTQC` | 10 |
| BigWig | `THREADS_BIGWIG` | 16 |
| ChIPQC BiocParallel | `THREADS_CHIPQC` | 20 |

Peak concurrent load: 8 × 16 = 128 threads (leaves headroom on a 140-logical-core server).

---

## Downstream: DiffBind

```r
library(DiffBind)

# Sharp marks (H3K27ac, H3K4me3, CTCF)
dba <- dba(sampleSheet = "diffbind/diffbind_samplesheet_hg38_narrow.csv")

# Broad marks (H3K27me3, H3K9me3, H3K36me3)
dba <- dba(sampleSheet = "diffbind/diffbind_samplesheet_hg38_broad.csv")
```

All DiffBind required columns are pre-filled: `SampleID`, `Factor`, `Condition`, `Treatment`, `Tissue`, `Replicate`, `bamReads`, `bamControl`, `Peaks`, `PeakCaller`.

---

## Known issues and workarounds

| Symptom | Root cause | Fix applied in |
|---|---|---|
| `THREADS_BOWTIE2: unbound variable` | Config uses `THREADS_ALIGN` | v3.0.4 |
| `trim_galore: error: unknown option -o` | v2.2.0 uses `--output_dir` | v3.0.4 |
| TrimGalore runs but produces no renamed output | `mv` in subshell swallowed errors; stem computed from temp path | v3.0.4 |
| Picard expects `*.sorted_stChr.bam` | Legacy naming from v2.1 | v3.0.4 |
| `samtools sort: failed to read header` | bowtie2 exited immediately on empty `-p` | v3.0.4 (same as THREADS fix) |
| FastQC trimmed finds 0 files | Downstream of TrimGalore rename bug | v3.0.4 |

---

## Citation

If you use fastq2tracks in your work, please cite:

```
Golebiewski M. fastq2tracks — ChIP-seq track-generation and QC workflow, v3.0.
https://github.com/MichalGd/fastq2tracks (2026)
```

---

*fastq2tracks is developed and maintained by Michal Golebiewski.  
Bug reports and contributions are welcome via [GitHub Issues](https://github.com/MichalGd/fastq2tracks/issues).*
