# 04 — Input files

[← Installation](03_installation.md) | [Next: Running →](05_running.md)

---

## Samplesheet

The samplesheet is a **comma-separated CSV** with one row per sequencing library. It is the
primary driver of the entire pipeline — every step reads it to know which samples to process,
which genome to use, how to pair IPs with controls, and what type of peaks to call.

> **v3.1.0 change:** The `chipqc_annotation` column (column 17 in v3.0.x) has been **removed**.
> The samplesheet now has **17 columns**. If you are upgrading from v3.0.x, remove this column
> from your existing samplesheets before running v3.1.0.

### Column reference

| # | Column | Type | Required | Description |
|---|---|---|---|---|
| 1 | `sample_id` | string | ✓ | Unique identifier for this library. Used as the key for all output file names. |
| 2 | `fastq_1` | path | ✓ | Absolute path to R1 FASTQ (`.fastq.gz`). |
| 3 | `fastq_2` | path | PE only | Absolute path to R2 FASTQ. Leave empty for SE samples. |
| 4 | `layout` | `SE`\|`PE` | ✓ | Single-end or paired-end. |
| 5 | `genome` | `hg38`\|`mm39` | ✓ | Reference genome for alignment and peak calling. |
| 6 | `assay` | string | ✓ | Assay type, e.g. `ChIP-seq`, `ATAC-seq`, `CUT&RUN`, `ChIPmentation`. |
| 7 | `factor` | string | ✓ | Antibody target or histone mark, e.g. `H3K27ac`, `CTCF`, `Input`. |
| 8 | `condition` | string | ✓ | Biological condition, e.g. `day0`, `treated`, `rest`. |
| 9 | `treatment` | string | ✓ | Fixation protocol or treatment, e.g. `FA`, `DSG_FA`, `none`. |
| 10 | `cell_type` | string | ✓ | Cell line or tissue, e.g. `NHEK`, `LymphocyteT`. |
| 11 | `replicate` | integer | ✓ | Biological replicate number: `1`, `2`, `3`, … |
| 12 | `tech_replicate` | integer | ✓ | Technical replicate (sequencing run) number. Use `1` if the library was sequenced once. |
| 13 | `is_control` | `TRUE`\|`FALSE` | ✓ | `TRUE` for input DNA or IgG controls; `FALSE` for IP samples. |
| 14 | `control_id` | string | IP only | The `sample_id` of the matched control. Must match exactly. Leave empty for control rows. |
| 15 | `macs2_mode` | `both`\|`narrow`\|`broad`\|`none` | ✓ | Peak calling mode. Use `none` for all control samples. |
| 16 | `blacklist` | path | ✓ | Absolute path to the blacklist BED file for this genome. |
| 17 | `output_prefix` | string | ✓ | Prefix for output files. Typically same as `sample_id`. |

### Rules and conventions

**Controls**
- Set `is_control=TRUE` and `macs2_mode=none`.
- Leave `control_id` empty.
- Controls still go through all preprocessing steps (trimming → alignment → dedup → filtering → tracks).

**IP samples**
- Set `is_control=FALSE`.
- `control_id` must exactly match the `sample_id` of the appropriate control.
- `macs2_mode=both` runs MACS3 in both narrow and broad mode (recommended for most marks).

**Technical replicates**
- Add one row per sequencing run with the same `sample_id`, `replicate`, and different
  `tech_replicate` (1, 2, …).
- The pipeline merges their FASTQs before trimming.

**Mixed genomes**
- hg38 and mm39 rows can coexist in the same samplesheet.
- Steps 8–11 iterate over each genome separately.

**macs2_mode guidance**

| Mark / Assay | Recommended mode |
|---|---|
| H3K27ac, H3K4me3, H3K4me1, p63, CTCF, SATB1 | `both` |
| H3K27me3, H3K9me3, H3K36me3 | `both` (broad peaks are primary) |
| ATAC-seq | `narrow` |
| CUT&RUN, CUT&TAG, ChIPmentation | `narrow` (or `both`) |
| Input / IgG | `none` |

---

### Example samplesheet

```csv
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
NHEK_Input_bioR1,/data/NHEK_Input_R1.fastq.gz,,SE,hg38,ChIP-seq,Input,day0,none,NHEK,1,1,TRUE,,none,/ref/blacklist_hg38.bed,NHEK_Input_bioR1
NHEK_H3K27ac_day0_bioR1,/data/NHEK_H3K27ac_d0_R1.fastq.gz,,SE,hg38,ChIP-seq,H3K27ac,day0,none,NHEK,1,1,FALSE,NHEK_Input_bioR1,both,/ref/blacklist_hg38.bed,NHEK_H3K27ac_day0_bioR1
Tco_rest_1_SATB1_bioR1,/data/Tco_SATB1_R1.fastq.gz,/data/Tco_SATB1_R2.fastq.gz,PE,hg38,ChIP-seq,SATB1,rest,DSG_FA,LymphocyteT,1,1,FALSE,Tco_rest_1_Input_bioR1,both,/ref/blacklist_hg38.bed,Tco_rest_1_SATB1_bioR1
Tco_rest_1_Input_bioR1,/data/Tco_Input_R1.fastq.gz,/data/Tco_Input_R2.fastq.gz,PE,hg38,ChIP-seq,Input,rest,DSG_FA,LymphocyteT,1,1,TRUE,,none,/ref/blacklist_hg38.bed,Tco_rest_1_Input_bioR1
```

> Use `python3 scripts/validate_samplesheet.py samplesheet.csv` to check for errors before running.

---

## Configuration file

`config/config.conf` is a bash script that is `source`d by the master script and all subscripts.
Copy the template and edit for your server.

### All configuration parameters

#### Project paths

| Variable | Description |
|---|---|
| `SAMPLESHEET` | Absolute path to the samplesheet CSV |
| `OUTPUT_DIR` | Absolute path to the root output directory (created if absent) |

#### Reference genomes and indices

| Variable | Description |
|---|---|
| `INDEX_HG38` | Bowtie2 index prefix for hg38 |
| `INDEX_MM39` | Bowtie2 index prefix for mm39 |
| `CHROM_SIZES_HUMAN` | hg38 chromosome sizes file |
| `CHROM_SIZES_MOUSE` | mm39 chromosome sizes file |

#### Blacklist regions

| Variable | Description |
|---|---|
| `BLACKLIST_HG38` | hg38 blacklist BED (ENCODE ENCFF356LFX recommended) |
| `BLACKLIST_MM39` | mm39 blacklist BED (Boyle lab recommended) |

#### Software

| Variable | Description |
|---|---|
| `PICARD_JAR` | Absolute path to `picard.jar` |
| `BEDGRAPH_TO_BIGWIG` | Absolute path to `bedGraphToBigWig` binary |
| `R_BIN` | Rscript binary (usually just `Rscript`) |

#### Thread counts

| Variable | Step | Recommended |
|---|---|---|
| `THREADS_PARALLEL_JOBS` | Max samples in parallel | 8 |
| `THREADS_ALIGN` | Bowtie2 `-p` per sample | 16 |
| `THREADS_SAMTOOLS` | samtools `-@` | 16 |
| `THREADS_TRIMGALORE` | TrimGalore `--cores` | 8 |
| `THREADS_FASTQC` | FastQC `-t` | 10 |
| `THREADS_BIGWIG` | bedtools / samtools for coverage | 16 |
| `THREADS_DEEPTOOLS` | deepTools workers for Step 10 QC | 16 |

> **`THREADS_DEEPTOOLS`** is required as of v3.1.0. Add it to any existing config:
> `echo "THREADS_DEEPTOOLS=16" >> config.conf`

Peak concurrent CPU usage: `THREADS_PARALLEL_JOBS × THREADS_ALIGN`

#### Cleanup

| Variable | Description | Default |
|---|---|---|
| `KEEP_INTERMEDIATE_BAMS` | Keep pre-dedup BAMs in `bams/` | `false` |
| `KEEP_TRIMMED_FASTQ` | Keep trimmed FASTQs in `trimmedFastq/` | `false` |

#### Deprecated variables (v3.0.x → v3.1.0)

The following variables were used by `run_chipqc.R` (Step 10 in v3.0.x) and are no longer
read by the pipeline. They may be safely removed from your config:

| Deprecated variable | Reason |
|---|---|
| `CHIPQC_ANNOTATION_HG38` | ChIPQC replaced by deepTools |
| `CHIPQC_ANNOTATION_MM39` | ChIPQC replaced by deepTools |
| `CHIPQC_BLACKLIST_HG38_RDS` | ChIPQC replaced by deepTools |
| `CHIPQC_BLACKLIST_MM39_RDS` | ChIPQC replaced by deepTools |
| `THREADS_CHIPQC` | Replaced by `THREADS_DEEPTOOLS` |

---

[← Installation](03_installation.md) | [Next: Running →](05_running.md)
