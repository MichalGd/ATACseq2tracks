# Script reference

This page describes every uploaded script included in `scripts/`. The script contents were preserved as uploaded.

| Script | Role | Main inputs | Main outputs |
|---|---|---|---|
| `fastq2tracks.2.1.sh` | Master end-to-end workflow | Raw FASTQ folder, output folder, max jobs, species | QC reports, trimmed FASTQ, BAMs, dedup BAMs, bedGraph, BigWig, reports |
| `fastq2tracks.2.1_blocked.sh` | Partial/restart workflow | Existing trimmed FASTQ and output tree, max jobs, species | Alignment, deduplication, coverage, reports |
| `fastqc_batch.1.0.sh` | Batch FastQC runner | FASTQ folder, FastQC output folder, max jobs | FastQC reports and batch logs |
| `trimgalore_batch.1.0.sh` | Batch Trim Galore runner | Raw paired FASTQ folder, trimmed output folder, max jobs | Trimmed paired FASTQ and logs |
| `bowtie2_human_batch.2.1.sh` | Human alignment batch wrapper | hg38 Bowtie2 index, trimmed FASTQ folder, BAM folder, max jobs | Human aligned sorted BAM files and logs |
| `bowtie2_mouse_batch.2.1.sh` | Mouse alignment batch wrapper | mm39 Bowtie2 index, trimmed FASTQ folder, BAM folder, max jobs | Mouse aligned sorted BAM files and logs |
| `bowtie2_dovetail_pairedEnd_Hsapiens.2.1.sh` | Single-sample human Bowtie2 mapper | Bowtie2 index, R1 trimmed FASTQ, output folder | Sorted BAM, BAI, Bowtie2 text log |
| `bowtie2_dovetail_pairedEnd_MMusculus.2.1.sh` | Single-sample mouse Bowtie2 mapper | Bowtie2 index, R1 trimmed FASTQ, output folder | Sorted BAM, BAI, Bowtie2 text log |
| `picard_deduplication_batch.2.1.sh` | Batch Picard deduplication wrapper | BAM folder, dedup output folder, max jobs | Deduplicated BAMs, indexes, Picard metrics |
| `picard_deduplication.2.1.sh` | Single-sample Picard deduplication | BAM folder, output folder, BAM filename | Deduplicated BAM, BAI, metrics text |
| `genomecoverage_batch.1.0.sh` | Batch coverage wrapper | Dedup BAM folder, genome assembly, max jobs, species | Raw bedGraph, normalized bedGraph, BigWig |
| `genomeCoverage_DNA_human.1.0.sh` | Human coverage generator | Dedup BAM, genome assembly | Human standard-chromosome bedGraph and BigWig |
| `genomeCoverage_DNA_mouse.1.0.sh` | Mouse coverage generator | Dedup BAM, genome assembly | Mouse standard-chromosome bedGraph and BigWig |
| `create_ucsc_tracks.1.0.sh` | UCSC track writer | Output folder, remote BigWig URL base | `ucsc_tracks.txt`, `bigwig_summary.txt` |
| `generate_pipeline_report.1.0.sh` | Pipeline execution report generator | Output folder, report name, format | `.Rmd`, `.html`, optional `.pdf` |
| `generate_multiqc_unified_report.1.0.sh` | Unified MultiQC HTML report generator | Pipeline output folder, report output folder, mode | Portable HTML report and optional assets folder |
| `readme.2.1.sh` | Change-log note | None | Not executable as shell code |

## Master workflow: `fastq2tracks.2.1.sh`

Usage:

```bash
./scripts/fastq2tracks.2.1.sh <inputFolder> <outputFolder> <max_jobs> <human|mouse>
```

Responsibilities:

1. Validate the species argument.
2. Create output folders.
3. Run raw FastQC and raw MultiQC.
4. Run Trim Galore.
5. Run FastQC and MultiQC on trimmed reads.
6. Select human or mouse Bowtie2 batch mapping.
7. Summarize alignment metrics with MultiQC.
8. Add read groups to BAM headers.
9. Run Picard duplicate removal.
10. Summarize deduplication metrics with MultiQC.
11. Generate raw and normalized coverage.
12. Move coverage files into final folders.
13. Generate UCSC track lines.
14. Remove selected large intermediate files.
15. Generate the pipeline report and unified MultiQC report.

## Partial workflow: `fastq2tracks.2.1_blocked.sh`

This variant has early steps commented out. It is useful when raw QC and trimming were already done and you want to restart from alignment or later. It assumes `trimmedFastq/` and the relevant parent output folders already exist.

## Batch-control design used across component scripts

Most batch scripts share the same pattern:

- Validate inputs.
- Create output and log folders.
- Redirect stdout and stderr to a timestamped main log.
- Trap SIGHUP, SIGINT, and SIGTERM for cleanup.
- Keep arrays of active job PIDs and sample names.
- Launch new jobs as soon as existing jobs finish.
- Skip samples when expected output files already exist.
- Write per-sample job logs and error logs.

## Lower-level mapper details

The human and mouse lower-level Bowtie2 scripts are structurally the same. They infer R2 from an R1 filename ending in either:

```text
_R1_001_val_1.fq.gz
_1_val_1.fq.gz
```

They run Bowtie2 with dovetail paired-end settings, convert SAM to BAM with Samtools, sort, index, move final outputs, and remove intermediate SAM/BAM files.

## Lower-level coverage details

The coverage scripts call `bedtools genomecov` twice: once for raw coverage and once for scaled coverage. They use `fetchChromSizes` and `bedGraphToBigWig` from Kent utilities.

## Report generators

`generate_pipeline_report.1.0.sh` writes an R Markdown file and renders it using R `rmarkdown`. It can request `html`, `pdf`, or `both` output. PDF output requires a LaTeX installation.

`generate_multiqc_unified_report.1.0.sh` reads MultiQC TSV output and copies or embeds PNG control plots. It was written to avoid portability issues with HTML image paths.
