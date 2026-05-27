# Reporting and QC

## QC checkpoints

The workflow has four main MultiQC checkpoints:

| Checkpoint | Folder | Purpose |
|---|---|---|
| Raw FASTQ QC | `multiQC/multiQC_unTrimmed/` | Inspect raw input read quality, GC content, sequence duplication, adapters, overrepresented sequences, and basic read metrics |
| Trimmed FASTQ QC | `multiQC/multiQC_trimmed/` | Confirm adapter removal and trimming effects |
| Alignment QC | `multiQC/multiQC_alignments/` | Summarize Bowtie2 alignment logs and mapping rates |
| Deduplication QC | `multiQC/multiQC_deduplication/` | Summarize Picard duplicate-removal metrics |

The PowerPoint deck included in this package highlights many of these QC areas, including FastQC/MultiQC summaries, unique and duplicate sequence counts, per-sequence GC content, overrepresented sequences, adapter content, status checks, alignment percentages, paired-end alignment scores, and deduplication statistics.

## Pipeline execution report

`scripts/generate_pipeline_report.1.0.sh` is the newly added pipeline-report generator. It creates an R Markdown file in the output folder, then renders it to `html`, `pdf`, or `both`.

Usage:

```bash
./scripts/generate_pipeline_report.1.0.sh <output_folder> <report_name> [html|pdf|both]
```

Example:

```bash
./scripts/generate_pipeline_report.1.0.sh /data/results pipeline_report_$(date +%Y%m%d) html
```

The report includes:

- System information.
- Hardware summary.
- Software version checks.
- Script-level resource usage tables.
- Theoretical resource use estimates.
- Output directory summaries.
- Batch log summaries.
- File inventory and storage summaries.

PDF mode requires LaTeX. The master workflow calls this generator in HTML mode.

## Unified MultiQC report

`scripts/generate_multiqc_unified_report.1.0.sh` creates a portable HTML report that combines key per-sample MultiQC tables and plots.

Usage:

```bash
./scripts/generate_multiqc_unified_report.1.0.sh <InFolder> <OutFolder> [selfcontained|assets]
```

Mode choices:

- `selfcontained`: embeds copied PNG plots as base64 data URIs in a single HTML report.
- `assets`: writes HTML plus an adjacent assets folder, which must be kept together.

Example:

```bash
./scripts/generate_multiqc_unified_report.1.0.sh   /data/results   /data/results/reports/multiqc_summary_$(date +%Y%m%d)   selfcontained
```

## Presentation deck

The source deck is stored at:

```text
presentations/fastq2tracks_LAbMeeting_29jan2026.pptx
```

An extracted text companion is available at:

```text
docs/PRESENTATION_TEXT_EXTRACT.md
```

The deck appears to document an earlier version of the workflow and includes screenshots or examples of MultiQC reports, BigWig summaries, UCSC tracks, HTML report pages, and RPubs output. Current script filenames differ in places, so this repository treats the deck as historical workflow context rather than the definitive source of execution commands.
