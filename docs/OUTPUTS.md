# Outputs

The workflow writes results under the output folder supplied to the master script.

## Output tree

```text
<outputFolder>/
├── reports/
├── fastQC/
│   ├── fastQC_unTrimmed/
│   └── fastQC_trimmed/
├── multiQC/
│   ├── multiQC_unTrimmed/
│   ├── multiQC_trimmed/
│   ├── multiQC_alignments/
│   └── multiQC_deduplication/
├── trimmedFastq/
├── bams/
├── dedupBams/
├── bedGraph/
├── NormBedGraph/
├── bigwig/
├── ucsc_tracks.txt
└── bigwig_summary.txt
```

## Main file types

| File type | Created by | Meaning |
|---|---|---|
| `*_fastqc.html`, `*_fastqc.zip` | FastQC | Per-FASTQ quality control reports |
| `multiqc_*.html` | MultiQC | Aggregated QC reports |
| `*_val_1.fq.gz`, `*_val_2.fq.gz` | Trim Galore | Trimmed paired-end FASTQ files |
| `*.sorted_stChr.bam`, `*.sorted_stChr.bai` | Bowtie2 and Samtools | Sorted alignment BAM and index |
| `*_stChrH.bam` | Samtools | BAM after read-group insertion |
| `*_dedup.bam`, `*_dedup.bam.bai` | Picard | Duplicate-removed BAM and index |
| `*_dedup_rep.txt` | Picard | Deduplication metrics |
| `*.bedGraph.gz` | Bedtools | Raw coverage bedGraph |
| `*_Snorm.bedGraph.gz` | Bedtools and AWK | RPM-normalized sorted bedGraph |
| `*_Snorm.bw` | Kent utilities | BigWig coverage track |
| `ucsc_tracks.txt` | `create_ucsc_tracks.1.0.sh` | UCSC custom-track definitions |
| `bigwig_summary.txt` | `create_ucsc_tracks.1.0.sh` | BigWig inventory |
| `pipeline_report_*.html` | `generate_pipeline_report.1.0.sh` | Pipeline execution report |
| `multiqc_unified_report_*.html` | `generate_multiqc_unified_report.1.0.sh` | Combined QC summary |

## Cleanup behavior

The master script removes selected large intermediate files near the end:

```bash
rm "$2"/bams/*.bam "$2"/bams/*.bai
rm $2/trimmedFastq/*.gz
```

This keeps final outputs smaller, but it means the original alignment BAMs and trimmed FASTQ files may not remain after a successful complete run. Comment out those lines if you need to preserve intermediates for debugging or reanalysis.

## BigWig and UCSC custom tracks

To view BigWig files in UCSC Genome Browser:

1. Host the final `.bw` files at a URL accessible to the browser.
2. Run `create_ucsc_tracks.1.0.sh` with the URL base.
3. Open `ucsc_tracks.txt`.
4. Copy the generated track lines into the UCSC custom-track input box.

The track generator uses generic styling: black bar plots, fixed y-axis limits, full visibility, and one track per BigWig file.
