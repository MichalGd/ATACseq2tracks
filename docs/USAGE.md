# Usage examples

## End-to-end human run

```bash
INPUT=/data/project/raw_fastq
OUTPUT=/data/project/fastq2tracks_human
MAX_JOBS=8
SPECIES=human

./scripts/fastq2tracks.2.1.sh "$INPUT" "$OUTPUT" "$MAX_JOBS" "$SPECIES"
```

## End-to-end mouse run

```bash
INPUT=/data/project/raw_fastq_mouse
OUTPUT=/data/project/fastq2tracks_mouse
MAX_JOBS=8
SPECIES=mouse

./scripts/fastq2tracks.2.1.sh "$INPUT" "$OUTPUT" "$MAX_JOBS" "$SPECIES"
```

## Run with reduced load

Use a smaller `max_jobs` value on a busy server.

```bash
./scripts/fastq2tracks.2.1.sh /data/fastq /data/results 2 human
```

## Run from a terminal that might disconnect

Most batch component scripts trap hangup and termination signals. For the master workflow, a conservative launch pattern is still useful:

```bash
nohup ./scripts/fastq2tracks.2.1.sh /data/fastq /data/results 4 human > fastq2tracks_master.log 2>&1 &

tail -f fastq2tracks_master.log
```

## Run only FastQC

```bash
./scripts/fastqc_batch.1.0.sh /data/fastq /data/results/fastQC/fastQC_unTrimmed 8
```

The FastQC script doubles the `max_jobs` parameter internally. A value of `8` means up to 16 concurrent FastQC jobs, each using 10 threads.

## Run only Trim Galore

```bash
./scripts/trimgalore_batch.1.0.sh /data/fastq /data/results/trimmedFastq 8
```

## Run only Bowtie2 alignment

Human:

```bash
./scripts/bowtie2_human_batch.2.1.sh   /path/to/bowtie2/hg38   /data/results/trimmedFastq   /data/results/bams   8
```

Mouse:

```bash
./scripts/bowtie2_mouse_batch.2.1.sh   /path/to/bowtie2/mm39   /data/results/trimmedFastq   /data/results/bams   8
```

## Run only Picard deduplication

The BAMs should already contain read-group tags. The master script adds them before calling Picard.

```bash
./scripts/picard_deduplication_batch.2.1.sh   /data/results/bams   /data/results/dedupBams   8
```

## Run only coverage generation

Human:

```bash
./scripts/genomecoverage_batch.1.0.sh /data/results/dedupBams hg38 8 human
```

Mouse:

```bash
./scripts/genomecoverage_batch.1.0.sh /data/results/dedupBams mm39 8 mouse
```

## Generate UCSC custom tracks

The URL should point to a web-accessible folder containing your BigWig files.

```bash
./scripts/create_ucsc_tracks.1.0.sh   /data/results   https://your.server.example.org/project/bigwig
```

## Generate reports after the pipeline

Pipeline execution report:

```bash
./scripts/generate_pipeline_report.1.0.sh   /data/results   pipeline_report_$(date +%Y%m%d)   html
```

Unified MultiQC report:

```bash
./scripts/generate_multiqc_unified_report.1.0.sh   /data/results   /data/results/reports/multiqc_summary_$(date +%Y%m%d)   selfcontained
```

## Restart from alignment onward

Use `fastq2tracks.2.1_blocked.sh` only after raw QC, trimming, and trimmed QC have already been completed, and after the expected output folders are already present.

```bash
./scripts/fastq2tracks.2.1_blocked.sh /data/fastq /data/results 4 human
```

Before running the blocked version, verify:

```bash
ls /data/results/trimmedFastq/*_val_1.fq.gz
ls /data/results/multiQC
```
