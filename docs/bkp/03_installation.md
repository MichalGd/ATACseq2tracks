# 03 — Installation

[← Quick start](02_quickstart.md) | [Next: Input files →](04_inputs.md)

---

## Conda environment

The conda environment file `environment.yml` installs all command-line dependencies.

```bash
conda env create -f environment.yml
conda activate fastq2tracks
```

| Tool | Minimum version | Role |
|---|---|---|
| bowtie2 | 2.5 | Alignment |
| samtools | 1.18 | BAM sorting, indexing, merging |
| bedtools | 2.31 | Blacklist filtering, coverage |
| trim_galore | 2.2.0 | Adapter trimming (SE + PE) |
| fastqc | 0.12 | Read quality control |
| macs2 | 2.2.9 | Peak calling |
| multiqc | 1.20 | Aggregated QC HTML reports |
| picard | 3.1 (`.jar`) | Duplicate marking |
| bedGraphToBigWig | UCSC kent | bigwig conversion |
| python3 | 3.9+ | Samplesheet validation |

> **Note on trim_galore:** Version 2.2.0 ("Oxidized Edition") uses `--output_dir` instead of `-o`. The pipeline scripts are written for this version. Earlier versions will fail at step 2.

---

## R packages

R and Bioconductor are **not** managed by conda and must be installed separately.

```r
# Install Bioconductor manager if needed
if (!require("BiocManager")) install.packages("BiocManager")

# Core packages
BiocManager::install(c(
    "ChIPQC",
    "DiffBind",
    "BiocParallel",
    "GenomicAlignments",
    "rtracklayer"
))

# Utility packages
install.packages(c("ggplot2", "dplyr"))
```

Minimum R version: **4.3**

---

## Picard

Download `picard.jar` from [Broad Institute](https://github.com/broadinstitute/picard/releases) and set its path in `config.conf`:

```bash
PICARD_JAR="/path/to/picard.jar"
```

---

## Reference files

The pipeline requires pre-built reference files. Set their paths in `config.conf`.

### Bowtie2 indices

```bash
# hg38
bowtie2-build /path/to/hg38.fa /path/to/indices/hg38/hg38

# mm39
bowtie2-build /path/to/mm39.fa /path/to/indices/mm39/mm39
```

Set in config:
```bash
INDEX_HG38="/path/to/indices/hg38/hg38"
INDEX_MM39="/path/to/indices/mm39/mm39"
```

### Chromosome sizes

```bash
fetchChromSizes hg38 > /path/to/hs38n.chrom.sizes
fetchChromSizes mm39 > /path/to/mm39n.chrom.sizes
```

(`fetchChromSizes` is part of the UCSC kent tools, included in the conda environment)

### Blacklist BED files

| Genome | Source | Download |
|---|---|---|
| hg38 | ENCODE ENCFF356LFX | `wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz` |
| mm39 | Boyle lab | [github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist) |

### ChIPQC annotation and blacklist RDS objects

These are R objects required by `run_chipqc.R`. See [Reference file preparation](10_reference_files.md) for the complete build script.

```bash
CHIPQC_ANNOTATION_HG38="/path/to/anno_hg38_chipqc.rds"
CHIPQC_ANNOTATION_MM39="/path/to/anno_mm39_chipqc.rds"
CHIPQC_BLACKLIST_HG38_RDS="/path/to/blacklist_hg38.rds"
CHIPQC_BLACKLIST_MM39_RDS="/path/to/blacklist_mm39.rds"
```

---

## Make scripts executable

```bash
cd /path/to/fastq2tracks
chmod +x fastq2tracks.sh scripts/*.sh scripts/*.py scripts/*.R
```

---

## Verify installation

```bash
conda activate fastq2tracks
bash fastq2tracks/scripts/smoke_test.sh \
    /path/to/my_project/config/samplesheet.csv \
    /path/to/my_project/config/config.conf
```

A successful pre-flight check prints:
```
============================================================
 RESULT: 53 OK | 0 WARNINGS | 0 FAILURES
============================================================
Pre-flight PASSED -- ready to run.
```

---

[← Quick start](02_quickstart.md) | [Next: Input files →](04_inputs.md)
