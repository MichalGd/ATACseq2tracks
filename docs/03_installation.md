# 03 — Installation

[← Quick start](02_quickstart.md) | [Next: Input files →](04_inputs.md)

---

## Conda environment

The conda environment file `environment.yml` installs all command-line dependencies.

```bash
conda env create -f environment.yml
conda activate ATACseq2tracks
```

| Tool | Minimum version | Role |
|---|---|---|
| bowtie2 | 2.5 | Alignment |
| samtools | 1.18 | BAM sorting, indexing, merging |
| bedtools | 2.31 | Blacklist filtering, coverage, FRiP |
| trim_galore | 0.6.10 | Adapter trimming (SE + PE) |
| fastqc | 0.12 | Read quality control |
| macs3 | 3.0.4 | Peak calling (narrow + broad) |
| multiqc | 1.20 | Aggregated QC HTML reports |
| picard | 3.1 (`.jar`) | Duplicate marking |
| bedGraphToBigWig | UCSC kent | bigwig conversion |
| deeptools | 3.5+ | Post-alignment QC (Step 10) |
| ataqv | 1.3+ | ATAC TSS enrichment and fragment QC (Step 10) |
| python3 | 3.9+ | Samplesheet validation, karyogram plots |
| matplotlib | 3.7+ | Chromosome coverage plots (Step 10) |
| numpy | 1.24+ | Numerical arrays for karyogram plots |
| pandas | 2.0+ | Table handling for karyogram plots |

> **Note on Trim Galore:** the Conda environment requires version 0.6.10 or newer.

> **Note on MACS3:** The pipeline calls `macs3` (not `macs2`). MACS3 v3.0.4+ is the
> Python 3.11-compatible successor to MACS2 with an identical command-line interface.
> Install via `pip install macs3` if not included in `environment.yml`.

---

## R packages

R and Bioconductor are required for DESeq2 consensus-peak size factors in Step 10 and for DiffBind in Steps 11–12. ChIPQC is not required.

```r
# Install Bioconductor manager if needed
if (!require("BiocManager")) install.packages("BiocManager")

# DiffBind and its dependencies
BiocManager::install(c(
    "DESeq2",
    "DiffBind",
    "BiocParallel",
    "GenomicAlignments",
    "rtracklayer"
))

# Utility packages
install.packages(c("ggplot2", "dplyr"))
```

Minimum R version: **4.3**

> Disabling DiffBind does not remove the R requirement when `GENERATE_DESEQ2_CONSENSUS_TRACKS=true`.

---

## Picard

Download `picard.jar` from [Broad Institute](https://github.com/broadinstitute/picard/releases)
and set its path in `config.conf`:

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

See [Reference file preparation](10_reference_files.md) for the full download commands.

> ChIPQC and its RDS annotation objects are not part of v3.2.0. Use an older Git tag if historical ChIPQC reproduction is required.

---

## Shared Linux server installation

On a multi-user Linux server, install the repository into a shared path such as `/opt/ATACseq2tracks` or `/usr/local/ATACseq2tracks`. Make the code tree readable and executable by all users, while keeping each user’s project config and output directories separate from the shared code installation.

```bash
sudo mkdir -p /opt/ATACseq2tracks
sudo chown -R root:your_group /opt/ATACseq2tracks
sudo chmod -R a+rX /opt/ATACseq2tracks
```

Each user should create their own project directory elsewhere, copy `config/config.conf` and `config/samplesheet_template.csv` into it, and run the pipeline from the shared installation.

## Make scripts executable

```bash
cd /path/to/ATACseq2tracks
chmod +x atacseq2tracks.sh scripts/*.sh scripts/*.py scripts/*.R
```

---

## Verify installation

```bash
conda activate ATACseq2tracks
bash ATACseq2tracks/scripts/smoke_test.sh \
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
