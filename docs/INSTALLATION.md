# Installation and configuration

## 1. Install command-line tools

The workflow expects a Linux or Unix-like shell environment with these tools available on `PATH` unless a script uses an absolute path.

Required tools:

- Bash
- FastQC
- Trim Galore and cutadapt
- Bowtie2
- Samtools
- Java
- Picard
- Bedtools
- MultiQC
- R
- Kent utilities: `fetchChromSizes` and `bedGraphToBigWig`

R packages used by the report scripts:

```r
install.packages(c(
  "rmarkdown",
  "knitr",
  "dplyr",
  "ggplot2",
  "readr",
  "lubridate",
  "stringr",
  "DT",
  "kableExtra",
  "base64enc"
))
```

`base64enc` is optional for the unified MultiQC report if the system `base64` command is available. `kableExtra` is needed by the PDF/LaTeX paths in the pipeline report script.

For PDF rendering, install TinyTeX or another LaTeX distribution:

```r
install.packages("tinytex")
tinytex::install_tinytex()
```

## 2. Prepare reference genomes

The master script currently uses these hard-coded Bowtie2 index prefixes:

```text
human: /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38
mouse: /home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39
```

Each prefix should point to a valid Bowtie2 index set. For example, an hg38 index prefix might have files such as:

```text
hg38.1.bt2
hg38.2.bt2
hg38.3.bt2
hg38.4.bt2
hg38.rev.1.bt2
hg38.rev.2.bt2
```

## 3. Decide how to handle hard-coded paths

The scripts were written for the original server and contain absolute paths. You have two deployment choices.

### Option A. Install into the expected paths

Use this only on the original server or a compatible clone of it.

```bash
sudo mkdir -p /home/micgdu/workflows/RNAseq/scripts
sudo cp scripts/*.sh /home/micgdu/workflows/RNAseq/scripts/
sudo chmod +x /home/micgdu/workflows/RNAseq/scripts/*.sh
```

Also make sure these paths exist or are adjusted:

```text
/home/micgdu/myenv/bin/activate
/home/micgdu/software/picard.jar
/home/micgdu/kentutils/fetchChromSizes
/home/micgdu/kentutils/bedGraphToBigWig
/home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38
/home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39
```

### Option B. Edit the paths for a new server

Search all hard-coded paths:

```bash
grep -R "/home/micgdu" -n scripts
```

Edit these common values:

| Current path or value | Meaning | Replace with |
|---|---|---|
| `/home/micgdu/workflows/RNAseq/scripts` | Location of this workflow's scripts | Absolute path to your cloned `scripts/` folder |
| `/home/micgdu/myenv/bin/activate` | Python environment containing MultiQC | Your environment activation script |
| `/home/micgdu/software/picard.jar` | Picard jar | Your Picard jar path or Picard command wrapper |
| `/home/micgdu/kentutils/fetchChromSizes` | UCSC/Kent chromosome-size utility | Your executable path |
| `/home/micgdu/kentutils/bedGraphToBigWig` | BigWig conversion utility | Your executable path |
| `/home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38` | Human Bowtie2 index prefix | Your hg38 index prefix |
| `/home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39` | Mouse Bowtie2 index prefix | Your mm39 index prefix |
| `http://your-server.com/data` | Placeholder BigWig hosting URL | Your real public URL |

## 4. Make scripts executable

```bash
chmod +x scripts/*.sh
```

`readme.2.1.sh` is not actually a shell script, but `chmod` is harmless. It is retained as an original change note.

## 5. Check syntax

```bash
bash tests/check_bash_syntax.sh
```

Expected result: all runnable scripts pass. `readme.2.1.sh` is reported as not runnable.

## 6. Test on a tiny dataset first

Before launching a full NovaSeq run, test with one small paired FASTQ sample and `max_jobs=1`.

```bash
./scripts/fastq2tracks.2.1.sh /data/test_fastq /data/test_fastq2tracks 1 human
```

Check these final outputs before scaling up:

```text
bigwig/*.bw
ucsc_tracks.txt
bigwig_summary.txt
reports/*.html
multiQC/*/*.html
```
