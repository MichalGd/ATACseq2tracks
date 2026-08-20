# 02 — Quick start

[← Overview](01_overview.md) | [Next: Installation →](03_installation.md)

---

This page gets you from zero to a running pipeline in five steps.
For full details see [Installation](03_installation.md) and [Input files](04_inputs.md).

---

## Step 1 — Install

```bash
git clone https://github.com/MichalGd/ATACseq2tracks.git ATACseq2tracks
cd ATACseq2tracks
conda env create -f environment.yml
conda activate ATACseq2tracks
chmod +x atacseq2tracks.sh scripts/*.sh scripts/*.py scripts/*.R
```

> **R/Bioconductor packages** for DiffBind, DESeq2ATAC and consensus-track
> normalization must be installed; see [Installation](03_installation.md#r-packages).

---

## Step 2 — Create a project directory

```bash
mkdir -p /path/to/my_project/config
cp config/config.conf              /path/to/my_project/config/config.conf
cp config/samplesheet_template.csv /path/to/my_project/config/samplesheet.csv
```

---

## Step 3 — Edit the config

Open `/path/to/my_project/config/config.conf` and set at minimum:

```bash
SAMPLESHEET="/path/to/my_project/config/samplesheet.csv"
OUTPUT_DIR="/path/to/my_project/analysis/"

INDEX_HG38="/path/to/bowtie2_index/hg38"
CHROM_SIZES_HUMAN="/path/to/hs38n.chrom.sizes"
BLACKLIST_HG38="/path/to/blacklist_hg38_ENCFF356LFX.bed"
GTF_HUMAN="/path/to/gencode.hg38.annotation.gtf"
CCRE_BED_HG38="/path/to/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
PICARD_JAR="/path/to/picard.jar"
BEDGRAPH_TO_BIGWIG="/path/to/bedGraphToBigWig"

# deepTools QC threads (Step 10)
THREADS_DEEPTOOLS=8
```

> ChIPQC RDS annotation variables are unsupported and should be removed from current configuration files. The retired `run_chipqc.R` implementation is available only through older Git history.

---

## Step 4 — Fill in the samplesheet

Edit `/path/to/my_project/config/samplesheet.csv`. Minimum required columns:

```
sample_id, fastq_1, fastq_2, layout, genome, assay, factor, condition, treatment,
cell_type, replicate, tech_replicate, is_control, control_id, macs2_mode,
blacklist, output_prefix
```

See [Input files](04_inputs.md) for the full column reference and example rows.

Validate before running:

```bash
python3 /path/to/ATACseq2tracks/scripts/validate_samplesheet.py \
    /path/to/my_project/config/samplesheet.csv
```

---

## Step 5 — Run

```bash
conda activate ATACseq2tracks

# Run from the PARENT folder of ATACseq2tracks/
cd /path/to/parent_folder

nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    > /path/to/my_project/run.log 2>&1 &

echo $! > /path/to/my_project/run.pid
tail -f /path/to/my_project/run.log
```

---

## Expected output locations

| Result | Path |
|---|---|
| CPM bigWigs | `analysis/bigwig/` |
| DESeq2 consensus bigWig/bedGraph tracks | `analysis/bigwig_deseq2_consensus/` |
| DESeq2 robust-CPM permissive/intermediate/stringent tracks | `analysis/bigwig_deseq2_robust_cpm/<policy>/` |
| Merged tracks | `analysis/bigwig_merged/` |
| Peaks (narrow + broad) | `analysis/peaks/per_replicate/` |
| Post-alignment QC | `analysis/qc_post_alignment/` |
| DiffBind samplesheets | `analysis/diffbind/` |
| DiffBind results | `analysis/diffbind_results/` |
| DESeq2ATAC broad/narrow analyses | `analysis/deseq2atac/{broad,narrow}/` |
| Consensus peak annotations | DiffBind analysis roots and `analysis/deseq2atac/{broad,narrow}/` |
| Pair-level differential summary | `analysis/reports/differential_accessibility_summary.{tsv,html}` |
| HTML pipeline report | `analysis/reports/pipeline_report_<date>.html` |

---

## If the run fails

Check the log for the failed step, then:

```bash
# Remove only the failed step's checkpoint and re-run
rm /path/to/my_project/analysis/.checkpoints/step<N>.done
nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    >> /path/to/my_project/run_resume.log 2>&1 &
```

See [Troubleshooting](09_troubleshooting.md) for common errors.

---

[← Overview](01_overview.md) | [Next: Installation →](03_installation.md)
