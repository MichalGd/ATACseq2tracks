# Changelog

---

## v3.1.0 — 2026-06-05 (current)

### Breaking changes
- **`chipqc_annotation` samplesheet column removed.** Column 17 (`chipqc_annotation`) no longer
  exists. The samplesheet now has 17 columns (previously 18). Update your samplesheets accordingly.
- **`CHIPQC_ANNOTATION_*` and `CHIPQC_BLACKLIST_*_RDS` config variables deprecated.**
  These variables are no longer read by the main pipeline. They may be removed from `config.conf`
  safely. Legacy build instructions are preserved in `docs/10_reference_files.md`.

### New features
- **Post-alignment QC module** (`scripts/post_alignment_qc_batch.sh` + `scripts/plot_chrom_coverage.py`)
  replaces `run_chipqc.R` as Step 10. Based on deepTools and standard command-line tools.
  No R or Bioconductor required for QC.
- **New outputs:** `qc_post_alignment/` directory with FRiP tables, fingerprint Lorenz curves,
  Spearman correlation heatmaps (bins + peaks), PCA plots (bins + peaks), signal heatmap and
  average profile over consensus peaks, and per-sample chromosome karyogram plots.
- **Chromosome karyogram plots** matching the ChIPQC "ChIP Peaks over Chromosomes" panel style —
  one row per chromosome with signal bars, cytogenetic order, and shared position axis.
- **Consensus peak set construction** — merged BED from all per-replicate peak files used for
  cross-sample FRiP and signal comparisons.
- **Zero-peak safety** — samples with no peaks are flagged `NO_PEAKS` in the summary table and
  are never silently dropped from reports.
- **Quality warning flags** — `qc_warnings.tsv` lists samples with `LOW_READS`, `HIGH_DUPLICATION`,
  `HIGH_MITO`, `NO_PEAKS`, `LOW_FRIP_NARROW`, etc.
- **New config variable** `THREADS_DEEPTOOLS` required. Add to `config.conf`:
  `THREADS_DEEPTOOLS=8`
- **ChIPmentation** explicitly added as supported assay type.
- **Documentation:** `docs/12_post_alignment_qc.md` added (deepTools QC module reference).
  All other `docs/` pages updated to reflect QC module change.

### Upgrade from v3.0.4
```bash
# 1. Remove chipqc_annotation column (col 17) from samplesheet
# 2. Add THREADS_DEEPTOOLS to config.conf
echo "THREADS_DEEPTOOLS=8" >> config.conf
# 3. Install deepTools if not in conda environment
mamba install -c bioconda deeptools>=3.5
# 4. Copy new scripts
cp scripts/post_alignment_qc_batch.sh  /path/to/ATACseq2tracks/scripts/
cp scripts/plot_chrom_coverage.py       /path/to/ATACseq2tracks/scripts/
chmod +x /path/to/ATACseq2tracks/scripts/post_alignment_qc_batch.sh
# 5. Remove step10 checkpoint to rerun QC
rm /path/to/analysis/.checkpoints/step10.done
```

---

## v3.0.4 — 2026-05-27

- Fixed `THREADS_BOWTIE2` → `THREADS_ALIGN` variable name throughout
- Fixed Trim Galore `-o` → `--output_dir` (v2.2.0 Oxidized Edition)
- Fixed Picard BAM glob (`*.sorted_stChr.bam` → `*.bam`)
- Fixed `_bioRN` suffix in BAM lookup keys across all downstream scripts (steps 6–13)
- Fixed missing closing `}` braces in `blacklist_filter_batch.sh`
- Fixed `ctrl_rep` undefined variable in `macs2_batch.sh`
- Fixed bare `replicate` variable in R scripts
- Replaced MACS2 with MACS3 v3.0.4 (`pip install macs3`; identical CLI, Python 3.11 compatible)

---

## v3.0.0 — 2026-05-19

- Initial unified pipeline release
- SE and PE read support, blacklist filtering, replicate and control tracking
- MACS2 narrow + broad peak calling, ChIPQC module, DiffBind samplesheet preparation
- Checkpoint system, `environment.yml`, samplesheet validator
- Unified `config.conf` replacing per-genome hardcoded scripts

---

## v2.1 (archived)

Original version: PE reads only, separate genome-specific batch scripts for human and mouse,
hardcoded paths. Archived to `scripts/legacy/`.

### Historical notes (from `readme.2.1.sh`)
- Picard Java heap reduced from 128 GB to 32 GB
- Bowtie2 threads changed from 15 to 12
- samtools commands updated to use multithreading
- Human and mouse Bowtie2 scripts separated by filename
