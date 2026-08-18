# Changelog

## v3.2.0 universal differential/annotation update - 2026-08-18

- Support any number and names of conditions and export all pairwise contrasts
  among conditions with at least two biological samples.
- Use all non-control biological samples for broad/narrow consensus construction
  and descriptive raw counts; exclude singleton conditions only from models.
- Fit one multi-condition model per method and peak type rather than refitting
  each pair independently.
- Add pair-level status and tested/significant/increased/decreased summaries in
  TSV and HTML, with explicit numerator/reference direction.
- Isolate broad/narrow failures, generate reports before propagating failures,
  and suppress automatic cleanup after differential-analysis failure.
- Add built-in GTF promoter/exon/intron/gene-body/nearest-TSS annotation and
  default-required genome-matched cCRE annotation for shared universes and
  results, with an explicit GTF-only opt-out.
- Add `utilities/prepare_encode4_hg38_ccre.sh` to download, strictly validate
  and record provenance for the native-GRCh38 expanded ENCODE4 registry.
- Enable the existing post-success intermediate cleanup policy by default;
  `ENABLE_AUTOMATIC_CLEANUP=false` remains the explicit opt-out and filtered
  quantitative BAMs remain retained.
- Retain two-condition root-level result names where practical and leave all
  upstream processing, track normalization and checkpoint architecture intact.

## v3.2.0 differential-accessibility update - 2026-08-13

- Change the DiffBind ATAC default to configurable `summits=100` while allowing
  `DIFFBIND_SUMMITS=200` to reproduce the previous window width.
- Add DESeq2ATAC as an independent DESeq2 analysis peer and extend it to run
  separate broad- and narrow-peak consensus analyses.
- Require configurable biological-sample support (default two), canonical
  chromosomes and blacklist exclusion in the DESeq2ATAC region universe.
- Count strict PE proper-pair fragments once or SE alignments once with
  `GenomicAlignments::summarizeOverlaps()`; no HTSeq dependency is added.
- Export compressed raw/normalized matrices and complete/FDR result tables,
  sample and library metadata, size factors, session information and PNG/PDF QC.
- Store broad and narrow outputs independently and add a peak-type comparison
  summary; document the exact support-filtered consensus construction and its
  difference from DiffBind summit recentering.
- Add a separate `step12a.done` checkpoint and attempt both differential modules
  before propagating either module's failure.
- Add focused static/self-tests, final-report integration and a legacy-method
  audit documenting deliberate modernizations and excluded study-specific code.

## v3.2.0 track-normalization update — 2026-08-06

- Replace the legacy `_RPM.bw` contract with fragment/read CPM `*_CPM.bw`; no CPM bedGraph is emitted.
- Generate reciprocal-DESeq2-consensus and DESeq2-robust-CPM coverage as both bigWig and bedGraph.
- Match `DESeq2::fpm(robust=TRUE)` through the consensus cohort's geometric mean of count-matrix column sums.
- Count paired-end fragments once and retain read-based semantics for single-end inputs.
- Use the same canonical autosome/X/Y universe for consensus construction, consensus counting and every track.
- Export per-sample normalization counts, factors, cohort constants and applied scales.
- Add focused normalization/output-contract tests while preserving the existing checkpoint architecture.

### Explicit single-end variant — 2026-08-06

- Add `SE_SIGNAL_MODE="read"` as the explicit supported SE normalization mode.
- Use each filtered SE alignment once, without artificial read extension, for coverage, CPM denominators and consensus-peak counts.
- Reject mixed PE/SE samplesheets so fragment- and read-based units cannot enter one normalization cohort.
- Add an SE ATAC-seq samplesheet example plus mocked SE execution and validation tests.
- Keep MACS3 SE shift/extension limited to peak calling; Tn5 insertion-site tracks remain outside this update.

Deployment and resume instructions are in
`docs/v3.2.0_TRACK_NORMALIZATION_UPDATE_2026-08-06.md`.

## v3.2.0 test-run hotfix — 2026-08-04

The public version remains **3.2.0** while the first complete server run is
used to find and correct integration defects. This deliberately does not open
v3.2.1 yet.

- Read the DESeq2 consensus `size_factor` by its TSV header instead of the
  former positional field, which incorrectly returned the genome (`mm39`).
- Read Picard `PERCENT_DUPLICATION` by its tab-delimited header instead of a
  whitespace-shifted field that was actually the estimated library size.
- Include `logs/picard/` in the dedicated Step 5 MultiQC input paths.
- Ensure BAM indexes exist before `idxstats`, and validate the BAM, scaling
  factor, chromosome selection, `bamCoverage` exit status, and bigWig output.
- Interpret broadPeak files as generic `bed` input for DiffBind; `broad` is not
  a supported DiffBind peak-caller identifier.
- Skip a header-only narrow or broad DiffBind companion sheet, while still
  failing if no runnable differential analysis exists.
- At the 2026-08-04 hotfix stage, require exactly two conditions with at least
  two biological replicates per condition. This restriction was superseded by
  the universal multi-condition update of 2026-08-18 above.
- Use all normalized sites for DiffBind PCA and heatmap QC so those plots are
  not contingent on finding an FDR-significant site.
- Export both all tested DiffBind sites and the FDR 0.05 subset, and explicitly
  print lattice/ggplot objects into their PNG devices.
- Extend preflight checks to the Python and R packages used downstream.
- Add regression fixtures for both header-aware metric parsers.

Deployment and resume instructions are in
`docs/v3.2.0_TEST_RUN_HOTFIX_2026-08-04.md`.

### Windows-input safety follow-up — 2026-08-04

- Detect CRLF, bare CR, UTF-8/UTF-16/UTF-32 BOM encodings and a trailing DOS
  Ctrl-Z marker before sourcing the configuration or parsing the samplesheet.
- Create a timestamped, byte-for-byte backup beside every affected input before
  modifying it.
- Atomically replace affected inputs with UTF-8/LF text and use those corrected
  files for the current run.
- Leave clean UTF-8/LF files unchanged and create no unnecessary backup.
- Reject unknown encodings without modifying or backing up the input.
- Add regression coverage for conversion, backup integrity and idempotence.

Deployment instructions are in
`docs/v3.2.0_WINDOWS_INPUT_HOTFIX_2026-08-04.md`.

---

## v3.2.0 — 2026-08-02 (release candidate)

### Scientific and output changes

- Generate filtered-BAM CPM/RPM bigWigs by default for UCSC and IGV.
- Define the consensus peak set using a configurable minimum number of distinct biological samples; default two.
- Estimate DESeq2 size factors from consensus-peak counts and generate separate DESeq2-consensus-scaled bigWigs without an additional RPM divisor.
- Use paired-end fragments (`BAMPE`) for paired-end MACS3 and explicit ATAC shift/extension parameters for single-end data.
- Add default `ataqv` ENCODE-style TSS enrichment, compressed JSON metrics and an optional interactive viewer.
- Add full-scan paired-end nucleosome-periodicity metrics plus 300-dpi PNG and vector PDF plots.
- Retain DiffBind differential accessibility analysis and keep raw integer counts separate from visualization tracks.
- Preserve biological replicates in DiffBind preparation using `sample_id + replicate` keys and force the DESeq2 analysis method explicitly.
- Add optional HOMER peak annotation and motif enrichment hooks, disabled by default.

### Reliability and consistency

- Propagate failed child jobs instead of reporting silent stage success.
- Distinguish raw and trimmed FastQC inputs.
- Add MAPQ/SAM-flag, mitochondrial, blacklist and pair-safe BAM filtering with attrition metrics.
- Reject mixed-genome runs and empty peak results by default.
- Replace empty checkpoints with signatures covering the samplesheet, configuration and workflow version.
- Centralize the runtime version in `VERSION`.
- Remove ChIPQC from the environment and preflight requirements; require DESeq2 explicitly.
- Disable automatic cleanup by default and document the retained/deleted file policy.

### Upgrade notes

Review `docs/update_3.2.0/APPLY_UPDATE.md` and merge site-specific paths into the new configuration template. Remove the Step 10 checkpoint before regenerating QC for an existing run.

---

---

## v3.1.0 — 2026-06-05 (previous)

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
