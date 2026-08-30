# 04 — Input files

[← Installation](03_installation.md) | [Next: Running →](05_running.md)

---

## Samplesheet

The samplesheet is a **comma-separated CSV** with one row per sequencing library. It is the
primary driver of the entire pipeline — every step reads it to know which samples to process,
which genome to use, how to pair IPs with controls, and what type of peaks to call.

> **v3.1.0 change:** The `chipqc_annotation` column (column 17 in v3.0.x) has been **removed**.
> The samplesheet now has **17 columns**. If you are upgrading from v3.0.x, remove this column
> from your existing samplesheets before running v3.1.0.

### Column reference

| # | Column | Type | Required | Description |
|---|---|---|---|---|
| 1 | `sample_id` | string | ✓ | Unique identifier for this library. Used as the key for all output file names. |
| 2 | `fastq_1` | path | ✓ | Absolute path to R1 FASTQ (`.fastq.gz`). |
| 3 | `fastq_2` | path | PE only | Absolute path to R2 FASTQ. Leave empty for SE samples. |
| 4 | `layout` | `SE`\|`PE` | ✓ | Single-end or paired-end. Every row in one run must use the same layout. |
| 5 | `genome` | `hg38`\|`mm39` | ✓ | Reference genome for alignment and peak calling. |
| 6 | `assay` | string | ✓ | Assay type, e.g. `ChIP-seq`, `ATAC-seq`, `CUT&RUN`, `ChIPmentation`. |
| 7 | `factor` | string | ✓ | Antibody target or histone mark, e.g. `H3K27ac`, `CTCF`, `Input`. |
| 8 | `condition` | string | ✓ | Biological condition, e.g. `day0`, `treated`, `rest`. |
| 9 | `treatment` | string | ✓ | Fixation protocol or treatment, e.g. `FA`, `DSG_FA`, `none`. |
| 10 | `cell_type` | string | ✓ | Cell line or tissue, e.g. `NHEK`, `LymphocyteT`. |
| 11 | `replicate` | integer | ✓ | Biological replicate number: `1`, `2`, `3`, … |
| 12 | `tech_replicate` | integer | ✓ | Technical replicate (sequencing run) number. Use `1` if the library was sequenced once. |
| 13 | `is_control` | `TRUE`\|`FALSE` | ✓ | `TRUE` for input DNA or IgG controls; `FALSE` for IP samples. |
| 14 | `control_id` | string | IP only | The `sample_id` of the matched control. Must match exactly. Leave empty for control rows. |
| 15 | `macs2_mode` | `both`\|`narrow`\|`broad`\|`none` | ✓ | Peak calling mode. Use `none` for all control samples. |
| 16 | `blacklist` | path | ✓ | Absolute path to the blacklist BED file for this genome. |
| 17 | `output_prefix` | string | ✓ | Prefix for output files. Typically same as `sample_id`. |

For the optional v4.2.0 dm6 calibration branch, add all three columns below.
Partial declarations are rejected; legacy sheets remain valid when the branch
is disabled.

| Column | Required in spike-in mode | Description |
|---|---|---|
| `spikein_genome` | yes | `dm6` |
| `spikein_to_host_ratio` | yes | Positive relative dm6:host amount added on the same basis for every sample; use `1` when equal |
| `spikein_stage` | yes | `pre_tagmentation_nuclei` |

### Rules and conventions

**Controls**
- Set `is_control=TRUE` and `macs2_mode=none`.
- Leave `control_id` empty.
- Controls still go through all preprocessing steps (trimming → alignment → dedup → filtering → tracks).

**IP samples**
- Set `is_control=FALSE`.
- `control_id` must exactly match the `sample_id` of the appropriate control.
- `macs2_mode=both` runs MACS3 in both narrow and broad mode. It is required for
  every non-control sample when `RUN_DESEQ2ATAC=true`, because DESeq2ATAC runs
  separate broad- and narrow-consensus analyses.

**Technical replicates**
- Add one row per sequencing run with the same `sample_id`, `replicate`, and different
  `tech_replicate` (1, 2, …).
- The pipeline merges their FASTQs before trimming.
- v4.3.0 writes the exact ordered merge mapping to
  `metadata/technical_merge_audit.tsv`; review it during `--plan`.
- Technical rows contribute one downstream biological library and must never be
  used to inflate biological replicate counts.

**Sequencing layout**
- Use a PE-only or SE-only samplesheet; mixed PE/SE runs are rejected during validation.
- PE normalization counts properly paired fragments once.
- SE normalization uses one filtered alignment per observation with `SE_SIGNAL_MODE="read"`; leave `fastq_2` empty.
- Use separate output directories because PE fragment signal and SE read signal are different units.

**Drosophila spike-in**
- Add dm6 nuclei/cells before tagmentation; a post-library mixture cannot
  calibrate upstream technical losses.
- Technical-replicate rows for one biological key must use identical spike-in
  metadata.
- Use `config/samplesheet_example_atac_spikein.csv` and enable
  `GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS=true`.

**Batch metadata (optional)**
- You may add an optional `batch` column to the samplesheet.
- The pipeline preserves this field in DiffBind sample sheets if present.
- Use `batch` when you have known technical batches, library prep blocks, or sequencing lanes.

**Conditions and differential eligibility**
- Any number and names of conditions are accepted; no condition list is hard-coded.
- All non-control biological samples contribute to broad/narrow consensus construction.
- A condition needs at least two biological samples to enter differential models.
- Conditions represented once still receive trimming, alignment, filtering, peaks,
  coverage tracks, QC, and consensus participation; only their statistical
  modeling and contrasts are skipped.

**Genome build**
- Use one genome build per samplesheet; mixed hg38/mm39 runs are rejected.
- Use a separate configuration and output directory for each genome build.

**macs2_mode guidance**

| Mark / Assay | Recommended mode |
|---|---|
| H3K27ac, H3K4me3, H3K4me1, p63, CTCF, SATB1 | `both` |
| H3K27me3, H3K9me3, H3K36me3 | `both` (broad peaks are primary) |
| ATAC-seq | `both` with DESeq2ATAC enabled; otherwise `narrow` |
| CUT&RUN, CUT&TAG, ChIPmentation | `narrow` (or `both`) |
| Input / IgG | `none` |

---

### Example samplesheets

The repository provides layout-specific ATAC-seq examples:

- `config/samplesheet_example_atac.csv` — paired-end;
- `config/samplesheet_example_atac_se.csv` — single-end.

Minimal SE rows look like this:

```csv
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
ATAC_SE_WT_R1,/data/ATAC_SE_WT_R1.fastq.gz,,SE,hg38,ATAC-seq,accessibility,WT,none,cell,1,1,FALSE,,narrow,/ref/blacklist_hg38.bed,ATAC_SE_WT_R1
ATAC_SE_WT_R2,/data/ATAC_SE_WT_R2.fastq.gz,,SE,hg38,ATAC-seq,accessibility,WT,none,cell,2,1,FALSE,,narrow,/ref/blacklist_hg38.bed,ATAC_SE_WT_R2
```

> Use `python3 scripts/validate_samplesheet.py samplesheet.csv` to check for errors before running.

---

## Configuration file

`config/config.conf` is a bash script that is `source`d by the master script and all subscripts.
Copy the template and edit for your server.

### All configuration parameters

#### Project paths

| Variable | Description |
|---|---|
| `SAMPLESHEET` | Absolute path to the samplesheet CSV |
| `OUTPUT_DIR` | Absolute path to the root output directory (created if absent) |

#### Reference genomes and indices

| Variable | Description |
|---|---|
| `INDEX_HG38` | Bowtie2 index prefix for hg38 |
| `INDEX_MM39` | Bowtie2 index prefix for mm39 |
| `INDEX_HG38_DM6` | Bowtie2 composite index prefix for competitively aligned hg38+dm6 reads |
| `INDEX_MM39_DM6` | Bowtie2 composite index prefix for competitively aligned mm39+dm6 reads |
| `CHROM_SIZES_HUMAN` | hg38 chromosome sizes file |
| `CHROM_SIZES_MOUSE` | mm39 chromosome sizes file |
| `GTF_HUMAN` | hg38 GTF used to derive a TSS BED when none is supplied |
| `GTF_MOUSE` | mm39 GTF used to derive a TSS BED when none is supplied |
| `TSS_BED_HG38` | Optional strand-aware hg38 BED6 TSS file |
| `TSS_BED_MM39` | Optional strand-aware mm39 BED6 TSS file |
| `CCRE_BED_HG38` | Native hg38 ENCODE4 cCRE BED/bed.gz; required for an hg38 run when cCRE annotation is enabled |
| `CCRE_BED_MM39` | mm39 cCRE BED; required for an mm39 run when cCRE annotation is enabled |
| `CCRE_SOURCE_HG38` | Provenance label copied into hg38 annotation tables |
| `CCRE_SOURCE_MM39` | Provenance label copied into mm39 annotation tables |
| `CHROM_SIZES_DM6` | dm6 chromosome sizes used to validate the declared reference |
| `BLACKLIST_DM6` | dm6 blacklist BED applied before the external-reference denominator is counted |

See [Peak annotation](16_peak_annotation.md) for the exact GTF category rules,
cCRE class meanings, expected reference construction and output columns.

#### Blacklist regions

| Variable | Description |
|---|---|
| `BLACKLIST_HG38` | hg38 blacklist BED (ENCODE ENCFF356LFX recommended) |
| `BLACKLIST_MM39` | mm39 blacklist BED (Boyle lab recommended) |

#### Software

| Variable | Description |
|---|---|
| `PICARD_JAR` | Absolute path to `picard.jar` |
| `BEDGRAPH_TO_BIGWIG` | Absolute path to `bedGraphToBigWig` binary |
| `R_BIN` | Rscript binary (usually just `Rscript`) |

#### Thread counts

| Variable | Step | Recommended |
|---|---|---|
| `THREADS_PARALLEL_JOBS` | Max samples in parallel | 8 |
| `QC_SAMPLE_PARALLEL_JOBS` | Concurrent post-alignment metric, karyogram and consensus-FRiP samples | 4 |
| `ATAQV_PARALLEL_JOBS` | Concurrent ataqv/periodicity samples | 4 |
| `SPIKEIN_PARALLEL_JOBS` | Concurrent competitive alignment/calibration samples | 2 |
| `TRACK_PARALLEL_JOBS` | Concurrent DESeq2 track-generation samples | 2 |
| `POOLED_MACS_PARALLEL_JOBS` | Concurrent pooled MACS3 groups | 2 |
| `MERGE_PARALLEL_JOBS` | Concurrent replicate-merge/CPM groups | 2 |
| `THREADS_ALIGN` | Bowtie2 `-p` per sample | 16 |
| `THREADS_SAMTOOLS` | samtools `-@` | 16 |
| `THREADS_TRIMGALORE` | TrimGalore `--cores` | 8 |
| `THREADS_FASTQC` | FastQC `-t` | 10 |
| `THREADS_BIGWIG` | bedtools / samtools for coverage | 16 |
| `THREADS_DEEPTOOLS` | deepTools workers for Step 10 QC | 16 |
| `THREADS_ATAQV` | ataqv TSS calculation | 8 |

> **`THREADS_DEEPTOOLS`** is required as of v3.1.0. Add it to any existing config:
> `echo "THREADS_DEEPTOOLS=16" >> config.conf`

Peak concurrent CPU usage depends on the stage. For example, ataqv can use
`ATAQV_PARALLEL_JOBS × THREADS_ATAQV`, while DESeq2 track generation can use
`TRACK_PARALLEL_JOBS × THREADS_BIGWIG`. Large BAM scans are frequently limited
by storage throughput, so benchmark two and four jobs before increasing these
values further.

#### Cleanup

| Variable | Description | Default |
|---|---|---|
| `ENABLE_AUTOMATIC_CLEANUP` | Delete selected intermediates after full success; set false to retain everything | `true` |
| `KEEP_INTERMEDIATE_BAMS` | Keep pre-dedup BAMs in `bams/` | `false` |
| `KEEP_TRIMMED_FASTQ` | Keep trimmed FASTQs in `trimmedFastq/` | `false` |
| `KEEP_DEDUP_BAMS` | Keep pre-blacklist deduplicated BAMs | `false` |
| `KEEP_FILTERED_BAMS` | Keep quantitative analysis BAMs | `true` |
| `KEEP_NORMALIZATION_POLICY_BAMS` | Keep temporary permissive/intermediate BAMs after full success | `false` |
| `KEEP_SPIKEIN_BAMS` | Keep composite, deduplicated and host/dm6 split spike-in BAMs after full success | `false` |
| `KEEP_RAW_BEDGRAPH` | Keep raw coverage bedGraphs in `bedGraph/` | `false` |

#### Consensus, tracks and ATAC QC

| Variable | Description | Default |
|---|---|---|
| `CONSENSUS_MIN_SAMPLES` | Distinct biological samples required to support a consensus interval | `2` |
| `ALLOW_SINGLE_SAMPLE_CONSENSUS` | Permit a one-sample fallback | `false` |
| `GENERATE_DESEQ2_CONSENSUS_TRACKS` | Generate consensus-count DESeq2-scaled tracks | `true` |
| `GENERATE_CPM_TRACKS` | Generate current-policy CPM coverage | `true` |
| `GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS` | Generate robust-CPM tracks from pre-dedup primary alignments at `PERMISSIVE_MIN_MAPQ` | `true` |
| `GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS` | Generate robust-CPM tracks from Picard-deduplicated primary alignments at `INTERMEDIATE_MIN_MAPQ` | `true` |
| `GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS` | Generate robust-CPM tracks from the current deduplicated `MIN_MAPQ` final BAM | `true` |
| `GENERATE_COVERAGE_BIGWIGS` | Generate bigWig format for every enabled family | `true` |
| `GENERATE_COVERAGE_BEDGRAPHS` | Generate bedGraph format for every enabled family | `true` |
| `PERMISSIVE_MIN_MAPQ` | MAPQ threshold for duplicates-retained sensitivity BAMs | `0` |
| `INTERMEDIATE_MIN_MAPQ` | MAPQ threshold for deduplicated sensitivity BAMs | `0` |
| `GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS` | Run independent stringent dm6 calibration; requires samplesheet declarations | `false` |
| `GENERATE_DROSOPHILA_CONTROL_TRACKS` | With spike-in mode enabled, generate raw and CPM dm6 control bigWig/bedGraph files | `true` |
| `SPIKEIN_MIN_MAPQ` | Minimum composite-alignment MAPQ for both host and dm6 | `30` |
| `SPIKEIN_SCALE_TARGET` | Fixed retained-dm6 target in the scale formula | `1000000` |
| `SPIKEIN_CANONICAL_CONTIGS` | Comma-separated dm6 chromosomes included in counts | `2L,2R,3L,3R,4,X` |
| `SPIKEIN_MIN_FRAGMENTS_FAIL` | Hard minimum retained dm6 fragments/reads | `1000` |
| `SPIKEIN_MIN_FRAGMENTS_WARN` | Low-count warning threshold | `10000` |
| `SPIKEIN_WARN_LOW_FRACTION` | Low dm6-fraction warning threshold | `0.001` |
| `SPIKEIN_WARN_HIGH_FRACTION` | High dm6-fraction warning threshold | `0.20` |
| `COVERAGE_FILTER_PARALLEL_JOBS` | Concurrent policy-BAM filtering jobs | `2` |
| `SE_SIGNAL_MODE` | SE normalization unit; currently only one retained `read` is supported | `read` |
| `RUN_ATAQV_QC` | Run TSS enrichment and ATAC fragment QC | `true` |
| `GENERATE_ATAQV_VIEWER` | Build the local interactive ataqv viewer | `true` |
| `FRAGMENT_PLOT_MAX_BP` | Maximum fragment length shown in periodicity plots | `1000` |
| `UCSC_BIGDATA_URL_BASE` | Optional public base URL for custom-track definitions | empty |
| `DIFFERENTIAL_CONDITION_ORDER` | Optional comma-separated direction/order for universal condition pairs | first samplesheet appearance |
| `DIFFERENTIAL_MIN_ABS_LOG2FC` | Additional absolute log2-fold-change threshold for significant-site summaries | `0` |
| `RUN_SIMPLE_PEAK_ANNOTATION` | Add GTF gene context/nearest TSS and consensus annotations | `true` |
| `RUN_CCRE_ANNOTATION` | Require and add genome-matched cCRE classes; set false for GTF-only annotation | `true` |
| `PEAK_ANNOTATION_PROMOTER_UPSTREAM` | Bases upstream of TSS in promoter definition | `2000` |
| `PEAK_ANNOTATION_PROMOTER_DOWNSTREAM` | Bases downstream of TSS in promoter definition | `500` |

The three robust policies share the one consensus BED produced in Step 10, but
each filtering policy receives its own consensus count matrix, DESeq2 size
factors and robust cohort constant. Setting a family switch to `false` also
prevents construction of the BAM branch used only by that family. `MAPQ=0`
branches are exploratory sensitivity outputs; the `MIN_MAPQ=30` stringent/current
branch remains the principal high-confidence coverage policy.

#### Deprecated variables (v3.0.x → v3.1.0)

The following variables were used by `run_chipqc.R` (Step 10 in v3.0.x) and are no longer
read by the pipeline. They may be safely removed from your config:

| Deprecated variable | Reason |
|---|---|
| `CHIPQC_ANNOTATION_HG38` | ChIPQC replaced by deepTools |
| `CHIPQC_ANNOTATION_MM39` | ChIPQC replaced by deepTools |
| `CHIPQC_BLACKLIST_HG38_RDS` | ChIPQC replaced by deepTools |
| `CHIPQC_BLACKLIST_MM39_RDS` | ChIPQC replaced by deepTools |
| `THREADS_CHIPQC` | Replaced by `THREADS_DEEPTOOLS` |

---

[← Installation](03_installation.md) | [Next: Running →](05_running.md)
