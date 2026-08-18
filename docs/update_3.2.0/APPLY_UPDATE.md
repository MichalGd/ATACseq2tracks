# Applying the ATACseq2tracks 3.2.0 repository update

The prepared directory
`development/ATACseq2tracks_v3.2.0_universal_multicondition_annotation_update_2026-08-18`
is the complete cumulative review snapshot. It has no `.git` directory and does
not modify a local clone or the public repository by itself. The public `main`
branch may still represent an earlier 3.2.0 state, so compare the complete trees
rather than assuming matching version labels imply matching contents.

## Safe update procedure

1. Select one canonical Git clone. Do not maintain both local duplicate clones independently.
2. Create a review branch and confirm the worktree has no unrelated changes.
3. Compare the canonical clone with this prepared snapshot.
4. Copy or merge files on the review branch. Do not blindly overwrite an active, site-specific project configuration.
5. Review additions, replacements and deletions from the current tree comparison;
   the 2026-08-02 manifest is a historical baseline, not a complete list of the
   later cumulative updates.
6. Create the Conda environment, run static checks and complete a representative acceptance run before merging or tagging.

Example comparison from the shared workspace:

```powershell
git -C D:\bioinformatics\github\ATACseq2tracks status --short
git diff --no-index D:\bioinformatics\github\ATACseq2tracks D:\bioinformatics\ATAC_VCS\development\ATACseq2tracks_v3.2.0_universal_multicondition_annotation_update_2026-08-18
```

On Linux, restore executable bits:

```bash
chmod +x atacseq2tracks.sh fastq2tracks.sh scripts/*.sh scripts/*.py tests/*.sh
```

## Required configuration review

Set or verify:

- `SAMPLESHEET`, `OUTPUT_DIR`;
- Bowtie2 index, chromosome sizes, GTF and blacklist paths for the single genome in the run;
- the matching cCRE BED and source label, unless
  `RUN_CCRE_ANNOTATION=false` is deliberately selected for GTF-only annotation;
- optionally `TSS_BED_HG38` or `TSS_BED_MM39`; otherwise TSS BED is generated from GTF;
- `PICARD_JAR`, `PICARD_TMP` and `BEDGRAPH_TO_BIGWIG`;
- `UCSC_BIGDATA_URL_BASE` when tracks must be loaded directly by UCSC;
- thread counts appropriate for the server.

Review these defaults deliberately:

```bash
MIN_MAPQ=30
REMOVE_MITO=true
CONSENSUS_PEAK_TYPE="narrow"
CONSENSUS_MIN_SAMPLES=2
ALLOW_SINGLE_SAMPLE_CONSENSUS=false
GENERATE_DESEQ2_CONSENSUS_TRACKS=true
DIFFBIND_SUMMITS=100
RUN_DESEQ2ATAC=true
DESEQ2ATAC_MIN_SAMPLES=2
DESEQ2ATAC_ALPHA=0.05
DESEQ2ATAC_BLOCK_COLUMN=""
DESEQ2ATAC_REFERENCE_CONDITION=""
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=true
RUN_ATAQV_QC=true
GENERATE_ATAQV_VIEWER=true
RUN_PEAK_ANNOTATION=false
RUN_MOTIF_ENRICHMENT=false
ENABLE_AUTOMATIC_CLEANUP=true
KEEP_FILTERED_BAMS=true
```

## Static validation

```bash
bash tests/check_bash_syntax.sh
python3 scripts/validate_samplesheet.py config/samplesheet_example_atac.csv
bash scripts/smoke_test.sh /path/to/samplesheet.csv /path/to/config.conf
```

ShellCheck is also recommended for all Bash entry points and modules.

## Acceptance run

Use a small paired-end ATAC-seq dataset with at least two biological samples and run Steps 0–14. Verify:

- all expected filtered BAMs pass `samtools quickcheck`;
- filtering attrition tables reconcile with alignment and deduplication counts;
- PE MACS3 logs show `BAMPE` and requested peak files are non-empty;
- every expected `_CPM.bw`, `_DESeq2Consensus.{bw,bedGraph}` and
  `_DESeq2RobustCPM.{bw,bedGraph}` is non-empty, and no CPM bedGraph is produced;
- `consensus_peaks.bed` applies the configured biological-sample support;
- `consensus_sizeFactors.tsv` contains one positive finite factor per non-control sample;
- each ATAC sample has a non-empty compressed `*.ataqv.json.gz` and TSS enrichment value;
- each PE sample has periodicity PNG/PDF plots and quantitative metrics;
- DiffBind produces the requested contrast results and diagnostic figures;
- DiffBind summaries report the configured 100-bp summit half-width;
- `deseq2atac/broad/` and `deseq2atac/narrow/` contain independently supported
  consensuses, readable compressed count and result tables, explicit contrast
  direction, and non-empty PNG/PDF figures;
- every complete DiffBind and DESeq2ATAC consensus has an annotation table, and
  annotation columns propagate to complete and significant pair tables;
- a multi-condition fixture exports all `k*(k-1)/2` pairs among eligible
  conditions, while singleton-condition samples remain in consensus construction
  but not in statistical models;
- PE and SE fixtures confirm one fragment and one read per observation,
  respectively; do not combine both layouts in one run;
- a no-hit DESeq2ATAC fixture completes successfully with a header-only
  significant table and an explicit zero-significant summary;
- each populated bigWig directory has `ucsc_tracks.txt` with the intended
  relative or public locations;
- a deliberately failed child job causes a non-zero stage exit;
- resume checkpoints skip completed work and rerun when their signature changes.

For an existing output directory, remove only `.checkpoints/step12.done` to
rerun DiffBind and `.checkpoints/step12a.done` to rerun DESeq2ATAC. Step 12a is
new and does not reuse the established Step 13/14 checkpoint names.

The new settings have code defaults. When resuming an existing checkpointed run
with those defaults, do not edit its project config merely to add the variables:
the global run signature includes the entire config and a config edit makes all
older checkpoints stale. Add settings to an existing config only when overriding
a default is necessary.

Compare tracks, peak counts, consensus counts and differential results with the previous workflow on a known dataset before using 3.2.0 for biological conclusions.

## Rollback

Apply the update on a Git branch. Roll back with a normal Git revert or by abandoning the review branch; do not use a destructive reset when unrelated changes are present.
