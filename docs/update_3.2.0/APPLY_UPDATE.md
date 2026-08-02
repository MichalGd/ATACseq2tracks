# Applying the ATACseq2tracks 3.2.0 repository update

This directory is a complete, reviewable repository snapshot based on the current GitHub `main` branch plus the proposed 3.2.0 update. It contains no `.git` directory and has not modified either local source clone or the GitHub repository.

## Safe update procedure

1. Select one canonical Git clone. Do not maintain both local duplicate clones independently.
2. Create a review branch and confirm the worktree has no unrelated changes to files in the update manifest.
3. Compare the canonical clone with this prepared snapshot.
4. Copy or merge files on the review branch. Do not blindly overwrite an active, site-specific project configuration.
5. Remove the four obsolete scripts listed in `UPDATE_MANIFEST.md` only after confirming no external job calls them.
6. Create the Conda environment, run static checks and complete a representative acceptance run before merging or tagging.

Example comparison from the shared workspace:

```powershell
git -C fastq2tracks.2.1 status --short
git diff --no-index fastq2tracks.2.1 development\ATACseq2tracks_repository_update_3.2.0
```

On Linux, restore executable bits:

```bash
chmod +x atacseq2tracks.sh fastq2tracks.sh scripts/*.sh scripts/*.py tests/*.sh
```

## Required configuration review

Set or verify:

- `SAMPLESHEET`, `OUTPUT_DIR`;
- Bowtie2 index, chromosome sizes, GTF and blacklist paths for the single genome in the run;
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
RUN_ATAQV_QC=true
GENERATE_ATAQV_VIEWER=true
RUN_PEAK_ANNOTATION=false
RUN_MOTIF_ENRICHMENT=false
ENABLE_AUTOMATIC_CLEANUP=false
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
- every expected `_RPM.bw` and `_DESeq2Consensus.bw` is non-empty;
- `consensus_peaks.bed` applies the configured biological-sample support;
- `consensus_sizeFactors.tsv` contains one positive finite factor per non-control sample;
- each ATAC sample has a non-empty compressed `*.ataqv.json.gz` and TSS enrichment value;
- each PE sample has periodicity PNG/PDF plots and quantitative metrics;
- DiffBind produces the requested contrast results and diagnostic figures;
- `ucsc_tracks.txt` contains the intended relative or public locations;
- a deliberately failed child job causes a non-zero stage exit;
- resume checkpoints skip completed work and rerun when their signature changes.

Compare tracks, peak counts, consensus counts and differential results with the previous workflow on a known dataset before using 3.2.0 for biological conclusions.

## Rollback

Apply the update on a Git branch. Roll back with a normal Git revert or by abandoning the review branch; do not use a destructive reset when unrelated changes are present.
