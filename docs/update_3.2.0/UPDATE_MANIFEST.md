# ATACseq2tracks 3.2.0 update manifest

This folder is a complete repository snapshot prepared from GitHub `main` as inspected on 2026-08-02, plus the proposed 3.2.0 changes. It is not a Git clone and intentionally contains no `.git` directory.

## Source baseline

- Repository: `https://github.com/MichalGd/ATACseq2tracks`
- Branch: `main`
- Baseline commit: `4a2eabb5caa43a70ac08cc0ad449d26aba002bed`
- Remote tree at inspection: 59 tracked files
- Both local baseline copies matched all 59 remote Git blob hashes.

## Replaced runtime files

- `atacseq2tracks.sh`, `fastq2tracks.sh`
- `environment.yml`
- `config/config.conf`, `config/config_temp.conf.template`
- `scripts/blacklist_filter.sh`
- `scripts/blacklist_filter_batch.sh`
- `scripts/bowtie2_batch.sh`
- `scripts/consensus_peak_size_factors.R`
- `scripts/create_ucsc_tracks.sh`
- `scripts/fastqc_batch.sh`
- `scripts/genomecoverage_batch.sh`
- `scripts/genomecoverage_single.sh`
- `scripts/macs2_batch.sh`
- `scripts/macs2_peaks.sh`
- `scripts/merge_replicates.sh`
- `scripts/picard_dedup_batch.sh`
- `scripts/post_alignment_qc_batch.sh`
- `scripts/prepare_diffbind.R`
- `scripts/diffbind_analysis.sh`
- `scripts/diffbind_analysis.R`
- `scripts/smoke_test.sh`
- `scripts/trimgalore_batch.sh`
- `scripts/validate_samplesheet.py`

## New runtime files

- `VERSION`
- `config/samplesheet_example_atac.csv`
- `scripts/ataqv_qc_batch.sh`
- `scripts/extract_ataqv_metrics.py`
- `scripts/peak_interpretation.sh`
- `scripts/plot_fragment_periodicity.py`
- `scripts/prepare_tss_bed.py`

## Retained files needed by the entry point

- `scripts/plot_chrom_coverage.py`
- `scripts/generate_pipeline_report.sh`
- all current documentation, examples and syntax tests, subject to the 3.2.0 documentation notes.

## Obsolete files to remove from the target repository

These are not used by the 3.2.0 entry point and are excluded from this prepared snapshot:

- `scripts/run_chipqc.R`
- `scripts/validate_samplesheet2.py`
- `scripts/bowtie2_align.sh`
- `scripts/picard_dedup.sh`

Deletion should be performed in a review branch after confirming no external automation calls these legacy helpers.

## Files intentionally not required for runtime

The historical PowerPoint presentation is retained only for repository history. It is not needed to execute the workflow.

## Application policy

Do not overwrite an active project-specific `config.conf`. Merge the new variables into it, verify every reference/software path, create the Conda environment, run static checks, and complete the acceptance test in `APPLY_UPDATE.md` before merging or tagging the release.
