# Known issues and deployment notes

This repository preserves the uploaded scripts exactly. The notes below identify issues to review before production use.

## 1. Hard-coded absolute paths

Many scripts call other scripts or tools through `/home/micgdu/...` paths. These must exist or be edited for a new server.

Check with:

```bash
grep -R "/home/micgdu" -n scripts
```

See `docs/INSTALLATION.md` for the replacement table.

## 2. `readme.2.1.sh` is not runnable shell code

Despite its `.sh` suffix, `readme.2.1.sh` is a change-log note. It fails `bash -n` by design and should not be executed.

## 3. `genomecoverage_batch.1.0.sh` uses an undocumented fourth argument

The usage block lists:

```bash
./genomecoverage_batch.1.0.sh <input_folder> <genome> [max_jobs]
```

But the script chooses the human or mouse lower-level coverage script based on `$4`. Use:

```bash
./scripts/genomecoverage_batch.1.0.sh /data/results/dedupBams hg38 8 human
./scripts/genomecoverage_batch.1.0.sh /data/results/dedupBams mm39 8 mouse
```

The master workflow already passes the fourth argument.

## 4. `fastq2tracks.2.1_blocked.sh` assumes existing folders and files

The blocked variant comments out early folder creation and raw QC/trimming steps. It is not a clean end-to-end workflow. Use it only for restarts where `trimmedFastq/`, `multiQC/`, and related folders already exist.

## 5. Output move order may misplace normalized bedGraph files

In `fastq2tracks.2.1.sh`, these lines occur in this order:

```bash
mv "$2"/dedupBams/*.bedGraph.gz "$2"/bedGraph
mv "$2"/dedupBams/*Snorm*.bedGraph.gz "$2"/NormBedGraph
```

Because `*_Snorm.bedGraph.gz` also matches `*.bedGraph.gz`, normalized bedGraph files may be moved to `bedGraph/` before the second command can move them to `NormBedGraph/`.

A future patched version should move `*Snorm*.bedGraph.gz` first, then move the remaining raw `*.bedGraph.gz` files.

## 6. Trim Galore naming support should be tested on your input style

`trimgalore_batch.1.0.sh` attempts to support both `*_1.fq.gz` and `*_R1_001.fastq.gz` R1 names. Its output verification is most consistent for the `*_1.fq.gz` pattern. Test alternative Illumina-style names on a tiny dataset before a full run.

## 7. Report generator dependencies

`generate_pipeline_report.1.0.sh` checks for several R packages, but PDF rendering paths also call `kableExtra::kable_styling`. Install `kableExtra` if you need PDF output. HTML output is the safer default.

## 8. Report text contains older script names

The newly added report generator still mentions older script names such as `DNAfastqBigWig_human_main_v5_31aug2025.sh` in the generated R Markdown narrative. The current master script in this repository is `fastq2tracks.2.1.sh`.

## 9. The master script deletes intermediate files

The master script removes BAM files in `bams/` and gzipped trimmed FASTQ files in `trimmedFastq/` near the end. Comment out the cleanup lines if you want to keep all intermediate files.

## 10. No GitHub license is selected

Add a `LICENSE` file before publishing publicly if reuse rights should be clear.
