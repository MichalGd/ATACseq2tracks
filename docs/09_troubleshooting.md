# 09 — Troubleshooting

[← DiffBind](08_diffbind.md) | [Next: Reference files →](10_reference_files.md)

---

## How to diagnose a failure

1. Check the main log for the step that failed:
   ```bash
   grep -A5 "Error\|FAIL\|error" /path/to/analysis/run.log
   ```
2. Check the per-sample log in the step's log directory:
   ```bash
   ls /path/to/analysis/filteredBams/blacklist_logs_*/
   cat /path/to/analysis/filteredBams/blacklist_logs_*/<sample_id>.log
   ```
3. Check which checkpoints exist:
   ```bash
   ls /path/to/analysis/.checkpoints/
   ```

---

## Known issues and fixes

### Step 2 — TrimGalore fails with "unknown option"

**Symptom:**
```
trim_galore: error: unrecognized arguments: -o
```
**Cause:** Old version of Trim Galore (< 2.2.0). The pipeline uses `--output_dir`.
**Fix:** Update to Trim Galore v2.2.0+ (Oxidized Edition) via conda or from source.

---

### Step 4 — `THREADS_BOWTIE2: unbound variable`

**Symptom:**
```
/scripts/bowtie2_batch.sh: line XX: THREADS_BOWTIE2: unbound variable
```
**Cause:** Old config file using `THREADS_BOWTIE2`. The pipeline uses `THREADS_ALIGN`.
**Fix:** Rename the variable in your `config.conf`:
```bash
# Old:
THREADS_BOWTIE2=16
# New:
THREADS_ALIGN=16
```

---

### Step 5 — `samtools sort: failed to read header`

**Cause:** Downstream of the `THREADS_BOWTIE2` / `THREADS_ALIGN` bug — Bowtie2 exited early, producing an empty or truncated SAM.
**Fix:** Same as above — fix the variable name in config, delete checkpoints 4 and 5, rerun.

---

### Step 5 — Picard finds no BAMs

**Symptom:**
```
[timestamp] Found 0 BAMs to process
```
**Cause:** Picard batch script was scanning for `*.sorted_stChr.bam` (v2.1 naming). Fixed in v3.0.4 to scan all `*.bam`.
**Fix:** Ensure you are running the v3.0.4 script.

---

### Steps 6–13 — All samples skipped ("BAM not found")

**Symptom:**
```
SKIP <sample_id> — dedup BAM not found: .../dedupBams/<sample_id>_dedup.bam
```
**Cause:** Picard dedup output includes `_bioR<replicate>` suffix (e.g. `<sample_id>_bioR1_dedup.bam`) but the lookup in downstream scripts used only `<sample_id>_dedup.bam`.
**Fix:** Fixed in v3.0.4. All scripts from step 6 onward now construct the BAM key as `${sample_id}_bioR${replicate}_dedup.bam`.

---

### Step 6 — `unexpected EOF while looking for matching '"'`

**Symptom:**
```
/scripts/blacklist_filter_batch.sh: line 45: unexpected EOF while looking for matching `"'
```
**Cause:** Missing closing `}` braces in `_load_config()` and `wait_for_slot()` functions.
**Fix:** Replace the script with the v3.0.4 version (closing braces are present).

---

### Step 9 — `ctrl_rep: unbound variable`

**Symptom:**
```
/scripts/macs2_batch.sh: line XX: ctrl_rep: unbound variable
```
**Cause:** `ctrl_rep` was used but never assigned. Fixed in v3.0.4 by adding a pre-built `sid_rep[]` associative array lookup.
**Fix:** Use the v3.0.4 `macs2_batch.sh`.

---

### Step 10 — `THREADS_DEEPTOOLS: unbound variable`

**Symptom:**
```
/scripts/post_alignment_qc_batch.sh: line XX: THREADS_DEEPTOOLS: unbound variable
```
**Cause:** `THREADS_DEEPTOOLS` is missing from your `config.conf`.
**Fix:** Add the following line to `config.conf`:
```bash
THREADS_DEEPTOOLS=8
```

---

### Step 10 — `bamCoverage: error: the following arguments are required`

**Symptom:**
```
bamCoverage: error: the following arguments are required: --bam/-b
```
**Cause:** deepTools version < 3.5 uses a different argument syntax for `bamCoverage`.
**Fix:** Update deepTools to >= 3.5:
```bash
conda activate ATACseq2tracks
mamba install deeptools>=3.5
```

---

### Step 10 — Karyogram plot fails ("No module named matplotlib")

**Symptom:**
```
ModuleNotFoundError: No module named 'matplotlib'
```
**Cause:** `matplotlib`, `numpy`, or `pandas` are not installed in the active conda environment.
**Fix:**
```bash
conda activate ATACseq2tracks
mamba install matplotlib numpy pandas
```

---

### Step 10 — All samples show `NO_PEAKS` in QC summary

**Symptom:** `qc_post_alignment/tables/qc_summary.tsv` shows `NO_PEAKS` for every sample; fingerprint and bins-level QC still run.
**Cause:** Peak calling (Step 9) failed or produced empty peak files, or peak files are in an unexpected location.
**Fix:**
1. Check `peaks/per_replicate/` to confirm `.narrowPeak` and `.broadPeak` files exist and are non-empty.
2. If peak calling failed, delete `step9.done` and rerun.
3. If peaks exist but paths are wrong, check the `KEY` construction logic in `post_alignment_qc_batch.sh`.

---

### Samplesheet warnings about `macs2_mode`

**Symptom:**
```
[WARN]  row 2: Control sample: macs2_mode should be 'none'
```
**Cause:** Control rows in the samplesheet have `macs2_mode` set to something other than `none`.
**Action:** These are warnings, not errors. The pipeline will still run. Update the samplesheet to set `macs2_mode=none` for all `is_control=TRUE` rows to suppress the warnings.

---

### MultiQC reports no results

**Symptom:**
```
multiqc | No analysis results found.
```
**Cause:** The input directory contained no FastQC or alignment log files — usually because a previous step failed silently (BAMs not found, trimming output misnamed).
**Fix:** Check the upstream step logs to ensure files were actually produced.

---

### R: package not found (DiffBind)

**Symptom:**
```
Error in library(DiffBind) : there is no package called 'DiffBind'
```
**Fix:** Install missing packages — see [Installation](03_installation.md#r-packages).

---

### Disk full mid-run

**Symptom:** Pipeline fails partway through with write errors.
**Fix:**
1. Free disk space
2. Delete any incomplete output files from the failed step
3. Remove the checkpoint for that step
4. Rerun

Pre-flight check (step 0) requires 50 GB free. For 30+ samples with PE reads, budget at least 500 GB.

---

## General checklist

- [ ] Running from the **parent** of `ATACseq2tracks/`, not from inside it
- [ ] `F2T_CONFIG` is exported OR `--config` is passed to the master script
- [ ] All FASTQ paths in the samplesheet are **absolute**, not relative
- [ ] `control_id` values match `sample_id` values **exactly** (case-sensitive)
- [ ] `macs2_mode=none` for all `is_control=TRUE` rows
- [ ] `tech_replicate` is `1` for samples with only one sequencing run
- [ ] Preflight reports both inputs as UTF-8/LF or creates a timestamped
  `*.windows-artifact-backup.*` file and reports that it normalized them.
- [ ] `THREADS_DEEPTOOLS` is set in `config.conf`
- [ ] deepTools >= 3.5 is installed in the conda environment

---

[← DiffBind](08_diffbind.md) | [Next: Reference files →](10_reference_files.md)
