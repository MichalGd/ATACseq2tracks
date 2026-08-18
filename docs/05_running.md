# 05 — Running the pipeline

[← Input files](04_inputs.md) | [Next: Pipeline steps →](06_pipeline_steps.md)

---

## Before you run

1. Conda environment activated: `conda activate ATACseq2tracks`
2. Samplesheet validated: `python3 scripts/validate_samplesheet.py samplesheet.csv`
3. Config file edited with correct paths
4. Working directory = **parent folder of `ATACseq2tracks/`**, not inside it

---

## Full pipeline run

```bash
cd /path/to/parent_of_ATACseq2tracks

nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    > /path/to/my_project/run.log 2>&1 &

# Save PID for later monitoring or kill
echo $! > /path/to/my_project/run.pid

# Follow the log
tail -f /path/to/my_project/run.log
```

---

## Monitoring progress

```bash
# Live log
tail -f /path/to/my_project/run.log

# Check which steps have completed
ls /path/to/my_project/analysis/.checkpoints/

# Check if the process is still running
ps -p $(cat /path/to/my_project/run.pid)
```

Each completed step prints:
```
[CHECKPOINT] Step N complete -- to re-run: rm .../stepN.done
```

---

## Stopping the pipeline

```bash
# Graceful stop (waits for current background jobs)
kill $(cat /path/to/my_project/run.pid)

# Force stop
kill -9 $(cat /path/to/my_project/run.pid)
```

Any step that had not yet written its `.done` file will be repeated on resume.

---

## Resuming an interrupted run

```bash
# Just re-run the same command — completed steps are skipped automatically
nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    >> /path/to/my_project/run_resume.log 2>&1 &
```

---

## Rerunning individual steps

Delete the checkpoint file for that step and rerun:

```bash
# Rerun step 9 (MACS2) only
rm /path/to/my_project/analysis/.checkpoints/step9.done

nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    >> /path/to/my_project/run_rerun.log 2>&1 &
```

To rerun multiple steps:
```bash
rm /path/to/my_project/analysis/.checkpoints/step{6,7,8,9,10,11}.done
```

DiffBind and DESeq2ATAC resume independently:

```bash
rm /path/to/my_project/analysis/.checkpoints/step12.done   # DiffBind only
rm /path/to/my_project/analysis/.checkpoints/step12a.done  # DESeq2ATAC only
```

Each differential module refits its broad/narrow model and rewrites that model's
condition-pair outputs when rerun; status files are diagnostic, not independent
contrast checkpoints. Broad/narrow and peer modules are failure-isolated. If a
differential module fails, the workflow still writes the consolidated TSV/HTML
summary, exits non-zero, and does not perform automatic cleanup.

To rerun everything from scratch:
```bash
rm -rf /path/to/my_project/analysis/.checkpoints/
```

---

## Running a second project

Each project needs its own config and samplesheet. The `ATACseq2tracks/` installation is shared.

```bash
mkdir -p /path/to/project_B/config
cp config/config.conf /path/to/project_B/config/config.conf
# Edit config_B with different SAMPLESHEET and OUTPUT_DIR

nohup bash ATACseq2tracks/atacseq2tracks.sh \
    --config /path/to/project_B/config/config.conf \
    > /path/to/project_B/run.log 2>&1 &
```

---

## Running paired-end and single-end projects

The same workflow installation handles both layouts, but they are separate analyses. Use one samplesheet, configuration and output directory for PE, and another set for SE:

```bash
# Paired-end project
bash ATACseq2tracks/atacseq2tracks.sh \
  --config /path/to/project_PE/config/config.conf

# Single-end project; config contains SE_SIGNAL_MODE="read"
bash ATACseq2tracks/atacseq2tracks.sh \
  --config /path/to/project_SE/config/config.conf
```

For SE rows, set `layout=SE` and leave `fastq_2` empty. Each filtered SE alignment is counted once for CPM, consensus-peak counts and DESeq2 factors. Do not merge PE and SE outputs into one DESeq2 normalization cohort.

---

## Running only certain steps manually

Each script in `scripts/` can be called independently. This is useful for testing or reprocessing a single sample.

```bash
# Example: run blacklist filtering for one sample
bash ATACseq2tracks/scripts/blacklist_filter.sh \
    /path/to/analysis/dedupBams/NHEK_H3K27ac_day0_bioR1_dedup.bam \
    /path/to/blacklist_hg38.bed \
    /path/to/analysis/filteredBams/

# Example: call peaks for one sample
bash ATACseq2tracks/scripts/macs2_peaks.sh \
    /path/to/filteredBams/NHEK_H3K27ac_day0_bioR1_dedup_blFilt.bam \
    /path/to/filteredBams/NHEK_Input_bioR1_dedup_blFilt.bam \
    /path/to/peaks/NHEK_H3K27ac_day0_bioR1 \
    both hg38 NHEK_H3K27ac_day0_bioR1
```

---

[← Input files](04_inputs.md) | [Next: Pipeline steps →](06_pipeline_steps.md)
