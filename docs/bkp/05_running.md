# 05 — Running the pipeline

[← Input files](04_inputs.md) | [Next: Pipeline steps →](06_pipeline_steps.md)

---

## Before you run

1. Conda environment activated: `conda activate fastq2tracks`
2. Samplesheet validated: `python3 scripts/validate_samplesheet.py samplesheet.csv`
3. Config file edited with correct paths
4. Working directory = **parent folder of `fastq2tracks/`**, not inside it

---

## Full pipeline run

```bash
cd /path/to/parent_of_fastq2tracks

nohup bash fastq2tracks/fastq2tracks.sh \
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
nohup bash fastq2tracks/fastq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    >> /path/to/my_project/run_resume.log 2>&1 &
```

---

## Rerunning individual steps

Delete the checkpoint file for that step and rerun:

```bash
# Rerun step 9 (MACS2) only
rm /path/to/my_project/analysis/.checkpoints/step9.done

nohup bash fastq2tracks/fastq2tracks.sh \
    --config /path/to/my_project/config/config.conf \
    >> /path/to/my_project/run_rerun.log 2>&1 &
```

To rerun multiple steps:
```bash
rm /path/to/my_project/analysis/.checkpoints/step{6,7,8,9,10,11}.done
```

To rerun everything from scratch:
```bash
rm -rf /path/to/my_project/analysis/.checkpoints/
```

---

## Running a second project

Each project needs its own config and samplesheet. The `fastq2tracks/` installation is shared.

```bash
mkdir -p /path/to/project_B/config
cp config/config.conf /path/to/project_B/config/config.conf
# Edit config_B with different SAMPLESHEET and OUTPUT_DIR

nohup bash fastq2tracks/fastq2tracks.sh \
    --config /path/to/project_B/config/config.conf \
    > /path/to/project_B/run.log 2>&1 &
```

---

## Running only certain steps manually

Each script in `scripts/` can be called independently. This is useful for testing or reprocessing a single sample.

```bash
# Example: run blacklist filtering for one sample
bash fastq2tracks/scripts/blacklist_filter.sh \
    /path/to/analysis/dedupBams/NHEK_H3K27ac_day0_bioR1_dedup.bam \
    /path/to/blacklist_hg38.bed \
    /path/to/analysis/filteredBams/

# Example: call peaks for one sample
bash fastq2tracks/scripts/macs2_peaks.sh \
    /path/to/filteredBams/NHEK_H3K27ac_day0_bioR1_dedup_blFilt.bam \
    /path/to/filteredBams/NHEK_Input_bioR1_dedup_blFilt.bam \
    /path/to/peaks/NHEK_H3K27ac_day0_bioR1 \
    both hg38 NHEK_H3K27ac_day0_bioR1
```

---

[← Input files](04_inputs.md) | [Next: Pipeline steps →](06_pipeline_steps.md)
