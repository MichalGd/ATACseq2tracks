# Contributing

Contributions to ATACseq2tracks are welcome. This project is primarily used for production
ChIP-seq / ATAC-seq / CUT&RUN analysis in our laboratory, so stability and reproducibility
are the top priorities.

---

## Guidelines

1. **Keep original scripts unchanged** unless a version bump is intentional.
   Document any code change in `CHANGELOG.md` with a version bump.

2. **Run syntax checks before committing:**
   ```bash
   bash -n scripts/*.sh
   python3 -m py_compile scripts/*.py
   ```

3. **Test on a small dataset** before running on production data.
   A paired-end dataset with 2 samples (1 IP + 1 input, ~1 M reads each) is sufficient
   to validate the full pipeline.

4. **Do not commit generated data files:**
   FASTQ, BAM, bedGraph, bigWig, peak files, and run-output folders must not be committed.
   The `.gitignore` already excludes the most common patterns.

5. **Update documentation** when adding or modifying pipeline steps, scripts, or config
   parameters. The relevant files are in `docs/`.

6. **Config variables:** any new config parameter must be:
   - Added to `config/config.conf` with a comment explaining its purpose
   - Validated in `scripts/smoke_test.sh`
   - Documented in `docs/04_inputs.md` (configuration file section)

---

## Reporting issues

Please open a GitHub issue with:
- The pipeline version (see badge in README or `CHANGELOG.md`)
- The step that failed (step number + script name)
- The relevant lines from the run log
- The output of `conda list | grep -E 'bowtie2|samtools|deeptools|macs|trim'`

---

## Adding a new assay type

New assay types (e.g. CUT&TAG-seq, ChIRP-seq) should be validated against the samplesheet
validator (`scripts/validate_samplesheet.py`) and tested end-to-end before documentation
is updated in `docs/01_overview.md`.
