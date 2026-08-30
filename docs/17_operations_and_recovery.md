# 17 — Operations and recovery

## Shared, activation-free command

The recommended shared installation exposes one stable command:

```bash
atacseq2tracks --config /absolute/path/project/config/config.conf
```

Users do not run `conda activate`. `bin/atacseq2tracks` constructs a controlled
`PATH` from the managed main and ataqv environments and clears inherited Python
and R library variables before executing the versioned workflow. The only
required command-line input is the configuration; `SAMPLESHEET` remains inside
that file. `CONDA_ENV_ACTIVATE` is retained only as an inert legacy key.

Administrators install a versioned release and atomically update the stable
links with:

```bash
sudo bash utilities/install_shared_release.sh
```

Defaults are `/opt/bioinformatics/workflows/ATACseq2tracks/releases/<version>`,
`/opt/bioinformatics/workflows/ATACseq2tracks/current`, and
`/usr/local/bin/atacseq2tracks`. Environment variables documented in the
installer can override those deployment locations.

## Safe configuration

The user configuration is parsed as data, not sourced as shell code. Only
uppercase keys present in the shipped template (plus the v4.3.0 operational
keys) are accepted. Duplicate/unknown keys, shell substitutions, pipes,
redirections and malformed assignments are rejected. Relative `SAMPLESHEET`
and `OUTPUT_DIR` paths resolve against the user configuration directory.

For compatibility, downstream v4.2.0 modules source a generated file containing
only shell-quoted literal assignments. The original file is never executed.
The resolved forms are retained under `metadata/`.

## Plan and preflight

```bash
atacseq2tracks --config /absolute/path/config.conf --plan
atacseq2tracks --config /absolute/path/config.conf --preflight-only
```

`--plan` validates the configuration and samplesheet, writes metadata, and
prints the ordered stage list without checking software or starting analyses.
`--preflight-only` additionally performs the complete software, package,
reference and permission checks and captures provenance, then exits.

## Checkpoint resume and named boundaries

A normal restart with the same command skips stages whose signature-matching
checkpoint exists. v4.3.0 signatures include the resolved configuration,
samplesheet and workflow implementation files. A code, samplesheet or config
change therefore invalidates old checkpoints instead of silently mixing runs.
Each `.done` checkpoint has a JSON sidecar with stage name, signature, UTC time
and elapsed seconds.

```bash
atacseq2tracks --config /absolute/path/config.conf --from-stage peaks
atacseq2tracks --config /absolute/path/config.conf --stop-after qc
```

`--from-stage` forces that named stage and all following stages to rerun; prior
stages still require matching checkpoints. `--stop-after` exits successfully
after the named stage writes its checkpoint. Available names are shown by
`atacseq2tracks --help`.

Automatic cleanup can make an earlier stage impossible to reconstruct. For
example, rerunning filtering-sensitivity tracks needs pre-dedup/dedup BAMs. In
that case start from FASTQs in a new output directory or retain the required
intermediates in advance.

## Reports without rerunning analysis

```bash
/opt/bioinformatics/workflows/ATACseq2tracks/current/utilities/regenerate_reports.sh \
  --config /absolute/path/config.conf
```

This calls only the established report generator. It does not align reads,
call peaks, recreate tracks or refit differential models.

## Operational metadata

`<OUTPUT_DIR>/metadata/` contains:

- `resolved_config.conf` and `resolved_config.tsv` — literal effective settings;
- `user_config_path.txt` — submitted configuration path;
- `validated_sequencing_units.tsv` — one row per samplesheet row;
- `biological_libraries.tsv` — one row per `sample_id + replicate`;
- `technical_merge_audit.tsv` — ordered FASTQs merged before trimming;
- `planned_stages.tsv` — named execution order;
- `resource_budget.tsv` — maximum configured thread demand by major stage;
- `software_versions.tsv` and `reference_manifest.tsv` — preflight provenance;
- `workflow_events.tsv` — run, start, skip, completion, failure and timing events;
- `cleanup_manifest.tsv` — exact paths and sizes removed after full success.

`TOTAL_CPU_BUDGET=0` disables a CPU ceiling. Set a positive value and choose
`RESOURCE_CHECK_MODE="warn"`, `"error"` or `"off"`. This is a static upper-bound
check; it is not a scheduler or a memory estimator.

## Limitations

- The Bash workflow has stage-level checkpoints, not per-contrast scheduling.
- Named recovery cannot recreate intermediates already removed by cleanup.
- Provenance records configured references and tool versions but does not hash
  multi-gigabyte FASTQs, BAMs or Bowtie2 indexes.
- Scientific batch correction remains the responsibility of the statistical
  design; operational reproducibility does not remove biological confounding.
