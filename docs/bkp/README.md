# fastq2tracks — Documentation Index

Welcome to the **fastq2tracks** documentation.

fastq2tracks is a checkpoint-based, samplesheet-driven pipeline that processes raw FASTQ files from ChIP-seq, ATAC-seq, CUT&RUN, and CUT&TAG assays all the way to bigwig tracks, MACS2 peaks (narrow and broad), ChIPQC quality reports, and DiffBind-ready samplesheets.

---

## Documentation pages

| Page | Contents |
|---|---|
| [01 — Overview](01_overview.md) | What the pipeline does, design principles, comparison with rnaseq2tracksP |
| [02 — Quick start](02_quickstart.md) | Fastest path from installation to first run |
| [03 — Installation](03_installation.md) | Conda environment, R packages, reference files |
| [04 — Input files](04_inputs.md) | Samplesheet format, column reference, config file parameters |
| [05 — Running the pipeline](05_running.md) | Launch commands, monitoring, resume, partial reruns |
| [06 — Pipeline steps](06_pipeline_steps.md) | Every step explained with inputs, outputs, and key logic |
| [07 — Outputs](07_outputs.md) | Output directory tree, file naming conventions, how to read results |
| [08 — Downstream: DiffBind](08_diffbind.md) | Using pipeline outputs for differential binding analysis |
| [09 — Troubleshooting](09_troubleshooting.md) | Common errors and fixes, all v3.0.4 bugs documented |
| [10 — Reference file preparation](10_reference_files.md) | Building ChIPQC RDS annotations and blacklist objects |

---

## Quick navigation

- **New user?** Start with [Quick start](02_quickstart.md)
- **Setting up a server?** See [Installation](03_installation.md) and [Reference file preparation](10_reference_files.md)
- **Preparing your samplesheet?** See [Input files](04_inputs.md)
- **Something failed?** See [Troubleshooting](09_troubleshooting.md)
- **Pipeline ran — what do I have?** See [Outputs](07_outputs.md)

---

Back to [main README](../README.md)
