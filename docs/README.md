# ATACseq2tracks documentation

Current prepared version: **4.2.0**.

| Page | Contents |
|---|---|
| [01 — Overview](01_overview.md) | Scope, outputs, workflow sequence and limitations |
| [02 — Quick start](02_quickstart.md) | Fastest path from installation to first run |
| [03 — Installation](03_installation.md) | Conda, R and reference requirements |
| [04 — Inputs](04_inputs.md) | Samplesheet and configuration variables |
| [05 — Running](05_running.md) | Launch, monitoring and checkpoint resume |
| [06 — Pipeline stages](06_pipeline_steps.md) | Stage-by-stage implementation |
| [07 — Outputs](07_outputs.md) | Output layout and file interpretation |
| [08 — DiffBind](08_diffbind.md) | DiffBind usage |
| [09 — Troubleshooting](09_troubleshooting.md) | Common failures and recovery |
| [10 — References](10_reference_files.md) | Bowtie2 indices, GTF, chromosome sizes and blacklists |
| [11 — Blacklist filtering](11_blacklist_filtering.md) | Filtering policy and sources |
| [12 — Post-alignment QC](12_post_alignment_qc.md) | deepTools, TSS enrichment and periodicity QC |
| [13 — Differential accessibility](13_differential_accessibility.md) | Statistical design and results |
| [14 — Replicates and design](14_replicates_and_experimental_design.md) | Biological/technical replication guidance |
| [DESeq2ATAC legacy-method review](15_deseq2atac_legacy_method_review.md) | Legacy script audit, retained behavior, modernizations and limitations |
| [16 — Peak annotation](16_peak_annotation.md) | GTF categories, cCRE classes, output columns, provenance and limitations |
| [v4.0.0 release](v4.0.0_RELEASE_2026-08-19.md) | Current release scope, configurable parallelism and validation |
| [v4.0.0 reporting/DiffBind hotfix](v4.0.0_REPORTING_DIFFBIND_HOTFIX_2026-08-20.md) | Canonical pre-count peak filtering and stable MultiQC deepTools presentation |
| [v4.1.0 filtering-policy coverage update](v4.1.0_READ_FILTERING_COVERAGE_UPDATE_2026-08-20.md) | Five default coverage families, independent switches and policy-specific normalization |
| [v4.2.0 Drosophila spike-in update](v4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md) | Competitive host+dm6 alignment, calibration formula, QC, outputs and limitations |
| [v3.2.0 track-normalization update](v3.2.0_TRACK_NORMALIZATION_UPDATE_2026-08-06.md) | PE fragment and SE read semantics, output contract, validation and resume instructions |
| [v3.2.0 DiffBind/DESeq2ATAC update](v3.2.0_DESEQ2ATAC_UPDATE_2026-08-13.md) | Configuration, deployment, checkpoints, validation and affected files |
| [v3.2.0 universal differential/annotation update](v3.2.0_UNIVERSAL_DIFFERENTIAL_ANNOTATION_UPDATE_2026-08-18.md) | Any-number condition pairs, singleton policy, pair summaries and GTF/cCRE annotation |
| [Original 3.2.0 manifest](update_3.2.0/UPDATE_MANIFEST.md) | Historical 2026-08-02 baseline plus links to cumulative addenda |
| [3.2.0 application guide](update_3.2.0/APPLY_UPDATE.md) | Safe review, update and acceptance testing for the cumulative snapshot |
| [Original 3.2.0 critical scope](update_3.2.0/MINIMAL_CRITICAL_IMPROVEMENTS.md) | Historical initial scope; superseded output details are labeled |

For first use, read the [main README](../README.md), then follow the [v4.2.0 application guide](../APPLY_V4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md).
