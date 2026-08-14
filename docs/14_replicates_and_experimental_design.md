# 14 — Replicates and experimental design

[← Differential accessibility](13_differential_accessibility.md) | [Next: Outputs →](07_outputs.md)

---

This page explains how `ATACseq2tracks` handles biological and technical replicates, and how to structure your experiment for reliable ATAC-seq analysis.

## Biological vs technical replicates

- **Biological replicates** are independent biological samples.
  - Examples: different animals, different donors, separate cell culture flasks.
  - They capture true biological variation and are required for statistical inference.

- **Technical replicates** are repeated measurements of the same biological sample.
  - Examples: the same library sequenced on two lanes, the same DNA prep split across two runs.
  - They capture technical variability and should be merged before modeling.

## Pipeline handling

`ATACseq2tracks` uses the samplesheet fields `sample_id`, `replicate`, and `tech_replicate` to:

- merge technical replicates early during trimming and FASTQ handling
- keep biological replicates separate for peak calling and quantification
- generate per-biological-replicate BAMs and tracks
- generate pooled/merged tracks for condition-level visualization only

### Technical replicate handling

Technical replicates must have the same:

- `sample_id`
- `replicate`
- `layout`,`genome`,`assay`,`factor`,`condition`,`treatment`,`cell_type`,`is_control`,`control_id`,`macs2_mode`,`blacklist`,`output_prefix`

They may differ only in `fastq_1`, `fastq_2`, and `tech_replicate`.

If your samplesheet contains multiple rows with identical `sample_id` and `replicate`, the pipeline treats them as technical replicates and merges the FASTQs before trimming.

### Biological replicate handling

Biological replicates should remain separate through:

- alignment
- duplicate marking
- blacklist filtering
- peak calling
- count matrix generation
- differential analysis

The pipeline produces one filtered BAM per biological replicate.

## When to merge

- **Merge technical replicates:** always merge at the library level before peak calling.
- **Never merge biological replicates** for statistical tests.
- **Merged biological replicate tracks** may be generated for visualization or publication figures, but they should not be treated as independent samples in the differential model.

## Replicate concordance metrics

`ATACseq2tracks` computes replicate-aware QC metrics including:

- FRiP per sample and per-consensus peak set
- peak counts and peak-width distributions
- `deepTools` correlation and PCA across consensus peaks
- TSS enrichment and fragment-size periodicity

Recommended thresholds from ENCODE and best practices:

- TSS enrichment: **>7** ideal for human; **>10** ideal for mouse
- FRiP: **>0.3** preferred; **>0.2** acceptable
- replicate Pearson r: **>0.95** for cell lines; tissue may be lower
- IDR rescue/self-consistency ratio: **<2** if IDR is applied

## Batch and experimental design

For differential accessibility analysis, encode batch structure explicitly in your metadata.

- Use a `batch` column in the samplesheet if batch is a known source of variation.
- Do not include batch covariates that are completely confounded with condition.
- Check PCA before modeling: if samples cluster by batch instead of biology, include batch in the design or redesign the experiment.

### Recommended study designs

| Scenario | Recommended model |
|---|---|
| Simple two-condition | `~ condition` |
| Condition + batch | `~ batch + condition` |
| Paired samples / donor | `~ donor + condition` |
| Multiple covariates | `~ batch + donor + condition + sex` |

> If batch is confounded with condition, statistical correction is not possible. Redesign the experiment or collect additional replicates.

## Notes for ATAC-seq experiments

- Paired-end data is preferred because fragment-size distribution is a core QC metric.
- Single-end data is supported through `SE_SIGNAL_MODE="read"`: one retained alignment is counted once and no physical fragment is inferred.
- Keep PE and SE libraries in separate samplesheets, output directories and DESeq2 cohorts because their signal units differ.
- SE MACS3 shift/extension is a peak-calling model and does not convert normalization signal into observed fragments.
- Tn5 insertion-site tracks are not implemented in this update.
- Use at least **2 biological replicates** for any comparison; **3 or more** is highly recommended.
- Technical replicates increase depth and reproducibility but do not increase inferential power.
