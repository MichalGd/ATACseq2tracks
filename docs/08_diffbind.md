# 08 — Downstream: DiffBind

[← Outputs](07_outputs.md) | [Next: Troubleshooting →](09_troubleshooting.md)

---

After the pipeline completes, two DiffBind-compatible samplesheets per genome are ready in `diffbind/`. These contain all information DiffBind needs — BAM paths, peak paths, sample metadata — pre-filled from the pipeline samplesheet.

---

## Samplesheet columns

| Column | Source |
|---|---|
| `SampleID` | `sample_id` from samplesheet |
| `Tissue` | `cell_type` |
| `Factor` | `factor` |
| `Condition` | `condition` |
| `Treatment` | `treatment` |
| `Replicate` | `replicate` |
| `bamReads` | Path to filtered BAM |
| `bamControl` | Path to matched control BAM |
| `Peaks` | Path to narrow or broad peak file |
| `PeakCaller` | `narrow` or `broad` |

---

## Basic DiffBind workflow

```r
library(DiffBind)

# Choose the appropriate samplesheet
# For sharp marks (H3K27ac, H3K4me3, CTCF, p63, SATB1):
ss <- "diffbind/diffbind_samplesheet_hg38_narrow.csv"

# For broad marks (H3K27me3, H3K9me3, H3K36me3):
# ss <- "diffbind/diffbind_samplesheet_hg38_broad.csv"

# Load and count reads in peaks
dba <- dba(sampleSheet = ss)
dba <- dba.count(dba, bUseSummarizeOverlaps = TRUE)

# Normalise
dba <- dba.normalize(dba)

# Set contrasts (by condition, treatment, or manually)
dba <- dba.contrast(dba, categories = DBA_CONDITION, minMembers = 2)

# Run differential analysis (DESeq2 by default)
dba <- dba.analyze(dba)

# Extract results
res <- dba.report(dba, th = 0.05)

# Export
export.bed(res, "diffbind_results.bed")
write.csv(as.data.frame(res), "diffbind_results.csv")
```

---

## Visualisation

```r
# Correlation heatmap
dba.plotHeatmap(dba)

# PCA
dba.plotPCA(dba, attributes = DBA_CONDITION, label = DBA_ID)

# MA plot
dba.plotMA(dba)

# Volcano plot
dba.plotVolcano(dba)
```

---

## Tips

- **Narrow vs broad samplesheet:** Always use the narrow samplesheet for sharp marks even if you ran `macs2_mode=both`. Use the broad samplesheet for H3K27me3, H3K9me3, H3K36me3.
- **Minimum replicates:** DiffBind requires at least 2 replicates per condition for statistical testing.
- **Normalisation:** `dba.normalize()` defaults to library-size normalisation. For ChIP-seq with spike-in controls, use `method = DBA_NORM_SPIKEIN`.
- **Missing peaks warning:** DiffBind will warn about samples with no peaks in a consensus set — check `chipqc_frip.csv` for those samples.

---

[← Outputs](07_outputs.md) | [Next: Troubleshooting →](09_troubleshooting.md)
