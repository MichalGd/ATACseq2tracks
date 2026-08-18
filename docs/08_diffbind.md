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

# Run differential analysis explicitly with DESeq2
dba <- dba.analyze(dba, method = DBA_DESEQ2)

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

## Automated DiffBind analysis

Steps 11 and 12 run DiffBind directly on the prepared samplesheets.

```bash
export F2T_CONFIG=/absolute/path/to/config.conf
bash scripts/diffbind_analysis.sh diffbind diffbind_results
```

This writes results and diagnostic plots under `diffbind_results/`.

The automated module accepts any number and names of conditions. Every valid
sample contributes to the all-sample consensus. Conditions with fewer than two
biological samples are listed in the eligibility table but excluded from the
statistical model; all remaining condition pairs are analyzed. The directory
layout is:

```text
diffbind_results/diffbind_samplesheet_<genome>_<broad|narrow>/
|-- diffbind_consensus_peaks.bed
|-- diffbind_consensus_peak_annotations.tsv.gz
|-- differential_accessibility_condition_eligibility.tsv
|-- differential_accessibility_comparisons.tsv
`-- comparisons/<comparison_id>/
    |-- diffbind_results_all_sites.csv
    |-- diffbind_results.csv
    |-- diffbind_significant_sites.bed
    |-- annotation_summary.tsv
    |-- diffbind_summary.txt
    `-- diffbind_{pca,heatmap,ma,volcano}.png
```

One `~Condition` model is fitted per broad/narrow analysis and every eligible
pair is added as an explicit contrast. Set `DIFFERENTIAL_CONDITION_ORDER` to
control numerator/reference direction; otherwise first samplesheet appearance
defines the order. A two-condition run also retains the established root-level
result filenames.

Automated ATAC counting uses `DIFFBIND_SUMMITS=100` by default, yielding an
approximately 201-bp summit-centred window. Set `DIFFBIND_SUMMITS=200` in the
configuration to reproduce the previous approximately 401-bp behavior. The
independent broad- and narrow-consensus DESeq2ATAC peer module is documented in
[Differential accessibility](13_differential_accessibility.md).
Annotation columns and reference provenance are documented in
[Peak annotation](16_peak_annotation.md).
