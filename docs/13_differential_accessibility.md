# 13 — Differential accessibility analysis

[← Post-alignment QC](12_post_alignment_qc.md) | [Next: Replicates and design →](14_replicates_and_experimental_design.md)

---

This page describes how `ATACseq2tracks` supports downstream differential accessibility analysis through DiffBind and consensus peak-based count matrices.

## What the pipeline delivers

The pipeline produces:

- `diffbind/diffbind_samplesheet_<genome>_narrow.csv`
- `diffbind/diffbind_samplesheet_<genome>_broad.csv`
- `diffbind_results/<diffbind_samplesheet>/*` (DiffBind analysis outputs)

The new `scripts/diffbind_analysis.sh` wrapper runs DiffBind on the prepared samplesheets and writes:

- `diffbind_results/<sample_sheet>/diffbind_results.csv`
- `diffbind_results/<sample_sheet>/diffbind_summary.txt`
- `diffbind_results/<sample_sheet>/diffbind_pca.png`
- `diffbind_results/<sample_sheet>/diffbind_heatmap.png`
- `diffbind_results/<sample_sheet>/diffbind_ma.png`
- `diffbind_results/<sample_sheet>/diffbind_volcano.png`
- `diffbind_results/<sample_sheet>/diffbind_consensus_peaks.bed`

---

## Recommended differential analysis strategy

1. Use the **same consensus peak set** for all biological replicates.
2. Count raw fragments per peak per biological replicate.
3. Use a negative-binomial model such as **DESeq2** or **edgeR**.
4. Include batch covariates when they are not confounded with condition.
5. Do not treat technical replicates as independent biological samples.

This pipeline uses DiffBind as a supported wrapper for peak-based differential analysis.

## Experimental design and design formulas

For a simple two-condition experiment with no batch structure:

```r
design(dds) <- ~ condition
```

For a design with a batch variable that is not confounded with condition:

```r
design(dds) <- ~ batch + condition
```

For paired or donor-matched designs:

```r
design(dds) <- ~ donor + condition
```

DiffBind will infer the `Condition` column from the prepared samplesheet.

## Normalization and count matrices

The pipeline produces normalized `BigWig` tracks for visualization, but the statistical model is built on raw counts.

Recommended normalization choices:

- `DESeq2` / `DiffBind` uses library-size (size factors) normalization.
- `edgeR` uses TMM normalization.
- `limma-voom` uses voom precision weights after TMM or library-size normalization.

If your experiment is expected to have global chromatin shifts, consider spike-in normalization or a background-normalized method such as `csaw`.

## Important considerations

- **Peak set consistency:** use the same consensus peak BED for all samples.
- **Technical replicates:** merge immediately before peak calling and count using the merged library.
- **Biological replicates:** keep separate through peak calling, counting, and statistical modeling.
- **Multiple testing:** report FDR-adjusted p-values (Benjamini–Hochberg) and log fold changes.
- **Fold-change shrinkage:** use shrinkage methods for stable log fold-change estimates if downstream interpretation or ranking is important.

## Running the analysis

To run the differential analysis stage from the pipeline:

```bash
bash ATACseq2tracks/atacseq2tracks.sh --config /path/to/config.conf
```

This now includes Step 12, which runs DiffBind on the prepared samplesheets.

To run only the DiffBind analysis later:

```bash
bash scripts/diffbind_analysis.sh diffbind diffbind_results
```

---

## Interpretation of results

Key outputs to inspect:

- `diffbind_results.csv`: differential accessibility table
- `diffbind_pca.png`: condition separation and batch structure
- `diffbind_heatmap.png`: sample similarity
- `diffbind_ma.png`: global signal and differential distribution
- `diffbind_volcano.png`: significance vs fold change

Always compare differential results with QC metrics such as FRiP, TSS enrichment, and replicate correlation before drawing biological conclusions.
