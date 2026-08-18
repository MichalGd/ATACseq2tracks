# 10 — Reference file preparation

[← Troubleshooting](09_troubleshooting.md) | [Back to Index →](README.md)

---

This page describes how to build the reference files required by ATACseq2tracks.

---

## Bowtie2 indices

```bash
# Download genome FASTA (example: hg38 from UCSC)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build Bowtie2 index (16 threads, ~45 min for hg38)
bowtie2-build --threads 16 hg38.fa /path/to/indices/hg38/hg38

# Same for mm39
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz
bowtie2-build --threads 16 mm39.fa /path/to/indices/mm39/mm39
```

---

## Chromosome sizes

```bash
conda activate ATACseq2tracks   # fetchChromSizes is included

fetchChromSizes hg38 > /path/to/ref/hs38n.chrom.sizes
fetchChromSizes mm39 > /path/to/ref/mm39n.chrom.sizes
```

---

## Blacklist BED files

### hg38 (ENCODE ENCFF356LFX)

```bash
wget -O /path/to/ref/blacklist_hg38_ENCFF356LFX.bed.gz \
    https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip /path/to/ref/blacklist_hg38_ENCFF356LFX.bed.gz
```

### mm39 (Boyle lab)

```bash
wget -O /path/to/ref/blacklist_mm39.bed.gz \
    https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm39-blacklist.v2.bed.gz
gunzip /path/to/ref/blacklist_mm39.bed.gz
```

---

## GTF and cCRE annotations

The built-in simple peak annotation requires the matching hg38 or mm39 GTF.
cCRE classification is enabled by default with `RUN_CCRE_ANNOTATION=true`, so
the selected genome's cCRE BED must exist and pass preflight. A server without
this reference can explicitly set `RUN_CCRE_ANNOTATION=false` for GTF-only
annotation.

For the exact annotation categories, precedence, table columns and biological
limitations, see [Peak annotation](16_peak_annotation.md). This page focuses on
preparing the reference files.

For human hg38, use the native GRCh38 ENCODE4 expanded Registry of cCREs. The
stable comprehensive BED is available from the
[ENCODE cCRE supplementary-data directory](https://users.moore-lab.org/ENCODE-cCREs/Supplementary-Data/)
as `Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz`; the current registry is also
described on the [SCREEN downloads page](https://screen.wenglab.org/downloads).
Reference: Moore JE, Pratt HE, Fan K, et al.,
["An expanded registry of candidate cis-regulatory elements"](https://www.nature.com/articles/s41586-025-09909-9),
*Nature* (2026), DOI `10.1038/s41586-025-09909-9`.

The repository includes a preparation utility analogous to the historical
mouse cCRE utility, but no liftOver is performed because this source is already
native GRCh38:

```bash
cd /home/micgdu/Analysis/workflows/ATACseq2tracks
bash utilities/prepare_encode4_hg38_ccre.sh
```

It downloads the official compressed six-column BED, requires the published
2,348,854 records, validates canonical chromosomes, coordinates, accessions and
cCRE classes, verifies gzip/checksum integrity, installs the file atomically,
and writes a provenance TSV. Its default output exactly matches the configured
path:

```text
/home/micgdu/Analysis/utilities/UCSC/CREs/human/hg38/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz
```

On another server, choose the destination explicitly:

```bash
bash utilities/prepare_encode4_hg38_ccre.sh \
  --output-dir /shared/references/hg38/ccre
```

Then set `CCRE_BED_HG38` to the reported BED path. Existing files are validated
and retained; use `--force` only to deliberately redownload and replace one.
The installed columns are `chrom`, `chromStart`, `chromEnd`, two ENCODE
accessions, and cCRE class. This is directly compatible with the workflow
parser, which uses BED column 4 as the element ID and detects the recognized
class vocabulary in the remaining columns; for this human file the class is in
column 6.

For mouse mm39, reproduce the existing
[ATAC-seq utility](https://github.com/MichalGd/ATAC-seq/blob/main/utilitiies/creating_ENCODE_cCRES_mm39_bigBed_track.sh):
download ENCODE3 `encodeCcreCombined.bb` for mm10, convert to BED, and lift it
to mm39. The expected server BED is:

```text
/home/micgdu/Analysis/utilities/UCSC/CREs/mouse/mm39/encodeCcreCombined_mm39_sorted.bed
```

Record this as `ENCODE3_mm10_liftOver_mm39`, not as native mm39 ENCODE4. Human
and mouse class totals are not directly comparable because the registry versions
and coordinate-generation methods differ.

```bash
RUN_CCRE_ANNOTATION=true
CCRE_BED_HG38="/home/micgdu/Analysis/utilities/UCSC/CREs/human/hg38/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
CCRE_SOURCE_HG38="ENCODE4_GRCh38_2026"
CCRE_BED_MM39="/home/micgdu/Analysis/utilities/UCSC/CREs/mouse/mm39/encodeCcreCombined_mm39_sorted.bed"
CCRE_SOURCE_MM39="ENCODE3_mm10_liftOver_mm39"
```

---

## Summary: required config.conf variables

```bash
# Bowtie2 indices
INDEX_HG38="/path/to/indices/hg38/hg38"
INDEX_MM39="/path/to/indices/mm39/mm39"

# Chromosome sizes
CHROM_SIZES_HUMAN="/path/to/ref/hs38n.chrom.sizes"
CHROM_SIZES_MOUSE="/path/to/ref/mm39n.chrom.sizes"

# Blacklist BED files
BLACKLIST_HG38="/path/to/ref/blacklist_hg38_ENCFF356LFX.bed"
BLACKLIST_MM39="/path/to/ref/blacklist_mm39.bed"

# Gene annotations (required for TSS QC and built-in peak annotation)
GTF_HUMAN="/path/to/ref/gencode.hg38.annotation.gtf"
GTF_MOUSE="/path/to/ref/gencode.mm39.annotation.gtf"

# deepTools QC threads (Step 10) — required as of v3.1.0
THREADS_DEEPTOOLS=16
```

The matching cCRE BED, Bowtie2 index, chromosome sizes, blacklist and GTF are
required by the default v3.2.0 configuration. Set `RUN_CCRE_ANNOTATION=false`
to remove only the cCRE requirement while retaining GTF annotation.

---

## Legacy: ChIPQC RDS annotation objects

> **These files are no longer required for the main pipeline as of v3.1.0.**
> Step 10 (post-alignment QC) now uses deepTools and does not need R annotation objects.
>
> The build instructions below are historical only. `run_chipqc.R` is not distributed in v3.2.0; use an older Git tag to reproduce v3.0.x ChIPQC results.

### Build script

Save as `build_chipqc_rds.R` and run with `Rscript build_chipqc_rds.R`:

```r
suppressPackageStartupMessages({
    library(BiocManager)
    library(rtracklayer)
    library(GenomicFeatures)
    library(ChIPQC)
})

OUT_DIR <- "/path/to/ref/AnnotationHub_cache"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── hg38 ──────────────────────────────────────────────────────────────────────
message("Building hg38 annotation...")

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

anno_hg38 <- ChIPQC:::GetAnnotation(
    annotation  = txdb_hg38,
    chromosomes = paste0("chr", c(1:22, "X", "Y"))
)
saveRDS(anno_hg38, file.path(OUT_DIR, "anno_hg38_chipqc.rds"))
message("  Saved: anno_hg38_chipqc.rds")

bl_hg38 <- rtracklayer::import("/path/to/ref/blacklist_hg38_ENCFF356LFX.bed",
                                format = "BED")
saveRDS(bl_hg38, file.path(OUT_DIR, "blacklist_hg38.rds"))
message("  Saved: blacklist_hg38.rds")

# ── mm39 ──────────────────────────────────────────────────────────────────────
message("Building mm39 annotation...")

if (!requireNamespace("TxDb.Mmusculus.UCSC.mm39.refGene", quietly = TRUE))
    BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")

library(TxDb.Mmusculus.UCSC.mm39.refGene)
txdb_mm39 <- TxDb.Mmusculus.UCSC.mm39.refGene

anno_mm39 <- ChIPQC:::GetAnnotation(
    annotation  = txdb_mm39,
    chromosomes = paste0("chr", c(1:19, "X", "Y"))
)
saveRDS(anno_mm39, file.path(OUT_DIR, "anno_mm39_chipqc.rds"))
message("  Saved: anno_mm39_chipqc.rds")

bl_mm39 <- rtracklayer::import("/path/to/ref/blacklist_mm39.bed",
                                format = "BED")
saveRDS(bl_mm39, file.path(OUT_DIR, "blacklist_mm39.rds"))
message("  Saved: blacklist_mm39.rds")

message("Done.")
```

Run:
```bash
conda activate ATACseq2tracks
Rscript build_chipqc_rds.R
```

This takes approximately 10–20 minutes per genome depending on server speed.

### Legacy config variables

Historical v3.0.x configurations used these variables:

```bash
CHIPQC_ANNOTATION_HG38="/path/to/ref/AnnotationHub_cache/anno_hg38_chipqc.rds"
CHIPQC_ANNOTATION_MM39="/path/to/ref/AnnotationHub_cache/anno_mm39_chipqc.rds"
CHIPQC_BLACKLIST_HG38_RDS="/path/to/ref/AnnotationHub_cache/blacklist_hg38.rds"
CHIPQC_BLACKLIST_MM39_RDS="/path/to/ref/AnnotationHub_cache/blacklist_mm39.rds"
```

---

[← Troubleshooting](09_troubleshooting.md) | [Back to Index →](README.md)
