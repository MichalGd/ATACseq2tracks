# 10 — Reference file preparation

[← Troubleshooting](09_troubleshooting.md) | [Back to Index →](README.md)

---

This page describes how to build the reference files required by fastq2tracks that are not available for direct download as pre-built objects.

---

## Bowtie2 indices

```bash
# Download genome FASTA (example: hg38 from UCSC)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Build Bowtie2 index
bowtie2-build --threads 16 hg38.fa /path/to/indices/hg38/hg38

# Same for mm39
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz
bowtie2-build --threads 16 mm39.fa /path/to/indices/mm39/mm39
```

---

## Chromosome sizes

```bash
conda activate fastq2tracks   # fetchChromSizes is included

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

## ChIPQC annotation and blacklist RDS objects

These are pre-built R objects that allow ChIPQC to run without downloading annotation packages at runtime. Build them once and reuse across projects.

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

# Install TxDb if not present
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Build annotation for ChIPQC (GRangesList of features)
anno_hg38 <- ChIPQC:::GetAnnotation(
    annotation  = txdb_hg38,
    chromosomes = paste0("chr", c(1:22, "X", "Y"))
)
saveRDS(anno_hg38, file.path(OUT_DIR, "anno_hg38_chipqc.rds"))
message("  Saved: anno_hg38_chipqc.rds")

# Convert blacklist BED to GRanges RDS
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

message("Done. Set these paths in config.conf:")
message("  CHIPQC_ANNOTATION_HG38=", file.path(OUT_DIR, "anno_hg38_chipqc.rds"))
message("  CHIPQC_ANNOTATION_MM39=", file.path(OUT_DIR, "anno_mm39_chipqc.rds"))
message("  CHIPQC_BLACKLIST_HG38_RDS=", file.path(OUT_DIR, "blacklist_hg38.rds"))
message("  CHIPQC_BLACKLIST_MM39_RDS=", file.path(OUT_DIR, "blacklist_mm39.rds"))
```

Run:
```bash
conda activate fastq2tracks
Rscript build_chipqc_rds.R
```

This takes approximately 10–20 minutes per genome depending on server speed.

---

## Summary: all reference file paths for config.conf

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

# ChIPQC RDS objects
CHIPQC_ANNOTATION_HG38="/path/to/ref/AnnotationHub_cache/anno_hg38_chipqc.rds"
CHIPQC_ANNOTATION_MM39="/path/to/ref/AnnotationHub_cache/anno_mm39_chipqc.rds"
CHIPQC_BLACKLIST_HG38_RDS="/path/to/ref/AnnotationHub_cache/blacklist_hg38.rds"
CHIPQC_BLACKLIST_MM39_RDS="/path/to/ref/AnnotationHub_cache/blacklist_mm39.rds"
```

---

[← Troubleshooting](09_troubleshooting.md) | [Back to Index →](README.md)
