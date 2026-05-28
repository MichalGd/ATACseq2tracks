# 11 — Blacklist Filtering

[← Index](README.md)

---

## What is a genomic blacklist?

A genomic **blacklist** (also called an *exclusion list*) is a curated set of genomic regions that consistently produce anomalous, artifactual, or unstructured signals in next-generation sequencing experiments — regardless of cell type, tissue, antibody, or experimental protocol. These regions are not simply low-complexity or repetitive in the conventional sense; they typically map uniquely, so standard mappability-based filters do not remove them. Instead, they arise from structural deficiencies in the genome assembly: unannotated repeat arrays, collapsed duplications, satellite DNA near centromeres and telomeres, and mitochondrial sequence homologs embedded in the nuclear genome.

Without blacklist filtering, these regions generate spuriously high read pileups that:

- Inflate library size estimates used for normalisation
- Produce false-positive peaks dominating MACS2/3 output
- Confound FRiP score calculations in ChIPQC
- Skew differential binding analysis in DiffBind

Blacklist filtering is therefore considered an **essential quality step** for ChIP-seq, ATAC-seq, CUT&RUN, and CUT&TAG data [1][2].

---

## How blacklist filtering is implemented in fastq2tracks

Blacklist filtering is **Step 6** of the pipeline, applied immediately after Picard deduplication (Step 5) and before all downstream analyses (tracks, peak calling, QC, DiffBind).

### Scripts involved

| Script | Role |
|---|---|
| `scripts/blacklist_filter_batch.sh` | Reads the samplesheet, dispatches one filtering job per sample in parallel |
| `scripts/blacklist_filter.sh` | Single-sample worker: performs the actual filtering |

### What happens for each sample

```
dedupBams/<sample_id>_bioR<N>_dedup.bam
              │
              │  bedtools intersect -v -abam <dedup.bam> -b <blacklist.bed>
              ▼
         [unsorted temp BAM]
              │
              │  samtools sort + samtools index
              ▼
filteredBams/<sample_id>_bioR<N>_dedup_blFilt.bam
              │
              │  samtools view -c (before and after)
              ▼
         Log: "before=X after=Y removed=Z reads"
```

### The core filtering command

```bash
bedtools intersect -v \
    -abam <sample>_dedup.bam \
    -b    <blacklist>.bed \
    > filtered_tmp.bam

samtools sort -@ ${THREADS_SAMTOOLS} filtered_tmp.bam -o <sample>_dedup_blFilt.bam
samtools index -@ ${THREADS_SAMTOOLS} <sample>_dedup_blFilt.bam
```

The `-v` flag inverts the match: **only reads that do NOT overlap any blacklisted region are retained**. This operates at the read level, not the peak level — reads overlapping blacklist regions are removed from the BAM file entirely before any downstream analysis begins.

### Blacklist path: per-sample from the samplesheet

The blacklist BED file path is read from **samplesheet column 16** (`blacklist`), allowing different samples in the same run to use different blacklists (e.g. hg38 and mm39 samples side-by-side):

```
sample_id, ..., blacklist, ...
NHEK_H3K27ac_day0_bioR1, ..., /ref/blacklist_hg38_ENCFF356LFX.bed, ...
Tco_rest_1_SATB1_bioR1,  ..., /ref/blacklist_hg38_ENCFF356LFX.bed, ...
```

If the `blacklist` column is empty for a sample, that sample is skipped with a warning — no filtering is applied.

### Parallelisation

Up to `THREADS_PARALLEL_JOBS` samples are filtered simultaneously. Each sample job writes its own log to:
```
filteredBams/blacklist_logs_<timestamp>/<sample_id>.log
```

The main log records before/after read counts for every sample.

---

## What the blacklist removes

Blacklisted regions fall into several categories, all characterised by producing anomalously high, cell-type-independent read signal:

| Region type | Why it causes artifacts |
|---|---|
| **Satellite repeats** (centromeric, pericentromeric) | Tandem repeat arrays with high sequence similarity; reads pile up due to multi-mapping and assembly gaps |
| **Telomeric repeats** | TTAGGG arrays at chromosome ends; incompletely assembled and generate spurious signal |
| **rDNA loci** (ribosomal DNA) | High-copy tandem arrays on acrocentric chromosomes; largely absent from reference assemblies |
| **Mitochondrial nuclear insertions (NUMTs)** | Fragments of mtDNA integrated in the nuclear genome; mtDNA is vastly over-represented in ChIP inputs |
| **Collapsed/duplicated assembly regions** | Genomic segments where the assembler merged near-identical sequences into one, causing artificial read depth doubling |
| **Low-mappability regions** | Short unique k-mers that nevertheless receive multi-mapper overflow |
| **"Hyper-chippable" regions** | Regions that immunoprecipitate non-specifically with virtually any antibody, regardless of the target |

In practical terms, for a typical ChIP-seq library aligned to hg38, blacklist removal affects **< 1% of reads** in high-quality experiments but can remove the top-scoring artifactual "peaks" that would otherwise dominate the MACS2/3 output and FRiP calculation.

---

## Blacklist sources and references

### hg38 — ENCODE unified exclusion list (ENCFF356LFX) ✓ **Used in fastq2tracks**

| Attribute | Detail |
|---|---|
| File accession | [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/) |
| Assembly | GRCh38 / hg38 |
| Regions | 910 |
| Lab | Anshul Kundaje, Stanford University |
| Project | ENCODE4 |
| Release date | 2020-05-05 |
| Download | `https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz` |

This is the **current recommended blacklist for hg38** and is referred to as the "Kundaje unified" or "GRCh38 unified blacklist". It was generated by analysing input control tracks from hundreds of ENCODE ChIP-seq experiments, identifying 1 kb windows with read depths or multi-mapping rates in the top 1% across all samples. The list is manually curated and supersedes earlier ENCODE DAC and Duke exclusion lists [1][3].

**Cite as:** Amemiya et al. 2019 (see below) + ENCODE accession ENCFF356LFX.

---

### mm39 / GRCm39 — excluderanges (Dozmorov lab) ✓ **Used in fastq2tracks**

| Attribute | Detail |
|---|---|
| Tool / package | [excluderanges](https://github.com/dozmorovlab/excluderanges) (Bioconductor) |
| Assembly | GRCm39 / mm39 |
| Method | Boyle-Lab Blacklist software applied to mm39 |
| Available via | BEDbase API / Bioconductor `excluderanges` package |
| Reference | Ogata et al. 2023 (see below) |

As of 2022–2023, no official ENCODE blacklist exists for mm39. The `excluderanges` package (Ogata et al. 2023) provides the most systematic option, generated using the Boyle-Lab Blacklist algorithm applied to the mm39 assembly [4]. Download via:

```r
# In R:
BiocManager::install("excluderanges")
library(excluderanges)
# Or download BED directly from BEDbase:
# bedbase_id: edc716833d4b5ee75c34a0692fc353d5
```

Or as a BED file:
```bash
wget -O mm39.excluderanges.bed.gz \
    "http://bedbase.org/api/bed/edc716833d4b5ee75c34a0692fc353d5/file/bed"
gunzip mm39.excluderanges.bed.gz
```

**Alternative for mm39:** Lift over the official mm10 Boyle-lab blacklist (`mm10-blacklist.v2.bed`) using UCSC liftOver. Note that ENCODE explicitly cautions that liftover is not ideal for blacklists, because improved assemblies may resolve regions that were problematic in an earlier assembly [2][5].

---

### The Boyle-Lab/ENCODE Blacklist algorithm (v2)

The blacklists for hg38, hg19, mm10, dm6, ce11 were generated by the **Boyle Lab** blacklist pipeline and published in 2019 [1]:

1. Input control tracks (fragmented genomic DNA / IgG) from all available ENCODE experiments are collected
2. The genome is tiled in 1 kb windows with 100 bp overlap
3. Windows are scored by input read depth and mappability; scores are quantile-normalised across all input samples
4. Windows in the **top 1% of read depth or multi-mapping rate** are flagged
5. Flagged regions are merged and manually curated

The key insight is that these regions are defined **cell-type agnostically** — a region is blacklisted only if it shows artifact signal across many independent experiments, not just in one sample.

---

## Primary references

**[1]** Amemiya HM, Kundaje A, Boyle AP. **The ENCODE Blacklist: Identification of Problematic Regions of the Genome.** *Scientific Reports* 9, 9354 (2019).
DOI: [10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z) · PMID: 31249361 · PMCID: PMC6597582
→ *Primary reference for the Boyle-Lab/ENCODE blacklist algorithm and all v2 blacklist files (hg38, mm10, dm6, ce11)*

**[2]** ENCODE Project Consortium. **An integrated encyclopedia of DNA elements in the human genome.** *Nature* 489, 57–74 (2012).
DOI: [10.1038/nature11247](https://doi.org/10.1038/nature11247)
→ *Original ENCODE project; source of input control experiments used to build the blacklist*

**[3]** ENCODE accession **ENCFF356LFX** — GRCh38 unified exclusion list regions.
[https://www.encodeproject.org/files/ENCFF356LFX/](https://www.encodeproject.org/files/ENCFF356LFX/)
Lab: Anshul Kundaje, Stanford. Released 2020-05-05. ENCODE4 project.

**[4]** Ogata JD, Mu W, Davis ES, Xue B, Harrell JC, Sheffield NC, Phanstiel DH, Love MI, Dozmorov MG. **excluderanges: exclusion sets for T2T-CHM13, GRCm39, and other genome assemblies.** *Bioinformatics* 39(4): btad198 (2023).
DOI: [10.1093/bioinformatics/btad198](https://doi.org/10.1093/bioinformatics/btad198) · PMID: 37067481
→ *Reference for mm39 blacklist; provides unified access to exclusion sets including GRCm39 via Bioconductor `excluderanges`*

**[5]** Carroll TS, Liang Z, Salama R, Stark R, de Santiago I. **Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data.** *Frontiers in Genetics* 5:75 (2014).
DOI: [10.3389/fgene.2014.00075](https://doi.org/10.3389/fgene.2014.00075)
→ *Demonstrates quantitative impact of blacklist removal on ChIP-seq quality metrics including cross-correlation and FRiP*

---

## Quick reference: blacklist files used in this workflow

| Genome | File | Regions | Config variable |
|---|---|---|---|
| hg38 | ENCFF356LFX | 910 | `BLACKLIST_HG38` |
| mm39 | excluderanges mm39 | ~2,300 | `BLACKLIST_MM39` |

Both are also specified **per sample** in samplesheet column 16 (`blacklist`), which takes precedence for mixed-genome runs.

---

[← Index](README.md) | [← Troubleshooting](09_troubleshooting.md)
