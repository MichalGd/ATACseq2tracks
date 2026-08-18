# 16 - Peak annotation

[Previous: DESeq2ATAC legacy-method review](15_deseq2atac_legacy_method_review.md) | [Back to index](README.md)

ATACseq2tracks adds two complementary annotation layers to differential-analysis
region universes: gene context derived from the configured GTF and regulatory
class derived from a genome-matched ENCODE candidate cis-regulatory element
(cCRE) reference. These annotations describe reference overlaps. They are not
used to fit the statistical models and do not change P-values or fold changes.

## Scope and propagation

With the default settings, annotation is applied once to **every region in each
complete all-sample consensus**, not only to significant regions:

- DiffBind broad and narrow consensuses are annotated independently;
- DESeq2ATAC broad and narrow consensuses are annotated independently;
- the resulting columns are joined by stable region ID to every complete and
  significant pairwise result table; and
- each comparison receives `annotation_summary.tsv`, which counts significant
  sites by primary gene context and primary cCRE class.

The complete universe-level tables are:

```text
diffbind_results/<analysis>/diffbind_consensus_peak_annotations.tsv.gz
deseq2atac/broad/deseq2atac_consensus_peak_annotations.tsv.gz
deseq2atac/narrow/deseq2atac_consensus_peak_annotations.tsv.gz
```

DiffBind pair tables are
`comparisons/<comparison_id>/diffbind_results_all_sites.csv` and
`diffbind_results.csv`. DESeq2ATAC pair tables are
`comparisons/<comparison_id>/deseq2atac_results_all.tsv.gz` and
`deseq2atac_results_significant.tsv.gz`. A two-condition analysis can also have
compatibility copies at the peak-type root; those copies contain the same
annotation columns.

Because each method and peak type constructs a different region universe, do
not assume that annotation rows correspond one-to-one across DiffBind broad,
DiffBind narrow, DESeq2ATAC broad and DESeq2ATAC narrow results.

## GTF-derived gene annotation

The workflow imports `gene` and `exon` features from the selected genome's GTF,
keeps canonical autosomes and X/Y, and harmonizes `chr` naming to the consensus
regions. A one-base transcription start site (TSS) is derived from the
strand-aware start of each `gene` feature. Overlap tests themselves ignore the
peak strand, as ATAC-seq peaks are unstranded.

### Primary gene-context category

Each region receives exactly one `gene_context` according to this precedence:

| Priority | Value | Implemented definition |
|---:|---|---|
| 1 | `promoter` | Overlaps any strand-aware promoter window around a gene TSS |
| 2 | `exon` | Does not overlap a promoter, but overlaps any GTF exon |
| 3 | `intron` | Overlaps a gene body, but no promoter or annotated exon |
| 4 | `other_gene_body` | Fallback gene-body category after the rules above |
| 5 | `distal_intergenic` | Does not overlap any imported gene body |

The default promoter interval is 2,000 bp upstream through 500 bp downstream
of a TSS, in transcriptional orientation. It is controlled by:

```bash
PEAK_ANNOTATION_PROMOTER_UPSTREAM=2000
PEAK_ANNOTATION_PROMOTER_DOWNSTREAM=500
```

Under the current classifier, any gene-body interval that is neither promoter
nor exon is assigned `intron`; consequently `other_gene_body` is a reserved
fallback and will normally not occur with a conventional gene/exon GTF. The
separate Boolean overlap columns remain available so users need not rely only
on the primary label.

### GTF output columns

| Column | Meaning |
|---|---|
| `gene_context` | Single precedence-based category above |
| `promoter_overlap` | Region overlaps one or more configured promoter windows |
| `exon_overlap` | Region overlaps one or more GTF exon features |
| `intron_overlap` | Region overlaps a gene body but no promoter or exon |
| `gene_body_overlap` | Region overlaps one or more GTF gene features |
| `nearest_gene_id` | Gene identifier for the nearest TSS |
| `nearest_gene_name` | Gene name when present, otherwise the gene identifier |
| `nearest_tss_distance_bp` | Peak-midpoint-to-TSS distance in transcriptional orientation |
| `promoter_definition` | Recorded window, for example `TSS-2000/+500bp` |

For `nearest_tss_distance_bp`, zero is the TSS; negative values are upstream and
positive values are downstream relative to transcription. The distance uses the
peak midpoint and is therefore not the minimum edge-to-TSS distance.

Gene context depends on the exact GTF release and on its `gene` and `exon`
features. Overlapping genes can make several Boolean flags true, while the
primary label follows the fixed precedence. “Nearest gene” is a proximity
statement, not evidence that the gene is regulated by the peak.

## cCRE regulatory annotation

The workflow overlaps each consensus region with the configured cCRE BED. It
does not infer a chromatin state from the current ATAC-seq samples. The source
registry has already classified cCREs from integrated biochemical evidence such
as chromatin accessibility, H3K4me3, H3K27ac, CTCF and transcription-factor
binding, together with TSS proximity. The labels are therefore best described
as **candidate regulatory signatures**, not a mutually exclusive chromatin-state
segmentation and not direct proof of activity in the assayed cell type.

### Recognized cCRE labels

The parser supports the current ENCODE4 vocabulary and the legacy ENCODE3 mouse
vocabulary used by the supplied mm39 reference:

| Label | Practical interpretation |
|---|---|
| `PLS` | Promoter-like signature |
| `pELS` | Proximal enhancer-like signature |
| `dELS` | Distal enhancer-like signature |
| `CA-H3K4me3` | Accessible, H3K4me3-associated non-PLS class in ENCODE4 |
| `CA-CTCF` | Accessible and CTCF-associated class in ENCODE4 |
| `CA-TF` | Accessible, TF-bound class lacking the defining promoter/enhancer/CTCF marks |
| `CA` | Chromatin-accessible class without the defining marks above |
| `TF` | TF-bound class with little accessibility or defining histone-mark signal |
| `DNase-H3K4me3` | Legacy ENCODE3 H3K4me3-associated, non-PLS class |
| `CTCF-only` | Legacy ENCODE3 CTCF-associated class |
| `CTCF-bound` | Additional CTCF-bound tag retained when present in a reference class field |

ENCODE4 defines eight principal classes (`PLS`, `pELS`, `dELS`,
`CA-H3K4me3`, `CA-CTCF`, `CA-TF`, `CA`, and `TF`). The expanded registry and
its evidence model are described by Moore, Pratt, Fan, et al.,
[“An expanded registry of candidate cis-regulatory elements”](https://www.nature.com/articles/s41586-025-09909-9),
*Nature* (2026), DOI `10.1038/s41586-025-09909-9`. The current human download
classes are also listed by [SCREEN](https://screen.wenglab.org/downloads).

### Multiple overlaps and the primary cCRE

A consensus region can overlap more than one cCRE, and one cCRE class field can
contain more than one label. The workflow preserves all unique overlaps in
semicolon-delimited fields. It also selects one deterministic primary element:

1. greatest number of overlapping base pairs with the consensus region;
2. if tied, class priority `PLS`, `pELS`, `dELS`, `CA-CTCF`, `CTCF-only`,
   `CTCF-bound`, `DNase-H3K4me3`, `CA-H3K4me3`, `CA-TF`, `CA`, `TF`;
3. if still tied, lexicographically smallest cCRE identifier.

The primary class is a reporting convenience, not a claim that other overlaps
are biologically unimportant. `enhancer_like` is `true` when **any** overlapping
cCRE carries `pELS` or `dELS`, regardless of which element was selected as
primary.

### cCRE output columns

| Column | Meaning |
|---|---|
| `ccre_overlap` | At least one reference cCRE overlaps the region |
| `ccre_primary_id` | Deterministically selected primary cCRE identifier |
| `ccre_primary_class` | Complete class label(s) of the primary cCRE |
| `ccre_all_ids` | All unique overlapping cCRE identifiers, semicolon-delimited |
| `ccre_all_classes` | All unique overlapping class values, semicolon-delimited |
| `enhancer_like` | Any overlap contains `pELS` or `dELS` |
| `ccre_source` | User/configured provenance label copied to every row |

The parser takes the cCRE identifier from BED column 4 and searches columns 4
onward for recognized class labels. This accommodates the native six-column
human BED and the wider legacy UCSC mouse BED.

## Genome-specific sources and construction

### Human hg38

The default is the native GRCh38 ENCODE4 expanded registry file
`Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz` from the official
[supplementary-data directory](https://users.moore-lab.org/ENCODE-cCREs/Supplementary-Data/).
It is not lifted from another assembly. The repository utility
`utilities/prepare_encode4_hg38_ccre.sh` downloads the file, validates its gzip
stream, expected 2,348,854 canonical records, coordinates, accessions and class
vocabulary, verifies the installed SHA-256 checksum, and writes a provenance
TSV. The enforced count is specific to this canonical BED download; broader
registry totals reported elsewhere must not be substituted for this file-level
validation contract.

### Mouse mm39

The configured mm39 file is a legacy ENCODE3 resource, not a native ENCODE4
mm39 registry. The historical
[`creating_ENCODE_cCRES_mm39_bigBed_track.sh`](https://github.com/MichalGd/ATAC-seq/blob/main/utilitiies/creating_ENCODE_cCRES_mm39_bigBed_track.sh)
downloads UCSC's mm10 `encodeCcreCombined.bb`, converts bigBed to BED, applies
the UCSC `mm10ToMm39` chain with `liftOver -bedPlus=9 -tab`, records unmapped
elements separately, and coordinate-sorts the mapped BED. UCSC documents the
source as the archival
[ENCODE3 mouse cCRE track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=mm10&g=encodeCcreCombined),
which contains `PLS`, `pELS`, `dELS`, `DNase-H3K4me3`, and `CTCF-only`
signatures. UCSC now recommends its ENCODE4 track for current analyses; the
workflow retains this lifted ENCODE3 default to match the established server
resource and records that limitation explicitly.

LiftOver can lose or alter intervals, and the mouse and human defaults use
different registry generations and class vocabularies. Their class counts are
not directly comparable. The source labels should therefore remain explicit:

```bash
CCRE_SOURCE_HG38="ENCODE4_GRCh38_2026"
CCRE_SOURCE_MM39="ENCODE3_mm10_liftOver_mm39"
```

## Configuration and opt-out behavior

```bash
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=true
CCRE_BED_HG38="/path/to/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
CCRE_BED_MM39="/path/to/encodeCcreCombined_mm39_sorted.bed"
```

`RUN_SIMPLE_PEAK_ANNOTATION=true` enables the built-in annotation for both
differential modules. With the default `RUN_CCRE_ANNOTATION=true`, preflight
requires a non-empty reference matching the run genome, a readable gzip stream
when compressed, and at least one valid BED interval. The human preparation
utility performs the stricter record-count, accession and class validation
described above; generic preflight does not certify biological provenance.
Set `RUN_CCRE_ANNOTATION=false` on a server without the cCRE resource; GTF gene
context and nearest-TSS annotation still run. Set
`RUN_SIMPLE_PEAK_ANNOTATION=false` only to disable the complete built-in GTF and
cCRE annotation layer. The separate HOMER module remains optional and is not
controlled by these settings.

## Interpretation limits

- A reference overlap is not evidence that a cCRE is active in the assayed
  cells or that a peak regulates its nearest gene.
- Cell type-agnostic registries combine evidence across many biosamples and can
  label elements that are inactive in the current experiment.
- `enhancer_like` means overlap with an ENCODE pELS/dELS signature; it is not a
  functional enhancer assay and does not provide an enhancer-to-gene link.
- Broad peaks can overlap several genes and cCREs; consult all-overlap columns
  and the genomic interval rather than only the primary labels.
- Annotation versions, genome assembly, chromosome filtering and liftOver
  provenance must be recorded when comparing analyses.

[Previous: DESeq2ATAC legacy-method review](15_deseq2atac_legacy_method_review.md) | [Back to index](README.md)
