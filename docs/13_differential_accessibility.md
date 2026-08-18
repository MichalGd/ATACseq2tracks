# 13 - Differential accessibility analysis

[Previous: Post-alignment QC](12_post_alignment_qc.md) | [Next: Replicates and design](14_replicates_and_experimental_design.md)

ATACseq2tracks runs two independent differential-accessibility modules. They
share filtered BAMs and sample metadata but construct and analyse their regions
separately. Agreement is useful sensitivity evidence; disagreement is not, by
itself, proof that either method failed.

## Universal condition and replicate policy

The condition names and their number are data-driven; names such as `stem`,
`prolif`, `Day4`, and `Day7` are examples only. For each method and peak type:

1. all non-control biological samples participate in all-sample consensus-peak
   construction, including samples from conditions represented only once;
2. a condition enters the statistical model only when it has at least two
   biological samples after technical rows are collapsed;
3. one multi-condition model is fitted across all eligible conditions;
4. every unordered pair of eligible conditions is extracted as a named contrast.

Thus, `k` eligible conditions produce `k * (k - 1) / 2` comparisons. A condition
with one sample still receives all upstream processing, peaks, tracks, QC and a
place in the consensus universe, but it is not used to estimate differential
model normalization, dispersion or contrasts. When fewer than two conditions
are eligible, differential testing is reported as `SKIPPED`, not failed.

By default, condition order follows first appearance in the samplesheet. For
each pair, the later condition is the numerator and the earlier condition is the
reference. Set a full or partial comma-separated order when a specific direction
is required:

```bash
DIFFERENTIAL_CONDITION_ORDER="stem,prolif,Day4,Day7"
DIFFERENTIAL_MIN_ABS_LOG2FC=0
```

The summary always records numerator, reference and replicate counts; condition
names are never inferred from a fixed list.

## DiffBind

Steps 11 and 12 prepare narrow and broad DiffBind sample sheets and run a
DESeq2-backed DiffBind analysis. The ATAC-specific default is:

```bash
DIFFBIND_SUMMITS=100
```

This passes `summits=100` to `DiffBind::dba.count()` and produces approximately
201-bp summit-centred counting windows. Set `DIFFBIND_SUMMITS=200` to reproduce
the previous approximately 401-bp behavior.

DiffBind first calls `dba.count(minOverlap=2, summits=...)` with every valid
biological sample. The resulting consensus is restricted to canonical autosomes
and X/Y and intervals overlapping the configured blacklist are removed. DiffBind
then recounts only model-eligible samples on that fixed, already recentered
all-sample consensus with `summits=FALSE` and `filter=0`. This keeps
singleton-condition samples in consensus construction without putting them in
the statistical model. Explicit `~Condition` contrasts are added for every
eligible pair and one DESeq2-backed DiffBind analysis is run per peak type.

DiffBind outputs remain under:

```text
diffbind/
diffbind_results/<diffbind_samplesheet>/
|-- diffbind_consensus_peaks.bed
|-- differential_accessibility_condition_eligibility.tsv
|-- differential_accessibility_comparisons.tsv
`-- comparisons/<comparison_id>/
```

Each comparison exports all tested sites, the significant subset, summary/status
text, and PCA, heatmap, MA and volcano plots. `DIFFBIND_ALPHA=0.05` is the
default. Root-level two-condition result names are retained for compatibility.

## DESeq2ATAC

Step 12a implements two separate analyses inspired by the legacy
[MichalGd/ATAC-seq](https://github.com/MichalGd/ATAC-seq) workflows:

1. a broad-peak consensus analysis using every sample's MACS3 `broadPeak` file;
2. a narrow-peak consensus analysis using every sample's MACS3 `narrowPeak` file.

The two analyses are independent. Each constructs its own consensus, count
matrix, DESeq2 size factors, model, results and diagnostic figures. Counts or
normalization factors are never shared between the broad and narrow analyses.

### Exact DESeq2ATAC consensus construction

The following procedure is run separately for broad and narrow peaks:

1. Collapse technical rows to distinct biological sample keys defined as
   `sample_id + replicate`.
2. Import the corresponding per-replicate MACS3 peak file for every non-control
   biological sample, irrespective of differential-model eligibility.
   Pooled-condition peak calls are not used.
3. Retain canonical autosomes and X/Y, harmonize chromosome naming, remove every
   peak overlapping the configured blacklist, and reduce overlapping peaks
   within each sample.
4. Disjoin all sample peak sets into non-overlapping atomic segments. For each
   segment, count how many distinct biological samples have a peak overlapping
   that segment.
5. Keep segments with support greater than or equal to
   `DESEQ2ATAC_MIN_SAMPLES` (default two).
6. Merge adjacent or overlapping retained segments and sort them to produce the
   final non-overlapping consensus BED.
7. Recalculate and report the supporting samples for every final region.

The threshold is enforced on the atomic segments in step 5. After adjacent
supported segments are merged, a sample listed for the final region may support
only part of that merged interval; the reported final support means overlap with
the region, not coverage of its complete width.

This is a **support-filtered consensus**, not a simple union. With the default
threshold, a region called in only one sample is excluded. A condition-specific
region is retained when it is called in at least two samples from that condition;
it does not need to be called in the other condition. Support is evaluated across
the whole cohort and is not a formal within-condition reproducibility or IDR test.

Broad consensus regions retain supported broad-peak boundaries and may therefore
be long and variable in width. Narrow consensus regions retain supported
narrowPeak-derived boundaries; they are not recentered to a fixed summit window.
After consensus construction, the module counts raw overlaps for all samples with
`GenomicAlignments::summarizeOverlaps()`, fits DESeq2 to non-negative integer
counts from model-eligible samples, and writes complete and FDR-filtered results
for every eligible condition pair.

Paired-end input uses strict `GAlignmentPairs` counting: each valid proper pair
is one fragment. Single-end input counts each retained alignment once. The
workflow already rejects mixed PE/SE runs, preventing fragment and read units
from entering one model.

No synthetic GTF or SAF file is required by `summarizeOverlaps`; the stable
BED4 consensus and its support table are the counting annotation. This replaces
the synthetic exon GTF used only to satisfy the legacy HTSeq-count interface.

The consensus-support default is:

```bash
DESEQ2ATAC_MIN_SAMPLES=2
```

Every non-control sample must request `macs2_mode=both` in the samplesheet so
both peak types exist. Preflight rejects incompatible runs before processing.

### How this differs from DiffBind consensus construction

DiffBind also analyzes narrow and broad input peak sets separately and explicitly
uses `minOverlap=2`, equivalent to support from two peak sets. It then calls
`dba.count(summits=100)`, finds a consensus
summit from read pileup and recenters every retained site to a fixed 201-bp
window, as defined by the
[DiffBind `dba.count` documentation](https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf).
DESeq2ATAC does not perform this summit recentering: its broad and narrow models
use the supported interval boundaries described above.

Consequently, even the DESeq2ATAC narrow analysis and DiffBind narrow analysis
do not test identical genomic intervals. Their results compare complete analysis
strategies, not DESeq2 and DiffBind statistics on one shared count matrix.

### Statistical model

The default model and contrast settings are:

```bash
DESEQ2ATAC_ALPHA=0.05
DESEQ2ATAC_BLOCK_COLUMN=""
DESEQ2ATAC_REFERENCE_CONDITION=""
```

The model is `~ condition` and supports any number of eligible conditions. One
model is fitted per broad/narrow universe and all pairwise contrasts are then
extracted from it. DESeq2 size factors and dispersion estimates use eligible
samples only; consensus support and the all-sample raw matrix use every sample.
`DESEQ2ATAC_REFERENCE_CONDITION` is retained for backward compatibility and,
when set, moves that condition to the start of the universal condition order.

Set `DESEQ2ATAC_BLOCK_COLUMN` only for a real paired or batch variable present
in the samplesheet. The module verifies that `~ block + condition` is full rank
and stops if the block is confounded with condition. It never guesses pairing
from replicate numbers.

DESeq2 uses `type="poscounts"` size factors because sparse ATAC peak matrices
often contain zeros. Independent filtering remains enabled. Consequently,
low-information regions can have `padj=NA`; these values are preserved in the
complete result table and never converted to significant values.

### Outputs

```text
deseq2atac/
|-- deseq2atac_peak_type_summary.tsv
|-- broad/
|   |-- deseq2atac_consensus_peaks.bed
|   |-- deseq2atac_consensus_peaks_with_support.tsv.gz
|   |-- deseq2atac_raw_counts.tsv.gz
|   |-- deseq2atac_normalized_counts.tsv.gz
|   |-- differential_accessibility_condition_eligibility.tsv
|   |-- differential_accessibility_comparisons.tsv
|   |-- deseq2atac_summary.txt
|   |-- plots/*.png and plots/*.pdf
|   `-- comparisons/<comparison_id>/
|       |-- deseq2atac_results_all.tsv.gz
|       |-- deseq2atac_results_significant.tsv.gz
|       |-- deseq2atac_summary.txt
|       `-- plots/*.{png,pdf}
`-- narrow/
    `-- (same shared-universe and per-comparison layout)
```

Each peak-type directory also contains sample/library metadata, size factors,
the serialized DESeq2 object, session information, and the full diagnostic set:
library-size, correlation, distance, PCA, dispersion, MA and volcano PNG/PDF
pairs. Shared PCA/correlation/dispersion plots use all model-eligible samples;
MA, volcano and significant-site plots are contrast-specific.

Each significant table is a valid compressed table with a header even when it
has zero rows. A zero-hit broad or narrow analysis exits successfully and states
this explicitly in its comparison summary.

## Built-in peak annotation

`RUN_SIMPLE_PEAK_ANNOTATION=true` annotates the shared broad and narrow universes
once and joins those columns to every complete and significant result table. It
does not replace the separate optional HOMER motif/annotation module.

The GTF supplies a precedence-based promoter/exon/intron/intergenic context,
non-exclusive overlap flags, nearest gene and signed TSS distance. The matching
cCRE BED adds every overlapping ID/class, a deterministic primary cCRE and an
`enhancer_like` flag for pELS/dELS overlaps. These are reference annotations,
not proof of regulatory activity or enhancer-to-gene links. Configure:

```bash
RUN_CCRE_ANNOTATION=true
CCRE_BED_HG38="/home/micgdu/Analysis/utilities/UCSC/CREs/human/hg38/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
CCRE_SOURCE_HG38="ENCODE4_GRCh38_2026"
CCRE_BED_MM39="/home/micgdu/Analysis/utilities/UCSC/CREs/mouse/mm39/encodeCcreCombined_mm39_sorted.bed"
CCRE_SOURCE_MM39="ENCODE3_mm10_liftOver_mm39"
```

With the default `RUN_CCRE_ANNOTATION=true`, a missing/empty reference, unreadable
gzip stream or file without a valid BED interval is a preflight failure. Set
`RUN_CCRE_ANNOTATION=false` explicitly to run GTF-only annotation; all consensus
peaks still receive gene context and nearest TSS fields.

See [Peak annotation](16_peak_annotation.md) for the exact precedence rules,
column definitions, multiple-overlap and primary-cCRE logic, ENCODE class
meanings, human and mouse resource construction, provenance and limitations.

## Counts versus visualization tracks

Neither module uses CPM, FPKM, DESeq2-scaled bigWigs or bedGraphs as statistical
input. Differential testing always starts from raw integer fragment/read counts.
Normalized count tables are descriptive outputs. Length-normalized FPKM values
from the historical scripts are not generated because every sample uses the
same region length for a given row and DESeq2 models raw counts.

## Independent recovery

DiffBind uses `.checkpoints/step12.done`; DESeq2ATAC uses
`.checkpoints/step12a.done`. The workflow attempts both peer modules even if one
fails; broad and narrow wrappers also continue independently. Browser metadata
and the pair-level TSV/HTML report are generated before a failed run exits.
Automatic cleanup is suppressed on failure. Individual comparison status files
distinguish `SUCCESS`, `SKIPPED`, and `FAILED`; zero significant sites is
`SUCCESS`. The current checkpoint granularity remains per module rather than per
contrast, so rerunning a failed module refits its model and rewrites comparisons.

To rerun only DESeq2ATAC:

```bash
rm /path/to/output/.checkpoints/step12a.done
bash atacseq2tracks.sh --config /path/to/config.conf
```

To run the module directly after upstream processing:

```bash
export F2T_CONFIG=/path/to/config.conf
bash scripts/deseq2atac_analysis.sh \
    /path/to/samplesheet.csv \
    /path/to/output/filteredBams \
    /path/to/output/peaks \
    /path/to/output/deseq2atac
```

The wrapper always runs `broad` and `narrow`; each peak type is failure-isolated.

## Interpretation

- Use FDR-adjusted P-values and the explicitly reported contrast direction.
- Inspect FRiP, TSS enrichment, library complexity, PCA and within-condition
  replicate agreement before interpreting significant or null results.
- Two replicates per group meet the software minimum but often provide limited
  power when replicate quality is heterogeneous.
- A consensus supported by two samples is a practical rule, not IDR.
- DESeq2 normalization assumes that most counted regions do not undergo a
  coordinated global accessibility change. Use spike-ins or another external
  reference when absolute global shifts are expected.
- Broad- and narrow-consensus DESeq2ATAC analyses answer related but different
  questions: broad regions favor domain-level changes, while narrow regions give
  more focal resolution and may avoid dilution by unchanged flanking signal.
- Both DESeq2ATAC analyses differ from summit-centred DiffBind and should not be
  expected to return identical adjusted P-values.
