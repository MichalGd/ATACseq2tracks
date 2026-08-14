# 13 - Differential accessibility analysis

[Previous: Post-alignment QC](12_post_alignment_qc.md) | [Next: Replicates and design](14_replicates_and_experimental_design.md)

ATACseq2tracks runs two independent differential-accessibility modules. They
share filtered BAMs and sample metadata but construct and analyse their regions
separately. Agreement is useful sensitivity evidence; disagreement is not, by
itself, proof that either method failed.

## DiffBind

Steps 11 and 12 prepare narrow and broad DiffBind sample sheets and run a
DESeq2-backed DiffBind analysis. The ATAC-specific default is:

```bash
DIFFBIND_SUMMITS=100
```

This passes `summits=100` to `DiffBind::dba.count()` and produces approximately
201-bp summit-centred counting windows. Set `DIFFBIND_SUMMITS=200` to reproduce
the previous approximately 401-bp behavior.

This update does not silently change DiffBind's existing normalization, contrast
or FDR behavior; only the summit width is newly configured.

DiffBind outputs remain under:

```text
diffbind/
diffbind_results/<diffbind_samplesheet>/
```

Each analysis exports all tested sites, the FDR 0.05 subset, its serialized
DiffBind object, summary text, and PCA, heatmap, MA and volcano plots.

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
   biological sample. Pooled-condition peak calls are not used.
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
After consensus construction, the module counts raw overlaps with
`GenomicAlignments::summarizeOverlaps()`, fits DESeq2 to non-negative integer
counts and writes complete and FDR-filtered results.

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

DiffBind also analyzes narrow and broad input peak sets separately and, because
the workflow does not override `minOverlap`, uses the DiffBind default support
of two peak sets. It then calls `dba.count(summits=100)`, finds a consensus
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

The model is `~ condition`. Exactly two conditions and at least two biological
replicates per condition are currently required, matching the supported
DiffBind comparison. When the reference is empty, the alphabetically first
condition is the denominator; the summary records the direction explicitly.

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
|   |-- deseq2atac_results_all.tsv.gz
|   |-- deseq2atac_results_significant.tsv.gz
|   |-- deseq2atac_summary.txt
|   `-- plots/*.png and plots/*.pdf
`-- narrow/
    |-- deseq2atac_consensus_peaks.bed
    |-- deseq2atac_consensus_peaks_with_support.tsv.gz
    |-- deseq2atac_raw_counts.tsv.gz
    |-- deseq2atac_normalized_counts.tsv.gz
    |-- deseq2atac_results_all.tsv.gz
    |-- deseq2atac_results_significant.tsv.gz
    |-- deseq2atac_summary.txt
    `-- plots/*.png and plots/*.pdf
```

Each peak-type directory also contains sample/library metadata, size factors,
the serialized DESeq2 object, session information, and the full diagnostic set:
library-size, correlation, distance, PCA, dispersion, MA and volcano PNG/PDF
pairs. A significant-site overview is present only when significant sites exist.

Each significant table is a valid compressed table with a header even when it
has zero rows. A zero-hit broad or narrow analysis exits successfully and states
this explicitly in its own `deseq2atac_summary.txt`.

## Counts versus visualization tracks

Neither module uses CPM, FPKM, DESeq2-scaled bigWigs or bedGraphs as statistical
input. Differential testing always starts from raw integer fragment/read counts.
Normalized count tables are descriptive outputs. Length-normalized FPKM values
from the historical scripts are not generated because every sample uses the
same region length for a given row and DESeq2 models raw counts.

## Independent recovery

DiffBind uses `.checkpoints/step12.done`; DESeq2ATAC uses
`.checkpoints/step12a.done`. The workflow attempts both peer modules even if one
fails, retains successful outputs, and exits non-zero after both attempts when a
failure occurred.

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

The wrapper always runs `broad` and `narrow`. The underlying R script retains an
optional final `broad|narrow` argument for focused testing; when omitted it
defaults to `broad` for compatibility with the preceding DESeq2ATAC update.

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
