# DESeq2ATAC legacy-method review

[Previous: Replicates and design](14_replicates_and_experimental_design.md) | [Next: Peak annotation](16_peak_annotation.md)

## Source reviewed

The DESeq2ATAC module was informed by the scripts in
[MichalGd/ATAC-seq](https://github.com/MichalGd/ATAC-seq), especially:

- `workflows/peakSetsAnalysis/reducedAllSamplesPeakSum_bed_gtf_workflow.sh`;
- `workflows/DifferentialAcessibilityAnalysis/countingATACseqReadsAllPeaksSetHTseq_workflow.sh`;
- `workflows/DifferentialAcessibilityAnalysis/makingCountsTables_DESeq2_workflow.sh`;
- `workflows/DifferentialAcessibilityAnalysis/DifferentialChromatinAcessibilityWithDEseq2_workflow.sh`;
- `workflows/DifferentialAcessibilityAnalysis/ATAC_DEresults_overview_dsitances_PCA_heatmapas_MA_Volcano_workflow.sh`.

The historical code is a study record rather than a portable software module.
It contains interactive R sessions, pasted console output, absolute server
paths, hard-coded sample order and contrasts, and manual file-cleaning commands.

## Behavior retained

- MACS broad peaks from all analysed biological samples define the broad-analysis
  region universe; a parallel narrowPeak analysis applies the same modernized
  consensus rule to provide a more focal alternative.
- Overlapping peak intervals are reduced to non-overlapping regions.
- Every sample is counted over the same regions.
- DESeq2 receives raw integer counts and fits condition contrasts.
- Raw and normalized matrices, size factors, full results and FDR-filtered
  results are exported.
- PCA, sample similarity/distance, MA and volcano diagnostics are generated.

## Behavior modernized

### Region support

The legacy script applies `reduce()` to the union of all broad peaks. A peak
called in only one sample can therefore enter the testing universe, and the
documented analysis tested more than 300,000 nonzero regions with a large
low-count fraction. For each peak type independently, DESeq2ATAC disjoins the
sample peak sets into atomic segments, calculates biological-sample support,
keeps atoms supported by at least two samples by default, and reduces adjacent
retained atoms. The threshold is configurable and reported for every consensus
region. Thus neither analysis is a one-sample union.

### Counting engine

The historical implementation uses HTSeq-count with generated GTF exons. PE BAMs
are name-sorted first, and diagnostic lines are removed manually from HTSeq
stdout. DESeq2ATAC instead uses the already installed Bioconductor
`GenomicAlignments::summarizeOverlaps()` Union mode:

- PE: coordinate-sorted filtered BAM, `singleEnd=FALSE`, strict proper-pair
  `GAlignmentPairs` counting, one fragment per pair;
- SE: `singleEnd=TRUE`, one retained alignment per observation;
- `inter.feature=TRUE`, matching HTSeq's conservative treatment of observations
  that overlap more than one feature.

This avoids temporary name-sorted BAMs, manual line deletion and a new HTSeq
dependency. The consensus intervals are non-overlapping, so ambiguity is uncommon
except for fragments spanning two nearby regions.

### Input filtering

The historical region union mixes inputs from separate studies and assumes their
preprocessing is compatible. DESeq2ATAC uses the workflow's filtered BAMs and
applies one genome build, canonical chromosomes, blacklist removal and one
PE-only or SE-only signal definition per run.

### Metadata and contrasts

Hard-coded column indices and condition vectors were replaced with validated
samplesheet metadata. Technical rows are collapsed by `sample_id + replicate`;
biological replicates remain separate. Any number and names of conditions are
accepted. All non-control samples contribute to consensus construction and the
descriptive raw matrix; conditions with fewer than two biological samples are
excluded only from the statistical model. One model is fitted per peak type and
all unordered pairs among eligible conditions are exported. The reference,
contrast direction and design formula are written to each comparison summary.

An optional block column is supported only when explicitly configured and when
the resulting design matrix is full rank. Replicate numbers are not assumed to
represent paired donors.

### DESeq2 input and normalization

The legacy scripts correctly construct DESeq2 from raw HTSeq counts, but also
produce a table labelled FPKM by dividing DESeq2-normalized counts by peak length.
DESeq2ATAC does not use or export FPKM as an inferential input. Per-region length
is constant between samples and is already accounted for by within-row
comparison; applying FPKM before DESeq2 would violate its count model.

The modern module uses `type="poscounts"` size-factor estimation for sparse ATAC
matrices, preserves independent-filtered `padj=NA` results, and writes compressed
large tables.

### Reproducibility and failure handling

The legacy analysis depends on interactive state and manual copying. DESeq2ATAC
is non-interactive, records session information, validates files and design,
produces deterministic names, treats zero significant sites as a successful
result, and has an independent checkpoint.

## Functionality intentionally not transferred

The EDC, Keratin I/II, p63 knockout, promoter and cCRE visualizations are tied to
specific studies, coordinates and sample names. They remain examples of future
optional biological-interpretation extensions and are not part of the generic
workflow module.

The historical scripts also calculate `apeglm`-shrunken fold changes for a few
hard-coded contrasts. This update reports standard DESeq2 maximum-likelihood
fold changes and does not add `apeglm` as a dependency. Generic shrinkage can be
added later only with an explicit coefficient/contrast policy; it is not needed
for the primary Wald-test P-values or adjusted P-values requested here.

The historical use of multiple public studies in one DESeq2 object is also not
reproduced automatically. Cross-study analysis requires harmonized preprocessing,
a valid study/batch design and careful assessment of confounding; it cannot be
made safe by concatenating count matrices alone.

## Remaining limitations

- Peak support of two samples is not a formal reproducibility method such as IDR.
- Consensus regions are derived from the analysed samples, so the tested universe
  remains cohort-dependent.
- Broad merged intervals can vary greatly in width and may dilute focal changes;
  the parallel narrow consensus provides a more localized alternative. Both
  retain supported peak-derived boundaries rather than fixed summit windows.
  Width-normalized values can assist descriptive plots but are not substitutes
  for raw-count modeling.
- DESeq2 size factors cannot establish absolute accessibility when most of the
  genome changes in one direction. External spike-ins or another calibrated
  reference are needed for absolute global shifts.
- Real PE and SE server fixtures are still required to validate runtime,
  performance and output appearance in the production environment.

[Previous: Replicates and design](14_replicates_and_experimental_design.md) | [Next: Peak annotation](16_peak_annotation.md)
