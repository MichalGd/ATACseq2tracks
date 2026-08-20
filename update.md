# ATACseq2tracks v4.2.0 release update

**Snapshot:** ATACseq2tracks 4.2.0 optional Drosophila spike-in release candidate

**Update date:** 2026-08-20

Version 4.2.0 preserves the five v4.1.0 coverage families and all upstream,
peak, QC, annotation and differential-analysis behavior. It adds one disabled-
by-default stringent dm6 calibration family based on a competitive host+dm6
alignment. Raw host coverage is scaled by a fixed retained-dm6 target, adjusted
for the declared pre-tagmentation dm6:host input ratio; it is not combined with
DESeq2 robust-CPM.

See
[`docs/v4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md`](docs/v4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md)
for methods, inputs, formula, outputs, QC and limitations, and
[`APPLY_V4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md`](APPLY_V4.2.0_DROSOPHILA_SPIKEIN_UPDATE_2026-08-20.md)
for the immutable all-user deployment procedure.

## Inherited v4.1.0 release notes

**Snapshot:** ATACseq2tracks 4.1.0 filtering-policy coverage release candidate

**Update date:** 2026-08-20

**Scope:** Five configurable default coverage families, policy-specific
consensus-count normalization and read-filtering sensitivity, layered on the
unchanged v4.0.0 multi-condition, annotation, QC and reporting workflow.

The inherited filtering-policy release is **4.1.0**. See
[`docs/v4.1.0_READ_FILTERING_COVERAGE_UPDATE_2026-08-20.md`](docs/v4.1.0_READ_FILTERING_COVERAGE_UPDATE_2026-08-20.md)
for filtering semantics, formulas, switches, outputs and validation, and
[`APPLY_V4.1.0_READ_FILTERING_COVERAGE_UPDATE_2026-08-20.md`](APPLY_V4.1.0_READ_FILTERING_COVERAGE_UPDATE_2026-08-20.md)
for versioned shared deployment.

The 2026-08-19 IRanges namespace correction and server-resume instructions are
documented in
[`docs/v3.2.0_IRANGES_NAMESPACE_HOTFIX_2026-08-19.md`](docs/v3.2.0_IRANGES_NAMESPACE_HOTFIX_2026-08-19.md).
The universal implementation details and validation handoff remain in
[`docs/v3.2.0_UNIVERSAL_DIFFERENTIAL_ANNOTATION_UPDATE_2026-08-18.md`](docs/v3.2.0_UNIVERSAL_DIFFERENTIAL_ANNOTATION_UPDATE_2026-08-18.md).

The bounded-parallelism implementation and server deployment notes are in
[`docs/v4.0.0_RELEASE_2026-08-19.md`](docs/v4.0.0_RELEASE_2026-08-19.md).

The current reporting and DiffBind input-filtering correction is documented in
[`docs/v4.0.0_REPORTING_DIFFBIND_HOTFIX_2026-08-20.md`](docs/v4.0.0_REPORTING_DIFFBIND_HOTFIX_2026-08-20.md).

## Universal comparison extension (2026-08-18)

- Any number and names of conditions are accepted.
- All non-control biological samples participate in consensus construction.
- Conditions with fewer than two biological samples are excluded only from
  statistical models and contrasts.
- Every pair among eligible conditions is exported from one model per method
  and peak type.
- DiffBind broad/narrow, DESeq2ATAC broad/narrow, and peer modules are
  failure-isolated.
- Pair-level TSV/HTML reporting distinguishes success, zero-hit, skipped and
  failed analyses and records both directions of significant-site counts.
- Built-in GTF/default cCRE annotation is added to shared universes and
  pair-level result tables.
- A native-GRCh38 ENCODE4 cCRE download/validation utility installs the default
  human reference with checksum and provenance metadata; no liftOver is used.
- Existing safe post-success cleanup is enabled by default; it can be disabled
  with `ENABLE_AUTOMATIC_CLEANUP=false`, and filtered BAMs remain retained.

## Implemented changes

- DESeq2ATAC now runs two independent analyses:
  - a MACS3 broadPeak-based analysis;
  - a MACS3 narrowPeak-based analysis.
- Each analysis independently constructs its consensus regions, counts fragments
  or reads, estimates DESeq2 size factors, fits the statistical model and writes
  results and diagnostic figures.
- Outputs are separated under:
  - `deseq2atac/broad/`;
  - `deseq2atac/narrow/`.
- `deseq2atac/deseq2atac_peak_type_summary.tsv` compares the number of consensus
  regions, nonzero tested regions and significant regions between both analyses.
- Preflight requires `macs2_mode=both` for every non-control biological sample
  when `RUN_DESEQ2ATAC=true`.
- The final differential-accessibility report includes separate broad- and
  narrow-consensus DESeq2ATAC entries.
- Existing workflow architecture, upstream processing, track generation,
  normalization, QC and dependencies remain unchanged.

## Consensus-peak construction

The following procedure is performed separately for broad and narrow peaks:

1. Technical rows are collapsed to distinct biological sample keys defined by
   `sample_id + replicate`.
2. The corresponding per-replicate MACS3 peak file is imported for every
   non-control biological sample. Pooled peak calls are not used.
3. Peaks are restricted to canonical autosomes and X/Y, chromosome naming is
   harmonized, blacklist-overlapping peaks are removed, and overlapping peaks
   within each sample are reduced.
4. All sample peak sets are disjoined into non-overlapping atomic segments.
5. Each segment receives a support count equal to the number of distinct
   biological samples whose peak set overlaps it.
6. Segments supported by at least `DESEQ2ATAC_MIN_SAMPLES` samples are retained;
   the default threshold is two.
7. Adjacent or overlapping retained segments are reduced into the final sorted,
   non-overlapping consensus regions.

This is a **support-filtered consensus**, not a simple union. A region detected
in only one sample is excluded by default. A condition-specific region remains
eligible when it is supported by at least two samples from that condition; it
does not need to be detected in both conditions.

Broad consensus regions can be long and variable in width. Narrow consensus
regions are more focal but retain narrowPeak-derived boundaries. They are not
recentered to a fixed summit window.

## Comparison with DiffBind

The workflow runs DiffBind separately on broad and narrow inputs. DiffBind uses
an explicit `minOverlap=2` all-sample consensus support rule and the configured
`DIFFBIND_SUMMITS=100`, which recenters retained sites around consensus summits
to approximately 201-bp windows.

DESeq2ATAC does not perform summit recentering. Therefore, DESeq2ATAC narrow and
DiffBind narrow do not test identical intervals. Differences between their
results reflect both region construction and downstream analysis, rather than
only a comparison of statistical packages.

### Clarification of `DIFFBIND_SUMMITS`

The current workflow deliberately sets `DIFFBIND_SUMMITS=100`; this is not an
accidental implementation error. DiffBind interprets `summits` as the number of
bases retained on each side of the consensus summit, so:

- `DIFFBIND_SUMMITS=100` produces 201-bp windows;
- `DIFFBIND_SUMMITS=200` produces 401-bp windows.

The parameter is configurable, but the current workflow default of 100 overrides
the upstream DiffBind package default of 200. A 100-bp half-width provides a more
focal ATAC-seq window, whereas 200 includes more flanking signal and may increase
counts while potentially diluting focal accessibility changes or combining nearby
sites. If the intended workflow policy is to follow the official DiffBind default
and restore the preceding 401-bp behavior, the configuration, script fallbacks,
tests and documentation should all be changed consistently to
`DIFFBIND_SUMMITS=200`. This snapshot does not make that change.

Reference: [DiffBind package manual](https://bioconductor.org/packages/release/bioc/manuals/DiffBind/man/DiffBind.pdf).

## Output layout

```text
deseq2atac/
|-- deseq2atac_peak_type_summary.tsv
|-- broad/
|   |-- deseq2atac_consensus_peaks.bed
|   |-- deseq2atac_consensus_peaks_with_support.tsv.gz
|   |-- deseq2atac_raw_counts.tsv.gz
|   |-- deseq2atac_normalized_counts.tsv.gz
|   |-- deseq2atac_consensus_peak_annotations.tsv.gz
|   |-- differential_accessibility_condition_eligibility.tsv
|   |-- differential_accessibility_comparisons.tsv
|   |-- deseq2atac_summary.txt
|   |-- plots/
|   `-- comparisons/<comparison_id>/
|       |-- deseq2atac_results_all.tsv.gz
|       |-- deseq2atac_results_significant.tsv.gz
|       |-- annotation_summary.tsv
|       |-- deseq2atac_summary.txt
|       `-- plots/
`-- narrow/
    `-- (same shared-universe and per-comparison layout)
```

Both peak-type directories also contain sample and library metadata, DESeq2 size
factors, the serialized analysis object, session information and PNG/PDF
diagnostics. In a two-condition analysis, selected root-level result copies are
also retained for backward compatibility; canonical results remain under
`comparisons/<comparison_id>/`.

See [`docs/16_peak_annotation.md`](docs/16_peak_annotation.md) for the exact GTF
and cCRE category definitions, columns, reference construction and limitations.

## Validation completed

- Bash syntax checks passed for all workflow and test scripts.
- Focused dual-analysis regression checks passed.
- Mocked broad and narrow wrapper execution passed.
- `scripts/deseq2atac_analysis.R` passed R syntax parsing.
- Modified files use LF line endings and contain no UTF-8 BOM.
- `atacseq2tracks.sh` and `fastq2tracks.sh` remain identical.
- `config/config.conf` and `config/config_temp.conf.template` remain identical.

The complete Bioconductor execution self-test and representative BAM analysis
must still be run in the server Conda environment:

```bash
conda activate ATACseq2tracks
cd /home/micgdu/Analysis/workflows/ATACseq2tracks

bash tests/check_bash_syntax.sh
Rscript scripts/deseq2atac_analysis.R --self-test
```

For an existing analysis, remove only the affected checkpoints before resuming:

```bash
rm -f /path/from_OUTPUT_DIR/.checkpoints/step12a.done \
      /path/from_OUTPUT_DIR/.checkpoints/step14.done
```

Legacy root-level DESeq2ATAC result files are not deleted or reused. New results
are written only to the `broad/` and `narrow/` subdirectories, and the wrapper
prints a warning if the preceding single-analysis layout is detected.
