# Minimal critical improvements for ATACseq2tracks 3.2.0

This update intentionally retains the Bash architecture. It is a focused safety and scientific-correctness update, not a Nextflow or Snakemake migration.

## Immediate priorities implemented in this overlay

1. **Stop silent success**

   FastQC, trimming, alignment, Picard, filtering, coverage and MACS3 batch jobs now propagate failures. Empty MACS3 peak outputs fail by default. Checkpoints contain a signature of the sample sheet, configuration and version instead of being empty marker files.

2. **Correct raw and trimmed FastQC selection**

   Raw QC explicitly uses sample-sheet FASTQs. Trimmed QC explicitly scans the trimmed FASTQ directory, preventing the raw files from being analysed twice.

3. **Apply a minimum usable-read filter**

   The filtering stage applies a configurable MAPQ threshold, excludes unmapped, secondary, supplementary, QC-failed and duplicate records, removes mitochondrial contigs, removes blacklist overlaps, and restores proper-pair integrity for PE data. Per-sample attrition metrics are written alongside the filtered BAM.

4. **Use paired-end fragments for MACS3**

   PE samples use `-f BAMPE`. SE samples use configurable ATAC shift and extension values. The sample-sheet `macs2_mode` is honoured instead of always running both modes. The executable is consistently `macs3`; legacy script and column names remain for compatibility.

5. **Generate browser bigWigs by default**

   Two visually comparable but scientifically distinct track families are produced:

   - `bigwig/*_RPM.bw`: deepTools CPM/RPM tracks generated from filtered BAMs.
   - `bigwig_deseq2_consensus/*_DESeq2Consensus.bw`: tracks scaled by the inverse DESeq2 size factor estimated from reads in consensus peaks.

   The DESeq2-consensus tracks are not additionally divided by total reads. They are deliberately not called RPM or RPKM. Neither track family is used as the quantitative input to differential analysis.

6. **Make consensus support explicit**

   The default is `CONSENSUS_MIN_SAMPLES=2`. “Sample” means a distinct biological sample key (`sample_id + replicate`) after technical-replicate rows are collapsed. Support is evaluated within one genome build and one selected peak type, narrow by default. There is no silent fallback to a single-sample union.

7. **Provide UCSC custom-track definitions**

   Each bigWig directory receives `ucsc_tracks.txt`. If `UCSC_BIGDATA_URL_BASE` is configured, the entries contain public URLs. If it is empty, relative filenames are emitted so the bigWigs and track file can be deployed together.

8. **Remove high-impact inconsistencies**

   The runtime version comes from `VERSION`; preflight checks MACS3 and DESeq2 rather than retired ChIPQC dependencies; the environment declares DESeq2 and removes ChIPQC; example paths are placeholders; and a bulk ATAC-seq example sample sheet is included.

9. **Add optional interpretation hooks**

   `peak_interpretation.sh` supports HOMER peak annotation and motif enrichment. Both are disabled by default and require explicit configuration and installed HOMER commands.

10. **Add ATAC-specific TSS and nucleosome-periodicity QC**

   `ataqv` now runs by default for each non-control ATAC-seq biological sample. A configured strand-aware BED6 TSS reference is used when supplied; otherwise it is generated deterministically from the genome GTF. Outputs include ENCODE-style TSS enrichment, the ataqv short-to-mononucleosomal ratio, compressed JSON metrics, and an optional local interactive viewer.

   Paired-end libraries additionally receive full-scan fragment-length histograms, 300-dpi PNG and vector PDF plots, and a compact metrics table. The table reports nucleosome-free (1–100 bp), mono- (180–247 bp), di- (315–473 bp), and tri-nucleosome (558–615 bp) fractions, NFR/mono ratio, periodic-fragment fraction, local mono/di peaks, and their spacing. Single-end libraries report periodicity as not applicable because insert sizes cannot be inferred reliably.

## Can consensus mean a peak present in at least two samples?

Yes. That is a reasonable minimal default when the “samples” are independent biological samples rather than technical lanes or duplicate peak files. It is best understood as a support-filtered consensus, not formal reproducibility testing.

The default retains intervals overlapped by peaks in at least two biological samples. With two or more biological replicates per condition, this can retain condition-specific peaks supported by two replicates. With only one replicate per condition, condition-specific singleton peaks can be excluded and valid differential inference is not available anyway.

This rule does not replace IDR, pooled-pseudoreplicate analysis or a more detailed replicate-aware consensus strategy. Those remain later enhancements.

## Important remaining limitations

- The workflow is still Bash and uses relatively simple CSV parsing.
- Technical FASTQs are still concatenated before trimming and alignment; raw FastQC remains per FASTQ, but per-lane alignment metrics are not retained.
- NRF/PBC/preseq and formal IDR are not added by this minimal update.
- The extra ataqv alignment/duplication/mitochondrial fields describe the already filtered analysis BAM; the pipeline's dedicated pre-filter attrition and Picard tables remain authoritative for those metrics.
- Differential modelling remains the existing DiffBind layer and still needs a later design/contrast overhaul.
- Pooled control BAMs are not inferred; standard bulk ATAC-seq is expected to run without an input control.
- DESeq2 size factors based on consensus peaks assume most counted accessible regions do not undergo a coordinated global shift. Experiments expecting global accessibility changes need spike-in or another explicit calibration strategy.
- Motif enrichment and annotation depend on separately installed and configured HOMER resources.

## Recommended review order

1. `config/config.conf`
2. `scripts/blacklist_filter.sh`
3. `scripts/macs2_peaks.sh` and `scripts/macs2_batch.sh`
4. `scripts/post_alignment_qc_batch.sh`
5. `scripts/consensus_peak_size_factors.R`
6. `scripts/genomecoverage_single.sh`
7. `atacseq2tracks.sh`
