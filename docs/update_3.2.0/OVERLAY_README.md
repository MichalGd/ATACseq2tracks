# ATACseq2tracks critical Bash update overlay

> **Historical initial overlay (2026-08-02):** This page predates the CPM,
> robust-CPM, universal multi-condition and built-in annotation updates. For
> current behavior use the [main README](../../README.md) and
> [documentation index](../README.md).

This reviewable overlay contains the proposed minimal urgent update for the current ATACseq2tracks repository. It intentionally does not introduce Nextflow or Snakemake and has not been applied to the active workflow directories.

Start with:

- [`MINIMAL_CRITICAL_IMPROVEMENTS.md`](MINIMAL_CRITICAL_IMPROVEMENTS.md) for scope, scientific choices and limitations;
- [`APPLY_UPDATE.md`](APPLY_UPDATE.md) for review, application and acceptance checks;
- [`CHANGELOG.md`](../../CHANGELOG.md) for the proposed release notes.

The principal default outputs added or corrected are:

- `bigwig/*_RPM.bw` for UCSC/IGV visualization;
- `bigwig_deseq2_consensus/*_DESeq2Consensus.bw` using DESeq2 size factors estimated from consensus-peak reads;
- `qc_post_alignment/peak_sets/consensus_peaks.bed`, requiring support in at least two distinct biological samples by default;
- `qc_post_alignment/atac_qc/` with ataqv TSS enrichment, compressed metrics, an optional interactive viewer, and paired-end nucleosome-periodicity metrics/plots;
- `ucsc_tracks.txt` in each bigWig directory.

Motif enrichment and peak annotation are optional and disabled by default.
