# 18 — Technical replicates

## Meaning of the identifiers

- `sample_id` identifies a biological sample name.
- `replicate` identifies an independent biological replicate.
- `tech_replicate` identifies repeated sequencing of the same prepared library,
  such as different flow-cell lanes.

The biological-library key is `<sample_id>_bioR<replicate>`. Each samplesheet
row must have a unique `sample_id + replicate + tech_replicate` combination.
Rows belonging to one biological-library key must have identical biological
metadata; only FASTQ paths and `tech_replicate` may differ.

## Exact processing behavior

Technical FASTQs are concatenated in samplesheet order within each biological
library before Trim Galore. The combined data are then trimmed and aligned once.
They are not treated as independent observations in consensus support,
normalization or differential models. Biological replicates remain separate.

This preserves the existing v4.2.0 behavior; v4.3.0 adds auditability, not a new
merge algorithm. The persistent mapping is written to:

```text
metadata/validated_sequencing_units.tsv
metadata/biological_libraries.tsv
metadata/technical_merge_audit.tsv
```

Review `technical_merge_audit.tsv` before a full run. It reports the downstream
library key, technical-replicate identifiers, number of units, ordered R1/R2
paths, merge stage and downstream observation count.

## Example

```csv
sample_id,fastq_1,fastq_2,layout,genome,assay,factor,condition,treatment,cell_type,replicate,tech_replicate,is_control,control_id,macs2_mode,blacklist,output_prefix
WT_1,/data/WT_1_L001_R1.fastq.gz,/data/WT_1_L001_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,WT,none,keratinocyte,1,1,FALSE,,both,/ref/hg38.blacklist.bed,WT_1
WT_1,/data/WT_1_L002_R1.fastq.gz,/data/WT_1_L002_R2.fastq.gz,PE,hg38,ATAC-seq,accessibility,WT,none,keratinocyte,1,2,FALSE,,both,/ref/hg38.blacklist.bed,WT_1
```

These rows produce one trimmed/aligned biological library, `WT_1_bioR1`, not
two replicates.

## Limitations and design guidance

- Do not label biological replicates as technical replicates; that would remove
  independent observations from dispersion estimation.
- Do not combine lanes from different library preparations under one biological
  key unless they truly represent the same prepared library.
- Lane-specific bias can be reviewed in raw FastQC, but downstream QC is for the
  combined library. Investigate a clearly failed lane before pooling it.
- Technical replication does not compensate for insufficient biological
  replication. Differential conditions still need at least two biological
  libraries to enter models.
