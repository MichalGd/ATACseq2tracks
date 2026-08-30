# 19 — Validation matrix and reference provenance

## Required release validation

| Test | PE | SE | dm6 off | dm6 on | Technical lanes | Expected evidence |
|---|---:|---:|---:|---:|---:|---|
| Static/focused tests | yes | yes | yes | mocked | yes | `tests/check_bash_syntax.sh` exits 0 |
| Plan | yes | yes | yes | yes | yes | metadata tables written; no analysis starts |
| Preflight | yes | yes | yes | yes | yes | zero failures and software/reference manifests |
| Small end-to-end run | yes | yes | yes | separate test | yes | all enabled non-empty outputs and report |
| Resume | yes | yes | yes | yes | yes | matching stages skip; event log records skips |
| Named boundary | yes | yes | yes | yes | yes | `--stop-after` and `--from-stage` behave as documented |
| Failure injection | yes | yes | no preference | no preference | no preference | stage exits non-zero and lacks completion checkpoint |
| Cleanup | yes | yes | yes | yes | yes | only configured intermediates removed and manifest complete |
| Report-only regeneration | yes | yes | yes | yes | yes | report changes without analysis-stage execution |

Scientific acceptance additionally requires the inherited v4.2.0 checks for BAM
integrity, PE/SE signal units, canonical chromosomes, normalization formulas,
peak support, QC tables, differential contrasts, annotations and spike-in QC.

## Reference provenance

Preflight writes `metadata/reference_manifest.tsv` from the resolved
configuration. It records each configured index prefix, GTF, chromosome-size,
blacklist, cCRE, Picard and Kent-utility path plus presence status.
`metadata/software_versions.tsv` records resolved command paths and reported
versions.

The manifest answers which configured resources were visible to the run. It
does not infer the biological release of an ambiguously named file or hash
multi-gigabyte references. Administrators should therefore keep a separate
shared-reference manifest containing source URL/accession, genome assembly,
download date, transformation commands and SHA-256 checksums. Preserve
`CCRE_SOURCE_HG38`/`CCRE_SOURCE_MM39` in each project configuration.

At minimum verify that all resources use the same host assembly, chromosome
naming agrees across FASTA/index/sizes/blacklist/GTF/cCRE, and composite dm6
indices have their preparation manifests when spike-in processing is enabled.
