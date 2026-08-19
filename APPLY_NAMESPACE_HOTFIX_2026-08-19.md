# Apply the ATACseq2tracks 3.2.0 namespace hotfix

This overlay corrects only the fatal IRanges namespace calls and adds their
preflight/regression coverage. Apply it to the same active workflow directory
used by the failed run.

## 1. Back up the active workflow

```bash
cd /home/micgdu/Analysis/workflows
cp -a ATACseq2tracks \
  "ATACseq2tracks_before_iranges_hotfix_$(date +%Y%m%d_%H%M%S)"
```

## 2. Extract the overlay

Place `ATACseq2tracks_v3.2.0_namespace_hotfix_2026-08-19_overlay.tar.gz`
in `/home/micgdu/Analysis/workflows`, then run:

```bash
cd /home/micgdu/Analysis/workflows
tar -xzf ATACseq2tracks_v3.2.0_namespace_hotfix_2026-08-19_overlay.tar.gz
```

The archive is rooted at `ATACseq2tracks/` and overwrites only the listed
hotfix files. It does not contain a replacement run config or samplesheet.

## 3. Validate the installation

```bash
cd /home/micgdu/Analysis/workflows/ATACseq2tracks
bash tests/check_bash_syntax.sh

Rscript -e '
stopifnot(requireNamespace("IRanges", quietly=TRUE))
exports <- getNamespaceExports("IRanges")
stopifnot(all(c("overlapsAny", "findOverlaps") %in% exports))
cat("IRanges namespace checks: OK\n")
'
```

The focused suite should report:

```text
OK   IRanges overlap generics dispatch correctly for GRanges
```

## 4. Preserve and clear only failed downstream results

After the failed workflow has exited:

```bash
OUT=/home/micgdu/Analysis/run7_ATAC/run7_ATACb_hs/config/atacseq2tracks
STAMP=$(date +%Y%m%d_%H%M%S)

[[ -d "$OUT/diffbind_results" ]] && \
  mv "$OUT/diffbind_results" "$OUT/diffbind_results_failed_namespace_$STAMP"
[[ -d "$OUT/deseq2atac" ]] && \
  mv "$OUT/deseq2atac" "$OUT/deseq2atac_failed_namespace_$STAMP"

rm -f "$OUT/.checkpoints/step12.done" \
      "$OUT/.checkpoints/step12a.done" \
      "$OUT/.checkpoints/step14.done"
```

Do not remove Steps 1--11 or 13. Do not delete filtered BAMs, peak calls, QC or
track outputs.

## 5. Resume

Activate the same environment and run the same command with the unchanged run
config. Because the active workflow path, config, samplesheet and `VERSION`
remain unchanged, compatible upstream checkpoints are reused. Differential
Steps 12/12a and report Step 14 rerun.

Review the Day4 replicate 1 and Day7 replicate 1 QC outliers before interpreting
the corrected statistical results. This hotfix deliberately does not make an
automatic biological sample-exclusion decision.
