# Apply the ATACseq2tracks v4.0.0 reporting and DiffBind hotfix

From the parent of the installed repository:

```bash
cd /home/micgdu/Analysis/workflows
cp -a ATACseq2tracks ATACseq2tracks_before_reporting_diffbind_hotfix_20260820
tar -xzf ATACseq2tracks_v4.0.0_reporting_diffbind_hotfix_2026-08-20_overlay.tar.gz
cd ATACseq2tracks
find . -type f -name '*.sh' -exec chmod u+x {} +
bash tests/check_bash_syntax.sh
```

For future analyses, launch the workflow normally. The existing completed
differential analysis does not need to be repeated.

To regenerate only the corrected unified report for the completed run:

```bash
export F2T_CONFIG=/home/micgdu/Analysis/run7_ATAC/run7_ATACb_hs/config/config.conf
OUT=/home/micgdu/Analysis/run7_ATAC/run7_ATACb_hs/config/atacseq2tracks
bash /home/micgdu/Analysis/workflows/ATACseq2tracks/scripts/generate_pipeline_report.sh \
  "$OUT" "$OUT/reports/pipeline_report_$(date +%Y%m%d)" html
```

This direct report command does not rerun alignments, peak calling or
differential models and does not alter workflow checkpoints.
