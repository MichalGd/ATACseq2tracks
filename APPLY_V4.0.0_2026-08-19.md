# Apply the ATACseq2tracks v4.0.0 release

Back up the installed repository, then overlay the archive from its parent
directory:

```bash
cd /home/micgdu/Analysis/workflows
cp -a ATACseq2tracks ATACseq2tracks_before_v4.0.0_20260819
tar -xzf ATACseq2tracks_v4.0.0_2026-08-19_overlay.tar.gz
cd ATACseq2tracks
find . -type f -name '*.sh' -exec chmod u+x {} +
bash tests/check_bash_syntax.sh
```

Existing run configs do not need editing. To override the defaults, add any of:

```bash
QC_SAMPLE_PARALLEL_JOBS=4
ATAQV_PARALLEL_JOBS=4
TRACK_PARALLEL_JOBS=2
POOLED_MACS_PARALLEL_JOBS=2
MERGE_PARALLEL_JOBS=2
```

Re-run the workflow with the same command. Existing valid checkpoints are
preserved. Remove only the checkpoint for a stage that you intentionally want
to rerun for performance validation.

## Publish the GitHub release

After copying the complete prepared repository into a review branch and
confirming the tests, commit it before creating the tag:

```bash
git switch -c release/v4.0.0
git add -A
git commit -m "ATACseq2tracks v4.0.0"
git push -u origin release/v4.0.0
```

Merge the reviewed branch into `main`, then tag the resulting `main` commit:

```bash
git switch main
git pull --ff-only origin main
git tag -a v4.0.0 -m "ATACseq2tracks v4.0.0"
git push origin main v4.0.0
```

Do not create the tag before the release commit is present on `main`.
