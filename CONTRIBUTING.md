# Contributing

This repository is primarily a preservation and documentation package for the uploaded workflow scripts.

Suggested contribution style:

1. Keep original uploaded scripts unchanged unless a version bump is intentional.
2. Document any code change in `CHANGELOG.md`.
3. Run `bash tests/check_bash_syntax.sh` before committing.
4. Test on a tiny paired-end dataset before running on production data.
5. Do not commit generated FASTQ, BAM, bedGraph, BigWig, or run-output folders.

For future improvements, consider replacing hard-coded paths with a configuration file and making the workflow runnable from its repository root.
