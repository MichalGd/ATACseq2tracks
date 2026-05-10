# Changelog notes

## Version 2.1 packaging

This GitHub-ready package was rebuilt from the uploaded workflow files, including the newly added `generate_pipeline_report.1.0.sh` and the PowerPoint presentation `fastq2tracks_LAbMeeting_29jan2026.pptx`.

The original scripts were copied into `scripts/` without changing their contents.

## Historical notes from `readme.2.1.sh`

The uploaded `readme.2.1.sh` records changes from an earlier `fastq2tracks` 1.5 version to 2.1. The main points are:

- Picard Java heap was reduced from 128 GB to 32 GB.
- Batch scripts were adjusted to current lower-level script names.
- Mouse Bowtie2 mapping logic was changed to avoid earlier chromosome filtering before Picard.
- Bowtie2 thread usage was changed from 15 to 12.
- Samtools commands were changed to use multithreading.
- Human and mouse Bowtie2 lower-level scripts were separated by filename, although their current logic is very similar.

`readme.2.1.sh` remains in `scripts/` as an original note, but it is not executable shell code.
