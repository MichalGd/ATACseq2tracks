ATACseq2tracks v3.2.0 configuration directory
================================================

Files
-----
config.conf                 Main safe configuration template.
config_temp.conf.template   Identical compatibility template.
samplesheet_template.csv    Column header template.
samplesheet_example_atac.csv
                            Recommended bulk ATAC-seq example.
samplesheet_example.csv     Historical multi-assay example retained from v3.1.

Setup
-----
1. Copy config.conf and a samplesheet outside the installed code directory.
2. Set SAMPLESHEET, OUTPUT_DIR, reference paths and software paths.
3. Configure one genome per run (hg38 or mm39).
4. Leave TSS_BED_* empty to derive a strand-aware TSS BED from the GTF.
5. Review cleanup, consensus, track and ATAC-QC defaults explicitly.
6. Run the main entry point with:

   bash /path/to/ATACseq2tracks/atacseq2tracks.sh \
        --config /path/to/project/config/config.conf

Do not store project-specific paths or credentials in the shared repository.
