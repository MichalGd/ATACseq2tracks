ATACseq2tracks v4.0.0 configuration directory
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
5. cCRE annotation is enabled by default. Verify the genome-matched
   CCRE_BED_* path. On a server without this reference, explicitly set
   RUN_CCRE_ANNOTATION=false for GTF-only annotation.
   Build/download the default human hg38 reference with:
   bash utilities/prepare_encode4_hg38_ccre.sh
6. Optionally set DIFFERENTIAL_CONDITION_ORDER; otherwise condition order is
   taken from first samplesheet appearance and all eligible pairs are analyzed.
7. Cleanup is enabled after full success. Set ENABLE_AUTOMATIC_CLEANUP=false
   before the run if every intermediate must be retained.
8. Run the main entry point with:

   bash /path/to/ATACseq2tracks/atacseq2tracks.sh \
        --config /path/to/project/config/config.conf

Do not store project-specific paths or credentials in the shared repository.
