ATACseq2tracks v4.3.2 configuration directory

For the shared installation, copy config.conf into a user-owned project folder,
set SAMPLESHEET and OUTPUT_DIR plus server references, then run:

  atacseq2tracks --config /absolute/path/project/config/config.conf --plan
  atacseq2tracks --config /absolute/path/project/config/config.conf --preflight-only
  atacseq2tracks --config /absolute/path/project/config/config.conf

The configuration is parsed as literal data. Do not use shell variables,
command substitutions, pipes or custom undocumented keys.
================================================

Files
-----
config.conf                 Main safe configuration template.
config_temp.conf.template   Identical compatibility template.
samplesheet_template.csv    Column header template.
samplesheet_example_atac.csv
                            Recommended bulk ATAC-seq example.
samplesheet_example_atac_spikein.csv
                            Optional dm6 spike-in declaration example.
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
8. All five coverage families and both bigWig/bedGraph formats are enabled by
   default. Set an individual GENERATE_*_TRACKS switch to false to disable that
   family. Disabling permissive/intermediate also skips its BAM filtering.
9. dm6 calibration remains off by default. To enable it, use the spike-in
   example sheet, configure INDEX_*_DM6/CHROM_SIZES_DM6/BLACKLIST_DM6 and set
   GENERATE_DROSOPHILA_SPIKEIN_STRINGENT_TRACKS=true. Raw and CPM dm6 UCSC
   control tracks are then produced by default; disable only those controls with
   GENERATE_DROSOPHILA_CONTROL_TRACKS=false.
10. Run the main entry point with:

   bash /path/to/ATACseq2tracks/atacseq2tracks.sh \
        --config /path/to/project/config/config.conf

Do not store project-specific paths or credentials in the shared repository.
