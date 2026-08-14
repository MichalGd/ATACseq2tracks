#!/usr/bin/env Rscript
# ATACseq2tracks v3.2.0 — explicit DESeq2 differential accessibility analysis
# Usage: Rscript scripts/diffbind_analysis.R <diffbind_samplesheet.csv> <out_dir> [summits]

suppressPackageStartupMessages({
    library(DiffBind)
    library(ggplot2)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
    stop("Usage: diffbind_analysis.R <diffbind_samplesheet.csv> <out_dir> [summits]")
}

ss_file <- args[1]
out_dir  <- args[2]
summits  <- if (length(args) == 3) suppressWarnings(as.integer(args[3])) else 100L
if (is.na(summits) || summits < 0L) {
    stop("summits must be a non-negative integer; received: ", ifelse(length(args) == 3, args[3], "NA"))
}
if (!file.exists(ss_file)) stop("Sample sheet not found: ", ss_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE)
if (nrow(ss) == 0) stop("Empty DiffBind sample sheet: ", ss_file)

sample_type <- tools::file_path_sans_ext(basename(ss_file))
summary_file <- file.path(out_dir, "diffbind_summary.txt")
log_file     <- file.path(out_dir, "diffbind_log.txt")
results_csv  <- file.path(out_dir, "diffbind_results.csv")
all_results_csv <- file.path(out_dir, "diffbind_results_all_sites.csv")
consensus_bed <- file.path(out_dir, "diffbind_consensus_peaks.bed")

sink(log_file, append = FALSE)
cat("DiffBind analysis for", ss_file, "\n")
cat("Loaded", nrow(ss), "samples\n")
cat("ATAC summit half-width:", summits, "bp (window width", 2L * summits + 1L, "bp)\n")

if (!"Condition" %in% colnames(ss)) {
    stop("DiffBind sample sheet must include a Condition column")
}

ss$Condition <- trimws(as.character(ss$Condition))
if (any(is.na(ss$Condition) | ss$Condition == "")) {
    stop("Condition values must be non-empty")
}
condition_counts <- table(ss$Condition)
if (length(condition_counts) != 2) {
    stop("This v3.2.0 test workflow requires exactly two conditions; found: ",
         paste(names(condition_counts), collapse = ", "))
}
if (any(condition_counts < 2)) {
    stop("Each condition requires at least two biological replicates: ",
         paste(names(condition_counts), condition_counts, sep = "=", collapse = ", "))
}

# Create the DiffBind object
cat("Building DiffBind object...\n")
dba_obj <- dba(
    sampleSheet = ss_file,
    config = list(
        AnalysisMethod = DBA_DESEQ2,
        doBlacklist = FALSE,
        doGreylist = FALSE
    )
)
cat("Samples loaded:", nrow(ss), "\n")

# Count reads over peaks using SummarizedOverlaps
cat("Counting reads over peak regions...\n")
dba_obj <- dba.count(
    dba_obj,
    summits = summits,
    bUseSummarizeOverlaps = TRUE
)

# Normalise
cat("Normalising with DiffBind defaults...\n")
dba_obj <- dba.normalize(dba_obj)

# Contrast by condition; use minimum 2 members per condition
cat("Creating contrasts by condition...\n")
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)

# Run differential analysis using the default method (DESeq2)
cat("Running differential analysis...\n")
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# Export both the complete quantitative result and the FDR-filtered subset.
cat("Extracting differential results...\n")
res <- dba.report(dba_obj, contrast = 1, method = DBA_DESEQ2, th = 0.05)
res_all <- dba.report(dba_obj, contrast = 1, method = DBA_DESEQ2, th = 1)

# Write results and optional consensus peak set
cat("Writing results...\n")
write.csv(as.data.frame(res), results_csv, row.names = FALSE)
write.csv(as.data.frame(res_all), all_results_csv, row.names = FALSE)
if (length(res) > 0) {
    rtracklayer::export(res, consensus_bed, format = "BED")
}
saveRDS(dba_obj, file.path(out_dir, "diffbind_analysis_object.rds"))

# Plot diagnostics
plot_png <- function(filename, code) {
    png(filename, width = 1400, height = 1000, res = 150)
    on.exit(dev.off(), add = TRUE)
    code
}

cat("Writing PCA plot...\n")
plot_png(file.path(out_dir, "diffbind_pca.png"), {
    # th=1 includes all normalized consensus sites, so this QC plot remains
    # informative even when no individual site passes the FDR threshold.
    print(dba.plotPCA(dba_obj, contrast = 1, method = DBA_DESEQ2, th = 1,
                      attributes = DBA_CONDITION, label = DBA_ID))
})

cat("Writing heatmap plot...\n")
plot_png(file.path(out_dir, "diffbind_heatmap.png"), {
    dba.plotHeatmap(dba_obj, contrast = 1, method = DBA_DESEQ2, th = 1,
                    correlations = FALSE)
})

cat("Writing MA plot...\n")
plot_png(file.path(out_dir, "diffbind_ma.png"), {
    dba.plotMA(dba_obj, contrast = 1, method = DBA_DESEQ2)
})

cat("Writing volcano plot...\n")
plot_png(file.path(out_dir, "diffbind_volcano.png"), {
    print(dba.plotVolcano(dba_obj, contrast = 1, method = DBA_DESEQ2,
                          bReturnSites = FALSE))
})

# Summary report
cat("Writing summary file...\n")
sink()
summary_text <- c(
    paste("Sample sheet:", ss_file),
    paste("Output directory:", normalizePath(out_dir)),
    paste("Samples:", nrow(ss)),
    paste("Conditions:", paste(unique(ss$Condition), collapse = ", ")),
    paste("Summit half-width (bp):", summits),
    paste("Counting window width (bp):", 2L * summits + 1L),
    paste("All tested sites:", length(res_all)),
    paste("FDR <= 0.05 sites:", length(res)),
    paste("All-sites results:", all_results_csv),
    paste("FDR-filtered results:", results_csv),
    paste("PCA plot:", file.path(out_dir, "diffbind_pca.png")),
    paste("Heatmap plot:", file.path(out_dir, "diffbind_heatmap.png")),
    paste("MA plot:", file.path(out_dir, "diffbind_ma.png")),
    paste("Volcano plot:", file.path(out_dir, "diffbind_volcano.png")),
    paste("Consensus BED:", ifelse(file.exists(consensus_bed), consensus_bed, "none"))
)
writeLines(summary_text, summary_file)

cat("DiffBind analysis complete.\n")
