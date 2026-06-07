#!/usr/bin/env Rscript
# ATACseq2tracks v3.1.x — DiffBind differential accessibility analysis
# Usage: Rscript scripts/diffbind_analysis.R <diffbind_samplesheet.csv> <out_dir>

suppressPackageStartupMessages({
    library(DiffBind)
    library(ggplot2)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: diffbind_analysis.R <diffbind_samplesheet.csv> <out_dir>")
}

ss_file <- args[1]
out_dir  <- args[2]
if (!file.exists(ss_file)) stop("Sample sheet not found: ", ss_file)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE)
if (nrow(ss) == 0) stop("Empty DiffBind sample sheet: ", ss_file)

sample_type <- tools::file_path_sans_ext(basename(ss_file))
summary_file <- file.path(out_dir, "diffbind_summary.txt")
log_file     <- file.path(out_dir, "diffbind_log.txt")
results_csv  <- file.path(out_dir, "diffbind_results.csv")
consensus_bed <- file.path(out_dir, "diffbind_consensus_peaks.bed")

sink(log_file, append = FALSE)
cat("DiffBind analysis for", ss_file, "\n")
cat("Loaded", nrow(ss), "samples\n")

if (!"Condition" %in% colnames(ss)) {
    stop("DiffBind sample sheet must include a Condition column")
}

if (length(unique(ss$Condition)) < 2) {
    stop("Need at least two distinct conditions for differential analysis")
}

# Create the DiffBind object
cat("Building DiffBind object...\n")
dba_obj <- dba(sampleSheet = ss_file)
cat("Samples loaded:", dba_obj$config$sampleCount, "\n")

# Count reads over peaks using SummarizedOverlaps
cat("Counting reads over peak regions...\n")
dba_obj <- dba.count(dba_obj, bUseSummarizeOverlaps = TRUE)
cat("Count matrix dimensions:", dim(dba_obj$counts), "\n")

# Normalise
cat("Normalising with DiffBind defaults...\n")
dba_obj <- dba.normalize(dba_obj)

# Contrast by condition; use minimum 2 members per condition
cat("Creating contrasts by condition...\n")
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)

# Run differential analysis using the default method (DESeq2)
cat("Running differential analysis...\n")
dba_obj <- dba.analyze(dba_obj)

# Extract results at FDR 0.05
cat("Extracting differential results...\n")
res <- dba.report(dba_obj, th = 0.05)

# Write results and optional consensus peak set
cat("Writing results...\n")
write.csv(as.data.frame(res), results_csv, row.names = FALSE)
if (length(res) > 0) {
    export.bed(res, consensus_bed)
}

# Plot diagnostics
plot_png <- function(filename, code) {
    png(filename, width = 1400, height = 1000, res = 150)
    code
    dev.off()
}

cat("Writing PCA plot...\n")
plot_png(file.path(out_dir, "diffbind_pca.png"), {
    dba.plotPCA(dba_obj, attributes = DBA_CONDITION, label = DBA_ID)
})

cat("Writing heatmap plot...\n")
plot_png(file.path(out_dir, "diffbind_heatmap.png"), {
    dba.plotHeatmap(dba_obj)
})

cat("Writing MA plot...\n")
plot_png(file.path(out_dir, "diffbind_ma.png"), {
    dba.plotMA(dba_obj)
})

cat("Writing volcano plot...\n")
plot_png(file.path(out_dir, "diffbind_volcano.png"), {
    dba.plotVolcano(dba_obj)
})

# Summary report
cat("Writing summary file...\n")
sink()
summary_text <- c(
    paste("Sample sheet:", ss_file),
    paste("Output directory:", normalizePath(out_dir)),
    paste("Samples:", nrow(ss)),
    paste("Conditions:", paste(unique(ss$Condition), collapse = ", ")),
    paste("Results rows:", ifelse(file.exists(results_csv), nrow(read.csv(results_csv)), 0)),
    paste("PCA plot:", file.path(out_dir, "diffbind_pca.png")),
    paste("Heatmap plot:", file.path(out_dir, "diffbind_heatmap.png")),
    paste("MA plot:", file.path(out_dir, "diffbind_ma.png")),
    paste("Volcano plot:", file.path(out_dir, "diffbind_volcano.png")),
    paste("Consensus BED:", ifelse(file.exists(consensus_bed), consensus_bed, "none"))
)
writeLines(summary_text, summary_file)

cat("DiffBind analysis complete.\n")
