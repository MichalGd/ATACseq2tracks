#!/usr/bin/env Rscript
# ATACseq2tracks v3.1.x — Consensus peak count matrix and DESeq2 size factors
# Usage: Rscript scripts/consensus_peak_size_factors.R <samplesheet.csv> <multiBamSummary_peaks.tab> <out_table_dir> <out_counts.tsv> <out_norm_counts.tsv>

suppressPackageStartupMessages({
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: consensus_peak_size_factors.R <samplesheet.csv> <multiBamSummary_peaks.tab> <out_table_dir> <out_counts.tsv> <out_norm_counts.tsv>")
}

ss_file <- args[1]
counts_file <- args[2]
out_table_dir <- args[3]
out_counts_file <- args[4]
out_norm_counts_file <- args[5]

if (!file.exists(ss_file)) stop("Samplesheet not found: ", ss_file)
if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!dir.exists(out_table_dir)) dir.create(out_table_dir, recursive = TRUE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE)
if (!all(c("sample_id", "rep", "is_control", "condition") %in% names(ss))) {
  stop("Samplesheet must contain sample_id, rep, is_control, and condition columns")
}

ss$sample_id <- trimws(ss$sample_id)
ss$rep <- trimws(as.character(ss$rep))
ss$is_control <- tolower(trimws(as.character(ss$is_control)))
ss$key <- paste0(ss$sample_id, "_bioR", ss$rep)
ss <- ss[!ss$is_control %in% c("true", "1", "yes"), ]
ss <- ss[!duplicated(ss$key), ]

counts_df <- read.delim(counts_file, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(counts_df) < 4) stop("Counts file must contain region columns followed by sample count columns")

region_cols <- c("chr", "chrom", "start", "end", "name", "region")
region_cols <- intersect(region_cols, tolower(names(counts_df)))
if (length(region_cols) < 3) {
  stop("Counts file must include at least chr, start and end columns")
}

# Find region columns by position, not name, to be robust to deepTools output column names.
chr_col <- 1
start_col <- 2
end_col <- 3
if (!all(c("chr", "start", "end") %in% tolower(names(counts_df)[1:3]))) {
  warning("Assuming first three columns are chr,start,end")
}
region_names <- paste0(counts_df[[chr_col]], ":", counts_df[[start_col]], "-", counts_df[[end_col]])

count_cols <- names(counts_df)[-(1:3)]
if (length(count_cols) < 1) stop("No sample count columns found in counts file")
counts <- counts_df[, count_cols, drop = FALSE]
rownames(counts) <- region_names

sample_keys <- ss$key
missing_keys <- setdiff(sample_keys, colnames(counts))
if (length(missing_keys) > 0) {
  warning("The following expected sample keys are missing from the count matrix: ", paste(missing_keys, collapse = ", "))
}
counts <- counts[, intersect(colnames(counts), sample_keys), drop = FALSE]
if (ncol(counts) == 0) stop("No matching sample columns found between the samplesheet and the counts matrix")

meta <- ss[match(colnames(counts), ss$key), c("sample_id", "key", "condition", "treatment", "cell_type", "assay", "genome")]
rownames(meta) <- meta$key

# Build DESeq2 object and estimate size factors
colData <- DataFrame(
  condition = factor(meta$condition),
  treatment = if ("treatment" %in% names(meta)) factor(meta$treatment) else NULL,
  cell_type = if ("cell_type" %in% names(meta)) factor(meta$cell_type) else NULL,
  assay = if ("assay" %in% names(meta)) factor(meta$assay) else NULL,
  row.names = meta$key
)

design_formula <- if (length(unique(meta$condition)) > 1) ~ condition else ~ 1

# ensure count matrix is integer
counts_matrix <- as.matrix(counts)
storage.mode(counts_matrix) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = colData, design = design_formula)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)

sf_table <- data.frame(
  sample_id = meta$sample_id,
  key = meta$key,
  size_factor = size_factors,
  stringsAsFactors = FALSE
)

size_factor_file <- file.path(out_table_dir, "consensus_sizeFactors.tsv")
write.table(sf_table, file = size_factor_file, sep = "\t", row.names = FALSE, quote = FALSE)

counts_out <- cbind(region = rownames(counts), counts)
write.table(counts_out, file = out_counts_file, sep = "\t", row.names = FALSE, quote = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_out <- cbind(region = rownames(norm_counts), norm_counts)
write.table(norm_counts_out, file = out_norm_counts_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("Written:", size_factor_file)
message("Written:", out_counts_file)
message("Written:", out_norm_counts_file)
