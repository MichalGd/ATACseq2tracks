#!/usr/bin/env Rscript
# ATACseq2tracks v3.2.0 - DESeq2 size factors from consensus-peak fragment/read counts
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: consensus_peak_size_factors.R <samplesheet.csv> <multiBamSummary.tab> <table_dir> <raw_counts.tsv> <normalized_counts.tsv>")
}
ss_file <- args[[1]]; count_file <- args[[2]]; table_dir <- args[[3]]
raw_out <- args[[4]]; normalized_out <- args[[5]]
if (!file.exists(ss_file)) stop("Samplesheet not found: ", ss_file)
if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_id", "replicate", "is_control", "genome")
missing <- setdiff(required, names(ss))
if (length(missing)) stop("Samplesheet is missing: ", paste(missing, collapse = ", "))
ss$sample_id <- trimws(ss$sample_id)
ss$replicate <- trimws(as.character(ss$replicate))
ss$is_control <- tolower(trimws(as.character(ss$is_control)))
ss$key <- paste0(ss$sample_id, "_bioR", ss$replicate)
ss <- ss[!ss$is_control %in% c("true", "1", "yes"), , drop = FALSE]
ss <- ss[!duplicated(ss$key), , drop = FALSE]
if (length(unique(ss$genome)) != 1L) stop("Size factors must be estimated separately for each genome build")

x <- read.delim(count_file, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
if (ncol(x) < 4L) stop("Expected chr, start, end and at least one sample column")
names(x) <- gsub("^#|'|\"", "", names(x))
region <- paste0(x[[1]], ":", x[[2]], "-", x[[3]])
sample_names <- gsub("^#|'|\"", "", names(x)[-(1:3)])
count_df <- x[, -(1:3), drop = FALSE]
names(count_df) <- sample_names

missing_samples <- setdiff(ss$key, names(count_df))
if (length(missing_samples)) stop("Samples missing from peak-count matrix: ", paste(missing_samples, collapse = ", "))
count_df <- count_df[, ss$key, drop = FALSE]
count_df[] <- lapply(count_df, function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (anyNA(v) || any(v < 0)) stop("Count matrix contains missing or negative values")
  as.integer(round(v))
})
rownames(count_df) <- region
count_df <- count_df[rowSums(count_df) > 0, , drop = FALSE]
if (!nrow(count_df)) stop("Every consensus peak has zero counts")

meta <- ss[match(names(count_df), ss$key), , drop = FALSE]
rownames(meta) <- meta$key
dds <- DESeqDataSetFromMatrix(as.matrix(count_df), DataFrame(row.names = meta$key), design = ~1)
# poscounts is robust to sparse ATAC peak matrices with many zeros.
dds <- estimateSizeFactors(dds, type = "poscounts")
sf <- sizeFactors(dds)
if (any(!is.finite(sf)) || any(sf <= 0)) stop("DESeq2 returned invalid size factors")

# Match DESeq2::fpm(dds, robust=TRUE) exactly. The robust effective library
# size is the DESeq2 size factor multiplied by the cohort geometric mean of
# the raw count-matrix column sums.
consensus_count_sum <- colSums(counts(dds))
if (any(!is.finite(consensus_count_sum)) || any(consensus_count_sum <= 0)) {
  stop("Every sample must have a positive consensus-peak count sum")
}
cohort_geometric_mean_column_sum <- exp(mean(log(consensus_count_sum)))
robust_effective_library_size <- sf * cohort_geometric_mean_column_sum
robust_cpm_scale <- 1e6 / robust_effective_library_size

# Guard the exported genome-wide scale against drift from DESeq2's definition.
robust_fpm <- fpm(dds, robust = TRUE)
expected_robust_fpm <- sweep(as.matrix(counts(dds)), 2, robust_cpm_scale, "*")
if (!isTRUE(all.equal(unname(robust_fpm), unname(expected_robust_fpm),
                      tolerance = 1e-10, check.attributes = FALSE))) {
  stop("Internal robust CPM scale does not reproduce DESeq2::fpm(robust=TRUE)")
}

sf_table <- data.frame(
  sample_id = meta$sample_id,
  key = meta$key,
  genome = meta$genome,
  size_factor = as.numeric(sf),
  track_scale_factor = as.numeric(1 / sf),
  deseq2_consensus_scale = as.numeric(1 / sf),
  consensus_count_sum = as.numeric(consensus_count_sum),
  cohort_geometric_mean_column_sum = rep(cohort_geometric_mean_column_sum, length(sf)),
  robust_effective_library_size = as.numeric(robust_effective_library_size),
  deseq2_robust_cpm_scale = as.numeric(robust_cpm_scale),
  basis = "DESeq2_poscounts_on_consensus_peaks",
  stringsAsFactors = FALSE
)
write.table(sf_table, file.path(table_dir, "consensus_sizeFactors.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(region = rownames(count_df), count_df), raw_out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(region = rownames(counts(dds, normalized = TRUE)), counts(dds, normalized = TRUE)),
            normalized_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("DESeq2 consensus and robust CPM scales written for ", ncol(count_df),
        " samples and ", nrow(count_df), " peaks")
