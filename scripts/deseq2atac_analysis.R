#!/usr/bin/env Rscript
# ATACseq2tracks - independent broad/narrow consensus DESeq2 differential accessibility
#
# Usage:
#   Rscript deseq2atac_analysis.R <samplesheet.csv> <bam_dir> <peaks_dir> \
#     <out_dir> <genome> <blacklist.bed> <min_support> <alpha> \
#     <block_column_or_empty> <reference_condition_or_empty> <min_mapq> [peak_type]
#
# Paired-end counting deliberately uses singleEnd=FALSE, fragments=FALSE.
# In GenomicAlignments this reads proper mate pairs as GAlignmentPairs and drops
# singletons/invalid pairs. The input BAMs have already undergone proper-pair,
# duplicate, MAPQ, mitochondrial and blacklist filtering. Single-end input is
# counted with singleEnd=TRUE, one retained alignment per observation.

suppressPackageStartupMessages({
    library(DESeq2)
    library(GenomicAlignments)
    library(GenomicRanges)
    library(Rsamtools)
    library(rtracklayer)
    library(ggplot2)
})

trim_character_columns <- function(x) {
    for (name in names(x)) {
        if (is.character(x[[name]])) x[[name]] <- trimws(x[[name]])
    }
    x
}

is_control_value <- function(x) {
    tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

canonical_names <- function(genome, with_chr = TRUE) {
    maximum <- switch(tolower(genome), hg38 = 22L, mm39 = 19L,
                      stop("Unsupported genome: ", genome))
    names <- c(as.character(seq_len(maximum)), "X", "Y")
    if (with_chr) paste0("chr", names) else names
}

canonicalize_ranges <- function(gr, genome) {
    observed <- as.character(seqnames(gr))
    allowed <- c(canonical_names(genome, TRUE), canonical_names(genome, FALSE))
    gr[observed %in% allowed]
}

harmonize_seqname_style <- function(gr, reference_names) {
    if (!length(gr)) return(gr)
    current <- as.character(seqnames(gr))
    if (length(intersect(unique(current), unique(reference_names)))) return(gr)
    reference_has_chr <- any(grepl("^chr", reference_names))
    current_has_chr <- any(grepl("^chr", current))
    if (reference_has_chr && !current_has_chr) {
        seqlevels(gr) <- paste0("chr", seqlevels(gr))
    } else if (!reference_has_chr && current_has_chr) {
        seqlevels(gr) <- sub("^chr", "", seqlevels(gr))
    }
    gr
}

sort_canonical_ranges <- function(gr, genome) {
    if (!length(gr)) return(gr)
    with_chr <- any(grepl("^chr", as.character(seqnames(gr))))
    desired <- canonical_names(genome, with_chr)
    keep <- intersect(desired, seqlevels(gr))
    seqlevels(gr, pruning.mode = "coarse") <- keep
    sort(gr, ignore.strand = TRUE)
}

write_tsv_gz <- function(x, path) {
    connection <- gzfile(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.table(x, connection, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE, na = "NA")
}

write_matrix_gz <- function(matrix, regions, path) {
    table <- cbind(regions, as.data.frame(matrix, check.names = FALSE))
    write_tsv_gz(table, path)
}

significant_rows <- function(result_table, alpha) {
    !is.na(result_table$padj) & result_table$padj <= alpha
}

render_pair <- function(stem, plot_function, width = 7, height = 6) {
    png_arguments <- list(filename = paste0(stem, ".png"), width = width,
                          height = height, units = "in", res = 300)
    if (capabilities("cairo")) png_arguments$type <- "cairo"
    do.call(png, png_arguments)
    tryCatch(plot_function(), finally = dev.off())
    pdf(paste0(stem, ".pdf"), width = width, height = height,
        useDingbats = FALSE)
    tryCatch(plot_function(), finally = dev.off())
}

plot_matrix <- function(matrix, title, palette) {
    if (nrow(matrix) < 2L) {
        plot.new(); title(main = title); text(0.5, 0.5, "Insufficient samples")
        return(invisible(NULL))
    }
    heatmap(matrix, Rowv = NA, Colv = NA, scale = "none", symm = TRUE,
            col = palette, margins = c(10, 10), main = title)
}

run_synthetic_self_test <- function() {
    set.seed(320)
    n <- 400L
    baseline <- pmax(20L, rnbinom(n, mu = 120, size = 15))
    a1 <- pmax(1L, round(baseline * 0.93))
    a2 <- pmax(1L, round(baseline * 1.07))
    no_effect <- cbind(A1 = a1, A2 = a2, B1 = a1, B2 = a2)
    metadata <- data.frame(condition = factor(c("A", "A", "B", "B"), levels = c("A", "B")))
    rownames(metadata) <- colnames(no_effect)

    fit_one <- function(counts) {
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                      design = ~ condition)
        dds <- estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq(dds, fitType = "mean", quiet = TRUE)
        as.data.frame(results(dds, contrast = c("condition", "B", "A"), alpha = 0.05))
    }

    zero_result <- fit_one(no_effect)
    stopifnot(sum(significant_rows(zero_result, 0.05)) == 0L)

    strong_effect <- no_effect
    strong_effect[seq_len(80L), 3:4] <- strong_effect[seq_len(80L), 3:4] * 12L
    effect_result <- fit_one(strong_effect)
    stopifnot(sum(significant_rows(effect_result, 0.05)) > 0L)

    probe <- data.frame(padj = c(NA_real_, 0.01, 0.2))
    stopifnot(identical(significant_rows(probe, 0.05), c(FALSE, TRUE, FALSE)))

    peak_a <- GRanges("chr1", IRanges(c(1, 101), c(50, 150)))
    peak_b <- GRanges("chr1", IRanges(c(25, 201), c(75, 250)))
    peak_c <- GRanges("chr1", IRanges(c(40, 301), c(60, 350)))
    test_sets <- list(A = peak_a, B = peak_b, C = peak_c)
    test_atoms <- disjoin(unlist(GRangesList(test_sets), use.names = FALSE))
    test_support <- rowSums(do.call(cbind, lapply(test_sets, function(peaks) {
        countOverlaps(test_atoms, peaks) > 0L
    })))
    test_consensus <- reduce(test_atoms[test_support >= 2L])
    stopifnot(length(test_consensus) == 1L, start(test_consensus) == 25L,
              end(test_consensus) == 75L)

    compressed_probe <- tempfile(fileext = ".tsv.gz")
    write_tsv_gz(data.frame(id = 1:2, value = c("a", "b")), compressed_probe)
    stopifnot(identical(read.delim(gzfile(compressed_probe))$value, c("a", "b")))
    unlink(compressed_probe)

    seqinfo_probe <- Seqinfo(seqnames = "chr1", seqlengths = 1000L)
    first_mate <- GAlignments(
        seqnames = Rle(factor("chr1", levels = "chr1")), pos = 10L,
        cigar = "20M", strand = Rle(strand("+")), seqinfo = seqinfo_probe
    )
    second_mate <- GAlignments(
        seqnames = Rle(factor("chr1", levels = "chr1")), pos = 60L,
        cigar = "20M", strand = Rle(strand("-")), seqinfo = seqinfo_probe
    )
    pair_probe <- GAlignmentPairs(first_mate, second_mate, isProperPair = TRUE)
    feature_probe <- GRanges("chr1", IRanges(1L, 100L), seqinfo = seqinfo_probe)
    pair_count <- assay(summarizeOverlaps(
        feature_probe, pair_probe, mode = "Union", ignore.strand = TRUE,
        inter.feature = TRUE
    ))[[1]]
    read_count <- assay(summarizeOverlaps(
        feature_probe, c(first_mate, second_mate), mode = "Union",
        ignore.strand = TRUE, inter.feature = TRUE
    ))[[1]]
    stopifnot(pair_count == 1L, read_count == 2L)
    cat("OK   DESeq2ATAC synthetic significant and zero-significant analyses\n")
    cat("OK   DESeq2ATAC support threshold and compressed output self-tests\n")
    cat("OK   DESeq2ATAC paired fragment once and single-end read counting self-tests\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && identical(args[[1]], "--self-test")) {
    run_synthetic_self_test()
    quit(save = "no", status = 0L)
}
if (!(length(args) %in% c(11L, 12L))) {
    stop(paste(
        "Usage: deseq2atac_analysis.R <samplesheet.csv> <bam_dir> <peaks_dir>",
        "<out_dir> <genome> <blacklist.bed> <min_support> <alpha>",
        "<block_column_or_empty> <reference_condition_or_empty> <min_mapq> [broad|narrow]"
    ))
}

samplesheet <- args[[1]]
bam_dir <- args[[2]]
peaks_dir <- args[[3]]
out_dir <- args[[4]]
genome <- tolower(args[[5]])
blacklist_file <- args[[6]]
min_support <- suppressWarnings(as.integer(args[[7]]))
alpha <- suppressWarnings(as.numeric(args[[8]]))
block_column <- trimws(args[[9]])
reference_condition <- trimws(args[[10]])
min_mapq <- suppressWarnings(as.integer(args[[11]]))
peak_type <- if (length(args) == 12L) tolower(trimws(args[[12]])) else "broad"

if (!file.exists(samplesheet)) stop("Samplesheet not found: ", samplesheet)
if (!dir.exists(bam_dir)) stop("BAM directory not found: ", bam_dir)
if (!dir.exists(peaks_dir)) stop("Peak directory not found: ", peaks_dir)
if (!file.exists(blacklist_file)) stop("Blacklist not found: ", blacklist_file)
if (is.na(min_support) || min_support < 1L) stop("min_support must be a positive integer")
if (is.na(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between zero and one")
if (is.na(min_mapq) || min_mapq < 0L) stop("min_mapq must be a non-negative integer")
if (!genome %in% c("hg38", "mm39")) stop("Unsupported genome: ", genome)
if (!(peak_type %in% c("broad", "narrow"))) {
    stop("peak_type must be broad or narrow; received: ", peak_type)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("DESeq2ATAC: loading sample metadata")
ss <- trim_character_columns(read.csv(samplesheet, stringsAsFactors = FALSE,
                                      check.names = FALSE))
required <- c("sample_id", "layout", "genome", "condition", "replicate", "is_control")
missing_columns <- setdiff(required, names(ss))
if (length(missing_columns)) {
    stop("Samplesheet is missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (block_column != "" && !block_column %in% names(ss)) {
    stop("Configured DESeq2ATAC block column is absent from samplesheet: ", block_column)
}
if (block_column != "" && !grepl("^[A-Za-z][A-Za-z0-9_.]*$", block_column)) {
    stop("DESeq2ATAC block column must be a syntactically simple column name")
}

ss <- ss[tolower(ss$genome) == genome & !is_control_value(ss$is_control), , drop = FALSE]
if (!nrow(ss)) stop("No non-control samples found for genome ", genome)
ss$key <- paste0(ss$sample_id, "_bioR", ss$replicate)

consistency_columns <- c("condition", "layout", "genome")
if (block_column != "") consistency_columns <- c(consistency_columns, block_column)
for (key in unique(ss$key)) {
    rows <- ss[ss$key == key, , drop = FALSE]
    for (column in consistency_columns) {
        values <- unique(rows[[column]])
        values <- values[!is.na(values) & values != ""]
        if (length(values) != 1L) {
            stop("Inconsistent or missing ", column, " values for biological sample ", key)
        }
    }
}
metadata <- ss[!duplicated(ss$key), , drop = FALSE]
rownames(metadata) <- metadata$key
metadata$condition <- trimws(as.character(metadata$condition))
metadata$layout <- toupper(trimws(as.character(metadata$layout)))

if (length(unique(metadata$layout)) != 1L || !unique(metadata$layout) %in% c("PE", "SE")) {
    stop("DESeq2ATAC requires one PE-only or SE-only cohort")
}
layout <- unique(metadata$layout)
condition_counts <- table(metadata$condition)
if (length(condition_counts) != 2L) {
    stop("DESeq2ATAC currently requires exactly two conditions; found: ",
         paste(names(condition_counts), collapse = ", "))
}
if (any(condition_counts < 2L)) {
    stop("Each DESeq2ATAC condition requires at least two biological replicates: ",
         paste(names(condition_counts), condition_counts, sep = "=", collapse = ", "))
}
if (nrow(metadata) < min_support) {
    stop("DESeq2ATAC minimum peak support exceeds the number of biological samples")
}

conditions <- sort(names(condition_counts))
if (reference_condition == "") reference_condition <- conditions[[1]]
if (!reference_condition %in% conditions) {
    stop("Reference condition is not present: ", reference_condition)
}
numerator_condition <- setdiff(conditions, reference_condition)
if (length(numerator_condition) != 1L) stop("Unable to determine contrast numerator")

metadata$condition <- factor(metadata$condition, levels = c(reference_condition, numerator_condition))
if (block_column != "") metadata[[block_column]] <- factor(metadata[[block_column]])

metadata$bam <- file.path(bam_dir, paste0(metadata$key, "_dedup_blFilt.bam"))
peak_extension <- if (peak_type == "broad") "broadPeak" else "narrowPeak"
metadata$peak_file <- file.path(
    peaks_dir, "per_replicate", metadata$key, peak_type,
    paste0(metadata$key, "_peaks.", peak_extension)
)
missing_bams <- metadata$key[!file.exists(metadata$bam) | file.info(metadata$bam)$size <= 0]
missing_peaks <- metadata$key[!file.exists(metadata$peak_file) | file.info(metadata$peak_file)$size <= 0]
if (length(missing_bams)) stop("Missing/empty filtered BAMs: ", paste(missing_bams, collapse = ", "))
if (length(missing_peaks)) {
    stop("Missing/empty ", peak_type, " peaks (set macs2_mode=both): ",
         paste(missing_peaks, collapse = ", "))
}

message("DESeq2ATAC ", peak_type, ": importing and filtering peaks")
blacklist <- canonicalize_ranges(import(blacklist_file), genome)
first_peak <- canonicalize_ranges(import(metadata$peak_file[[1]]), genome)
if (!length(first_peak)) stop("First ", peak_type, " peak file has no canonical intervals")
blacklist <- harmonize_seqname_style(blacklist, as.character(seqnames(first_peak)))

peak_sets <- setNames(vector("list", nrow(metadata)), metadata$key)
for (index in seq_len(nrow(metadata))) {
    key <- metadata$key[[index]]
    peaks <- canonicalize_ranges(import(metadata$peak_file[[index]]), genome)
    peaks <- harmonize_seqname_style(peaks, as.character(seqnames(first_peak)))
    if (length(blacklist)) {
        peaks <- peaks[!overlapsAny(peaks, blacklist, ignore.strand = TRUE)]
    }
    strand(peaks) <- "*"
    peaks <- reduce(peaks, ignore.strand = TRUE)
    if (!length(peaks)) stop("No ", peak_type, " peaks remain after filtering for ", key)
    peak_sets[[key]] <- peaks
}

message("DESeq2ATAC ", peak_type, ": constructing replicate-supported consensus regions")
all_peaks <- unlist(GRangesList(peak_sets), use.names = FALSE)
atoms <- disjoin(all_peaks, ignore.strand = TRUE)
support_matrix <- do.call(cbind, lapply(peak_sets, function(peaks) {
    countOverlaps(atoms, peaks, ignore.strand = TRUE) > 0L
}))
colnames(support_matrix) <- names(peak_sets)
atom_support <- rowSums(support_matrix)
supported_atoms <- atoms[atom_support >= min_support]
if (!length(supported_atoms)) {
    stop("No ", peak_type, "-peak interval has support from at least ",
         min_support, " samples")
}
consensus <- reduce(supported_atoms, ignore.strand = TRUE)
consensus <- sort_canonical_ranges(consensus, genome)

consensus_support_matrix <- do.call(cbind, lapply(peak_sets, function(peaks) {
    countOverlaps(consensus, peaks, ignore.strand = TRUE) > 0L
}))
colnames(consensus_support_matrix) <- names(peak_sets)
consensus_support <- rowSums(consensus_support_matrix)
supporting_samples <- apply(consensus_support_matrix, 1L, function(values) {
    paste(colnames(consensus_support_matrix)[values], collapse = ",")
})
if (any(consensus_support < min_support)) stop("Internal consensus-support validation failed")

peak_ids <- sprintf("peak_%06d", seq_along(consensus))
names(consensus) <- peak_ids
region_table <- data.frame(
    peak_id = peak_ids,
    peak_type = peak_type,
    chrom = as.character(seqnames(consensus)),
    start = start(consensus) - 1L,
    end = end(consensus),
    width = width(consensus),
    sample_support = as.integer(consensus_support),
    supporting_samples = supporting_samples,
    stringsAsFactors = FALSE
)
bed_table <- region_table[, c("chrom", "start", "end", "peak_id")]
write.table(bed_table, file.path(out_dir, "deseq2atac_consensus_peaks.bed"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write_tsv_gz(region_table, file.path(out_dir, "deseq2atac_consensus_peaks_with_support.tsv.gz"))

message("DESeq2ATAC ", peak_type, ": counting ",
        ifelse(layout == "PE", "properly paired fragments", "single-end reads"))
flag_filter <- if (layout == "PE") {
    scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE,
                isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE,
                isDuplicate = FALSE)
} else {
    scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
}
scan_param <- ScanBamParam(flag = flag_filter, mapqFilter = min_mapq)
bam_files <- BamFileList(metadata$bam, yieldSize = 1000000L)
names(bam_files) <- metadata$key
if (!identical(names(bam_files), metadata$key)) {
    stop("BAM columns do not match samplesheet biological-sample order")
}

if (layout == "PE") {
    # fragments=FALSE is the strict proper-pair mode in GenomicAlignments:
    # readGAlignmentPairs drops singletons and counts each mate pair once.
    counted <- summarizeOverlaps(
        features = consensus, reads = bam_files, mode = "Union",
        ignore.strand = TRUE, inter.feature = TRUE,
        singleEnd = FALSE, fragments = FALSE, param = scan_param,
        BPPARAM = BiocParallel::SerialParam()
    )
} else {
    counted <- summarizeOverlaps(
        features = consensus, reads = bam_files, mode = "Union",
        ignore.strand = TRUE, inter.feature = TRUE,
        singleEnd = TRUE, fragments = FALSE, param = scan_param,
        BPPARAM = BiocParallel::SerialParam()
    )
}

raw_counts_all <- assay(counted)
if (ncol(raw_counts_all) != nrow(metadata)) {
    stop("Count-matrix column count does not match samplesheet metadata")
}
colnames(raw_counts_all) <- metadata$key
rownames(raw_counts_all) <- peak_ids
if (anyNA(raw_counts_all) || any(raw_counts_all < 0) ||
    any(abs(raw_counts_all - round(raw_counts_all)) > .Machine$double.eps^0.5)) {
    stop("Count matrix is not a complete non-negative integer matrix")
}
storage.mode(raw_counts_all) <- "integer"
write_matrix_gz(raw_counts_all, region_table,
                file.path(out_dir, "deseq2atac_raw_counts.tsv.gz"))

nonzero <- rowSums(raw_counts_all) > 0L
if (!any(nonzero)) stop("All DESeq2ATAC consensus-region counts are zero")
raw_counts <- raw_counts_all[nonzero, , drop = FALSE]
tested_regions <- region_table[nonzero, , drop = FALSE]

col_data <- data.frame(
    sample_id = metadata$sample_id,
    key = metadata$key,
    condition = metadata$condition,
    replicate = metadata$replicate,
    layout = metadata$layout,
    genome = metadata$genome,
    row.names = metadata$key,
    stringsAsFactors = FALSE
)
if (block_column != "") col_data[[block_column]] <- metadata[[block_column]]
design_formula <- if (block_column == "") {
    ~ condition
} else {
    as.formula(paste("~", block_column, "+ condition"))
}
design_matrix <- model.matrix(design_formula, data = col_data)
if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("DESeq2ATAC design is not full rank; block and condition may be confounded")
}

message("DESeq2ATAC ", peak_type, ": fitting ", deparse(design_formula), " and contrast ",
        numerator_condition, " vs ", reference_condition)
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = col_data,
                              design = design_formula)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE)
result <- results(
    dds,
    contrast = c("condition", numerator_condition, reference_condition),
    alpha = alpha
)
result_table <- cbind(
    tested_regions[match(rownames(result), tested_regions$peak_id), , drop = FALSE],
    as.data.frame(result, stringsAsFactors = FALSE)
)
result_table$direction <- ifelse(
    is.na(result_table$log2FoldChange), NA_character_,
    ifelse(result_table$log2FoldChange > 0,
           paste0("higher_in_", numerator_condition),
           ifelse(result_table$log2FoldChange < 0,
                  paste0("higher_in_", reference_condition), "unchanged"))
)
is_significant <- significant_rows(result_table, alpha)
significant_table <- result_table[is_significant, , drop = FALSE]

write_tsv_gz(result_table, file.path(out_dir, "deseq2atac_results_all.tsv.gz"))
write_tsv_gz(significant_table, file.path(out_dir, "deseq2atac_results_significant.tsv.gz"))
normalized_counts <- counts(dds, normalized = TRUE)
write_matrix_gz(normalized_counts, tested_regions,
                file.path(out_dir, "deseq2atac_normalized_counts.tsv.gz"))

size_factor_table <- data.frame(
    key = names(sizeFactors(dds)),
    peak_type = peak_type,
    condition = as.character(col_data[names(sizeFactors(dds)), "condition"]),
    size_factor = as.numeric(sizeFactors(dds)),
    stringsAsFactors = FALSE
)
write.table(size_factor_table, file.path(out_dir, "deseq2atac_size_factors.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(col_data, file.path(out_dir, "deseq2atac_sample_metadata.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

canonical_signal_count <- function(bam, sample_layout) {
    targets <- scanBamHeader(bam)[[1]]$targets
    allowed <- c(canonical_names(genome, TRUE), canonical_names(genome, FALSE))
    targets <- targets[names(targets) %in% allowed]
    if (!length(targets)) stop("No canonical chromosomes in BAM: ", bam)
    whole_chromosomes <- GRanges(names(targets), IRanges(1L, as.integer(targets)))
    records <- sum(countBam(bam, param = ScanBamParam(
        which = whole_chromosomes, flag = flag_filter, mapqFilter = min_mapq
    ))$records)
    if (sample_layout == "PE") {
        if (records %% 2L != 0L) stop("Odd proper-pair alignment count in ", bam)
        records / 2L
    } else records
}

signal_counts <- vapply(metadata$bam, canonical_signal_count, numeric(1),
                        sample_layout = layout)
if (any(!is.finite(signal_counts)) || any(signal_counts <= 0)) {
    stop("Every sample must contain at least one canonical filtered signal unit")
}
assigned_counts <- colSums(raw_counts_all)
library_table <- data.frame(
    key = metadata$key,
    peak_type = peak_type,
    condition = as.character(metadata$condition),
    layout = layout,
    signal_unit = ifelse(layout == "PE", "fragment", "read"),
    canonical_filtered_signal_count = as.numeric(signal_counts),
    assigned_to_consensus = as.numeric(assigned_counts[metadata$key]),
    fraction_assigned = as.numeric(assigned_counts[metadata$key]) / signal_counts,
    size_factor = as.numeric(sizeFactors(dds)[metadata$key]),
    stringsAsFactors = FALSE
)
write.table(library_table, file.path(out_dir, "deseq2atac_library_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("DESeq2ATAC ", peak_type, ": generating diagnostic figures")
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
transformed <- assay(vsd)
correlation <- cor(transformed, method = "pearson")
sample_distance <- as.matrix(dist(t(transformed)))
write.table(cbind(sample = rownames(correlation), correlation),
            file.path(out_dir, "deseq2atac_sample_correlation.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(cbind(sample = rownames(sample_distance), sample_distance),
            file.path(out_dir, "deseq2atac_sample_distances.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

library_plot_data <- rbind(
    data.frame(key = library_table$key, metric = "Canonical filtered signal units",
               value = library_table$canonical_filtered_signal_count),
    data.frame(key = library_table$key, metric = "DESeq2 size factor",
               value = library_table$size_factor)
)
library_plot <- ggplot(library_plot_data, aes(x = key, y = value, fill = metric)) +
    geom_col(show.legend = FALSE) + facet_wrap(~metric, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL,
         title = paste("DESeq2ATAC", peak_type, "library sizes and normalization")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
render_pair(file.path(plot_dir, "library_sizes_and_size_factors"),
            function() print(library_plot), width = 8, height = 7)

heat_palette <- colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(101)
render_pair(file.path(plot_dir, "sample_correlation"),
            function() plot_matrix(correlation,
                                   paste("Pearson correlation -", peak_type, "VST counts"),
                                   heat_palette),
            width = 8, height = 7)
distance_palette <- colorRampPalette(c("#FFFFFF", "#2166AC"))(101)
render_pair(file.path(plot_dir, "sample_distance"),
            function() plot_matrix(sample_distance,
                                   paste("Sample distance -", peak_type, "VST counts"),
                                   distance_palette),
            width = 8, height = 7)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"))
pca_data$key <- rownames(pca_data)
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = key)) +
    geom_point(size = 3) + geom_text(vjust = -0.7, size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", percent_variance[[1]], "% variance")) +
    ylab(paste0("PC2: ", percent_variance[[2]], "% variance")) +
    ggtitle(paste("DESeq2ATAC", peak_type, "PCA - VST counts")) + theme_bw(base_size = 11)
render_pair(file.path(plot_dir, "pca"), function() print(pca_plot), width = 8, height = 6)
write.table(pca_data, file.path(out_dir, "deseq2atac_pca_data.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

render_pair(file.path(plot_dir, "dispersion_estimates"),
            function() plotDispEsts(dds), width = 7, height = 6)
render_pair(file.path(plot_dir, "ma"),
            function() plotMA(result, alpha = alpha, ylim = c(-5, 5)), width = 7, height = 6)

volcano_data <- result_table[!is.na(result_table$pvalue) &
                             !is.na(result_table$log2FoldChange), , drop = FALSE]
volcano_data$minus_log10_pvalue <- -log10(pmax(volcano_data$pvalue, .Machine$double.xmin))
volcano_data$significance <- ifelse(significant_rows(volcano_data, alpha),
                                    paste0("FDR <= ", alpha), "Not significant")
volcano_plot <- ggplot(volcano_data, aes(log2FoldChange, minus_log10_pvalue,
                                        color = significance)) +
    geom_point(alpha = 0.55, size = 0.8) +
    scale_color_manual(values = c("Not significant" = "grey65",
                                  setNames("#B2182B", paste0("FDR <= ", alpha)))) +
    labs(x = paste0("log2 fold change: ", numerator_condition, " / ", reference_condition),
        y = "-log10(raw P-value)",
        title = paste("DESeq2ATAC", peak_type, "volcano plot"), color = NULL) +
    theme_bw(base_size = 11)
render_pair(file.path(plot_dir, "volcano"), function() print(volcano_plot), width = 7, height = 6)

if (nrow(significant_table) > 0L) {
    direction_counts <- table(significant_table$direction)
    render_pair(file.path(plot_dir, "significant_site_overview"), function() {
        barplot(direction_counts, col = c("#2166AC", "#B2182B"), las = 2,
                ylab = "Significant consensus regions",
                main = paste0("DESeq2ATAC ", peak_type, " sites at FDR <= ", alpha))
    }, width = 7, height = 6)
}

saveRDS(dds, file.path(out_dir, "deseq2atac_analysis_object.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "deseq2atac_session_info.txt"))

zero_message <- if (nrow(significant_table) == 0L) {
    paste0("No regions passed FDR <= ", alpha,
           "; the analysis completed successfully and all tested regions are retained.")
} else {
    paste(nrow(significant_table), "regions passed the configured FDR threshold.")
}
summary_lines <- c(
    paste("DESeq2ATAC", peak_type, "differential accessibility summary"),
    paste("Status: SUCCESS"),
    paste("Peak type:", peak_type),
    paste("Genome:", genome),
    paste("Layout:", layout),
    paste("Signal unit:", ifelse(layout == "PE", "properly paired fragment", "single-end read")),
    paste("Biological samples:", nrow(metadata)),
    paste("Conditions:", paste(names(condition_counts), condition_counts, sep = "=", collapse = ", ")),
    paste("Design:", deparse(design_formula)),
    paste("Contrast:", numerator_condition, "vs", reference_condition),
    paste("Positive log2 fold change means: higher accessibility in", numerator_condition),
    paste("Minimum sample support:", min_support),
    paste("Consensus regions:", nrow(region_table)),
    paste("Nonzero tested regions:", nrow(result_table)),
    paste("FDR threshold:", alpha),
    paste("Significant regions:", nrow(significant_table)),
    paste("Independent-filtered/NA adjusted P-values:", sum(is.na(result_table$padj))),
    paste("Size-factor method: DESeq2 poscounts"),
    paste("Counting method: GenomicAlignments summarizeOverlaps Union; inter.feature=TRUE"),
    zero_message
)
writeLines(summary_lines, file.path(out_dir, "deseq2atac_summary.txt"))
message("DESeq2ATAC ", peak_type, " complete: ", nrow(result_table), " tested; ",
        nrow(significant_table), " significant at FDR <= ", alpha)
