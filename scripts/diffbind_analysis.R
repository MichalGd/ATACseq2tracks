#!/usr/bin/env Rscript
# Universal all-pair DiffBind differential accessibility analysis.
# All valid samples contribute to the DiffBind consensus. Conditions with fewer
# than two biological samples are retained in consensus construction but are
# excluded from the differential model and from pairwise contrasts.

suppressPackageStartupMessages({
    library(DiffBind)
    library(ggplot2)
    library(GenomicRanges)
    library(rtracklayer)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
source(file.path(script_dir, "peak_annotation_helpers.R"))

as_bool <- function(value) tolower(trimws(as.character(value))) %in% c("true", "1", "yes")

safe_token <- function(value) {
    token <- iconv(as.character(value), to = "ASCII//TRANSLIT", sub = "_")
    token <- gsub("[^A-Za-z0-9._-]+", "_", token)
    token <- gsub("^_+|_+$", "", token)
    ifelse(token == "", "condition", substr(token, 1L, 60L))
}

resolve_condition_order <- function(observed, configured = "") {
    observed <- unique(trimws(as.character(observed)))
    requested <- trimws(unlist(strsplit(configured, ",", fixed = TRUE)))
    requested <- requested[requested != ""]
    if (anyDuplicated(requested)) stop("DIFFERENTIAL_CONDITION_ORDER contains duplicate names")
    unknown <- setdiff(requested, observed)
    if (length(unknown)) stop("Unknown condition(s) in DIFFERENTIAL_CONDITION_ORDER: ",
                              paste(unknown, collapse = ", "))
    c(requested, setdiff(observed, requested))
}

comparison_plan <- function(condition_order) {
    if (length(condition_order) < 2L) {
        return(data.frame(comparison_id = character(), numerator = character(),
                          reference = character(), stringsAsFactors = FALSE))
    }
    pairs <- utils::combn(condition_order, 2L)
    data.frame(
        comparison_id = sprintf(
            "%03d_%s_vs_%s", seq_len(ncol(pairs)),
            vapply(pairs[2L, ], safe_token, character(1)),
            vapply(pairs[1L, ], safe_token, character(1))
        ),
        numerator = pairs[2L, ],
        reference = pairs[1L, ],
        stringsAsFactors = FALSE
    )
}

write_tsv_gz <- function(x, path) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection), add = TRUE)
    write.table(x, connection, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE, na = "NA")
}

render_png <- function(path, code) {
    png(path, width = 1400, height = 1000, res = 150)
    tryCatch(code(), finally = dev.off())
}

summary_columns <- c(
    "module", "peak_type", "comparison_id", "numerator", "reference",
    "numerator_replicates", "reference_replicates", "consensus_regions",
    "tested_sites", "significant_sites", "higher_in_numerator",
    "higher_in_reference", "alpha", "min_abs_log2fc", "status",
    "results_all", "results_significant", "summary_file", "message"
)

empty_summary_row <- function() {
    as.data.frame(setNames(replicate(length(summary_columns), NA_character_, simplify = FALSE),
                           summary_columns), stringsAsFactors = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 15L) {
    stop(paste(
        "Usage: diffbind_analysis.R <diffbind_samplesheet.csv> <out_dir>",
        "[summits] [alpha] [min_abs_log2fc] [condition_order_csv] [genome]",
        "[broad|narrow] [annotate] [gtf] [ccre_bed] [ccre_source]",
        "[promoter_upstream] [promoter_downstream] [blacklist_bed]"
    ))
}

value_or <- function(index, default) if (length(args) >= index) args[[index]] else default
ss_file <- args[[1]]
out_dir <- args[[2]]
summits <- suppressWarnings(as.integer(value_or(3L, "100")))
alpha <- suppressWarnings(as.numeric(value_or(4L, "0.05")))
min_abs_log2fc <- suppressWarnings(as.numeric(value_or(5L, "0")))
configured_order <- value_or(6L, "")
genome <- tolower(value_or(7L, ""))
peak_type <- tolower(value_or(8L, tools::file_path_sans_ext(basename(ss_file))))
annotation_enabled <- as_bool(value_or(9L, "false"))
gtf_file <- value_or(10L, "")
ccre_file <- value_or(11L, "")
ccre_source <- value_or(12L, "")
promoter_upstream <- suppressWarnings(as.integer(value_or(13L, "2000")))
promoter_downstream <- suppressWarnings(as.integer(value_or(14L, "500")))
blacklist_file <- value_or(15L, "")

if (!genome %in% c("hg38", "mm39")) {
    inferred <- regmatches(basename(ss_file), regexpr("hg38|mm39", basename(ss_file)))
    if (length(inferred) && inferred %in% c("hg38", "mm39")) genome <- inferred
}

if (!file.exists(ss_file)) stop("Sample sheet not found: ", ss_file)
if (is.na(summits) || summits < 0L) stop("summits must be a non-negative integer")
if (is.na(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between zero and one")
if (is.na(min_abs_log2fc) || min_abs_log2fc < 0) {
    stop("min_abs_log2fc must be a non-negative number")
}
if (!genome %in% c("hg38", "mm39")) stop("DiffBind requires genome=hg38 or genome=mm39")
if (annotation_enabled && !file.exists(gtf_file)) stop("GTF annotation not found: ", gtf_file)
if (blacklist_file != "" && !file.exists(blacklist_file)) {
    stop("Blacklist annotation not found: ", blacklist_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
comparisons_dir <- file.path(out_dir, "comparisons")
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)
comparison_summary_file <- file.path(out_dir, "differential_accessibility_comparisons.tsv")
eligibility_file <- file.path(out_dir, "differential_accessibility_condition_eligibility.tsv")
module_summary_file <- file.path(out_dir, "diffbind_summary.txt")
log_file <- file.path(out_dir, "diffbind_log.txt")

sink(log_file, split = TRUE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!nrow(ss)) stop("Empty DiffBind sample sheet: ", ss_file)
if (!"Condition" %in% names(ss)) stop("DiffBind sample sheet requires Condition")
ss$Condition <- trimws(as.character(ss$Condition))
if (any(is.na(ss$Condition) | ss$Condition == "")) stop("Condition values must be non-empty")
if (anyDuplicated(ss$SampleID)) stop("DiffBind SampleID values must be unique biological samples")

observed_order <- unique(ss$Condition)
condition_counts <- table(factor(ss$Condition, levels = observed_order))
eligible_conditions <- observed_order[as.integer(condition_counts) >= 2L]
condition_order <- resolve_condition_order(observed_order, configured_order)
eligible_order <- condition_order[condition_order %in% eligible_conditions]
eligibility <- data.frame(
    condition = observed_order,
    biological_replicates = as.integer(condition_counts[observed_order]),
    included_in_consensus = TRUE,
    included_in_differential_model = observed_order %in% eligible_conditions,
    reason = ifelse(observed_order %in% eligible_conditions, "eligible",
                    "fewer_than_two_biological_replicates"),
    stringsAsFactors = FALSE
)
write.table(eligibility, eligibility_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("DiffBind all-sample consensus analysis for", ss_file, "\n")
cat("Samples:", nrow(ss), "\n")
cat("Conditions:", paste(names(condition_counts), condition_counts, sep = "=", collapse = ", "), "\n")
cat("Eligible model conditions:", ifelse(length(eligible_order), paste(eligible_order, collapse = ", "), "none"), "\n")
cat("Summit half-width:", summits, "bp\n")

# Count all samples first. This is the authoritative DiffBind all-sample
# consensus, including samples from conditions that are later model-ineligible.
dba_all <- dba(
    sampleSheet = ss_file,
    config = list(AnalysisMethod = DBA_DESEQ2, doBlacklist = FALSE,
                  doGreylist = FALSE)
)
dba_all <- dba.count(
    dba_all, minOverlap = 2, summits = summits,
    bUseSummarizeOverlaps = TRUE
)
all_consensus <- dba.peakset(dba_all, bRetrieve = TRUE)
all_consensus <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(all_consensus),
    ranges = GenomicRanges::ranges(all_consensus), strand = "*"
)
all_consensus <- annotation_canonicalize(all_consensus, genome)
if (blacklist_file != "") {
    blacklist <- annotation_canonicalize(rtracklayer::import(blacklist_file), genome)
    blacklist <- annotation_harmonize_seqnames(
        blacklist, unique(as.character(GenomicRanges::seqnames(all_consensus)))
    )
    all_consensus <- all_consensus[
        !GenomicRanges::overlapsAny(all_consensus, blacklist, ignore.strand = TRUE)
    ]
}
if (!length(all_consensus)) {
    stop("No canonical, non-blacklisted DiffBind consensus regions remain")
}
region_ids <- sprintf("diffbind_peak_%06d", seq_along(all_consensus))
names(all_consensus) <- region_ids
consensus_table <- data.frame(
    chrom = as.character(seqnames(all_consensus)),
    start = start(all_consensus) - 1L,
    end = end(all_consensus),
    region_id = region_ids,
    stringsAsFactors = FALSE
)
write.table(consensus_table, file.path(out_dir, "diffbind_consensus_peaks.bed"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
saveRDS(dba_all, file.path(out_dir, "diffbind_all_sample_count_object.rds"))

annotations <- NULL
if (annotation_enabled) {
    annotations <- annotate_peak_ranges(
        all_consensus, region_ids, genome, gtf_file, ccre_file, ccre_source,
        promoter_upstream, promoter_downstream
    )
    write_tsv_gz(annotations, file.path(out_dir, "diffbind_consensus_peak_annotations.tsv.gz"))
}

plan <- comparison_plan(eligible_order)
summary_rows <- list()
if (!nrow(plan)) {
    row <- empty_summary_row()
    row$module <- "DiffBind"
    row$peak_type <- peak_type
    row$consensus_regions <- length(all_consensus)
    row$alpha <- alpha
    row$min_abs_log2fc <- min_abs_log2fc
    row$status <- "SKIPPED"
    row$summary_file <- module_summary_file
    row$message <- "Fewer than two conditions have at least two biological replicates"
    write.table(row, comparison_summary_file, sep = "\t", row.names = FALSE,
                col.names = TRUE, quote = FALSE, na = "NA")
    writeLines(c(
        "DiffBind differential accessibility summary",
        "Status: SKIPPED",
        paste("Peak type:", peak_type),
        paste("All-sample consensus regions:", length(all_consensus)),
        paste("Conditions:", paste(names(condition_counts), condition_counts, sep = "=", collapse = ", ")),
        "Reason: fewer than two conditions have at least two biological replicates",
        "All samples still participated in consensus construction."
    ), module_summary_file)
    cat("DiffBind statistical analysis skipped: no eligible condition pair\n")
    while (sink.number() > 0L) sink()
    quit(save = "no", status = 0L)
}

# Recount only model-eligible samples on the fixed all-sample consensus. Passing
# summits=FALSE prevents a second recentering; filter=0 preserves the universe.
eligible_ss <- ss[ss$Condition %in% eligible_order, , drop = FALSE]
eligible_ss$Condition <- factor(eligible_ss$Condition, levels = eligible_order)
eligible_ss_file <- tempfile(pattern = "diffbind_eligible_", fileext = ".csv")
write.csv(eligible_ss, eligible_ss_file, row.names = FALSE, quote = TRUE)
dba_model <- dba(
    sampleSheet = eligible_ss_file,
    config = list(AnalysisMethod = DBA_DESEQ2, doBlacklist = FALSE,
                  doGreylist = FALSE)
)
dba_model <- dba.count(
    dba_model, peaks = all_consensus, summits = FALSE, filter = 0,
    bUseSummarizeOverlaps = TRUE
)
dba_model <- dba.normalize(dba_model)

for (index in seq_len(nrow(plan))) {
    contrast <- c("Condition", plan$numerator[[index]], plan$reference[[index]])
    if (index == 1L) {
        dba_model <- dba.contrast(
            dba_model, design = "~Condition", contrast = contrast,
            reorderMeta = list(Condition = eligible_order)
        )
    } else {
        dba_model <- dba.contrast(dba_model, contrast = contrast)
    }
}
contrast_table <- dba.show(dba_model, bContrasts = TRUE)
if (nrow(contrast_table) != nrow(plan)) {
    stop("DiffBind did not create the expected number of pairwise contrasts")
}
dba_model <- dba.analyze(dba_model, method = DBA_DESEQ2)
saveRDS(dba_model, file.path(out_dir, "diffbind_analysis_object.rds"))
write.table(contrast_table, file.path(out_dir, "diffbind_contrast_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

comparison_failures <- 0L
consensus_keys <- paste(seqnames(all_consensus), start(all_consensus), end(all_consensus), sep = ":")
for (index in seq_len(nrow(plan))) {
    comparison <- plan[index, , drop = FALSE]
    comparison_out <- file.path(comparisons_dir, comparison$comparison_id)
    dir.create(comparison_out, recursive = TRUE, showWarnings = FALSE)
    status_file <- file.path(comparison_out, "status.txt")
    summary_file <- file.path(comparison_out, "diffbind_summary.txt")
    results_all_file <- file.path(comparison_out, "diffbind_results_all_sites.csv")
    results_sig_file <- file.path(comparison_out, "diffbind_results.csv")
    row <- empty_summary_row()
    row$module <- "DiffBind"
    row$peak_type <- peak_type
    row$comparison_id <- comparison$comparison_id
    row$numerator <- comparison$numerator
    row$reference <- comparison$reference
    row$numerator_replicates <- as.integer(condition_counts[comparison$numerator])
    row$reference_replicates <- as.integer(condition_counts[comparison$reference])
    row$consensus_regions <- length(all_consensus)
    row$alpha <- alpha
    row$min_abs_log2fc <- min_abs_log2fc
    row$results_all <- results_all_file
    row$results_significant <- results_sig_file
    row$summary_file <- summary_file

    error_message <- tryCatch({
        result_all_gr <- dba.report(dba_model, contrast = index, method = DBA_DESEQ2,
                                    th = 1, fold = 0)
        result_sig_gr <- dba.report(dba_model, contrast = index, method = DBA_DESEQ2,
                                    th = alpha, fold = min_abs_log2fc)
        result_all <- as.data.frame(result_all_gr, stringsAsFactors = FALSE)
        result_sig <- as.data.frame(result_sig_gr, stringsAsFactors = FALSE)
        result_all$region_id <- region_ids[match(
            paste(result_all$seqnames, result_all$start, result_all$end, sep = ":"),
            consensus_keys
        )]
        result_sig$region_id <- region_ids[match(
            paste(result_sig$seqnames, result_sig$start, result_sig$end, sep = ":"),
            consensus_keys
        )]
        if (anyNA(result_all$region_id) || anyNA(result_sig$region_id)) {
            stop("DiffBind result coordinates do not match the fixed all-sample consensus")
        }
        if (!is.null(annotations)) {
            result_all <- annotation_join(result_all, annotations, "region_id")
            result_sig <- annotation_join(result_sig, annotations, "region_id")
        }
        write.csv(result_all, results_all_file, row.names = FALSE)
        write.csv(result_sig, results_sig_file, row.names = FALSE)
        if (length(result_sig_gr)) {
            rtracklayer::export(result_sig_gr,
                                file.path(comparison_out, "diffbind_significant_sites.bed"),
                                format = "BED")
        }

        render_png(file.path(comparison_out, "diffbind_pca.png"), function() {
            print(dba.plotPCA(dba_model, contrast = index, method = DBA_DESEQ2,
                              th = 1, attributes = DBA_CONDITION, label = DBA_ID))
        })
        render_png(file.path(comparison_out, "diffbind_heatmap.png"), function() {
            dba.plotHeatmap(dba_model, contrast = index, method = DBA_DESEQ2,
                            th = 1, correlations = FALSE)
        })
        render_png(file.path(comparison_out, "diffbind_ma.png"), function() {
            dba.plotMA(dba_model, contrast = index, method = DBA_DESEQ2)
        })
        render_png(file.path(comparison_out, "diffbind_volcano.png"), function() {
            print(dba.plotVolcano(dba_model, contrast = index, method = DBA_DESEQ2,
                                  th = alpha, fold = min_abs_log2fc,
                                  bReturnSites = FALSE))
        })

        fold_values <- if ("Fold" %in% names(result_sig)) result_sig$Fold else numeric()
        higher_numerator <- sum(fold_values > 0, na.rm = TRUE)
        higher_reference <- sum(fold_values < 0, na.rm = TRUE)
        row$tested_sites <- nrow(result_all)
        row$significant_sites <- nrow(result_sig)
        row$higher_in_numerator <- higher_numerator
        row$higher_in_reference <- higher_reference
        row$status <- "SUCCESS"
        row$message <- if (nrow(result_sig)) "completed" else "completed_with_zero_significant_sites"
        if (!is.null(annotations)) {
            write_peak_annotation_summary(result_sig, rep(TRUE, nrow(result_sig)),
                                          file.path(comparison_out, "annotation_summary.tsv"))
        }
        writeLines(c(
            "DiffBind differential accessibility comparison",
            "Status: SUCCESS",
            paste("Peak type:", peak_type),
            paste("Comparison ID:", comparison$comparison_id),
            paste("Contrast:", comparison$numerator, "vs", comparison$reference),
            paste("Positive Fold means: higher accessibility in", comparison$numerator),
            paste("Numerator biological replicates:", row$numerator_replicates),
            paste("Reference biological replicates:", row$reference_replicates),
            paste("All-sample consensus regions:", length(all_consensus)),
            paste("All tested sites:", nrow(result_all)),
            paste("FDR threshold:", alpha),
            paste("Minimum absolute log2 fold change:", min_abs_log2fc),
            paste("Significant sites:", nrow(result_sig)),
            paste("Higher in numerator:", higher_numerator),
            paste("Higher in reference:", higher_reference)
        ), summary_file)
        writeLines("SUCCESS", status_file)
        NULL
    }, error = function(error) conditionMessage(error))

    if (!is.null(error_message)) {
        comparison_failures <- comparison_failures + 1L
        row$status <- "FAILED"
        row$message <- error_message
        writeLines(c("FAILED", error_message), status_file)
        writeLines(c(
            "DiffBind differential accessibility comparison",
            "Status: FAILED",
            paste("Peak type:", peak_type),
            paste("Comparison ID:", comparison$comparison_id),
            paste("Contrast:", comparison$numerator, "vs", comparison$reference),
            paste("Error:", error_message)
        ), summary_file)
    }
    summary_rows[[length(summary_rows) + 1L]] <- row
}

comparison_summary <- do.call(rbind, summary_rows)
write.table(comparison_summary, comparison_summary_file, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE, na = "NA")
overall_status <- if (comparison_failures) "FAILED" else "SUCCESS"
writeLines(c(
    "DiffBind differential accessibility summary",
    paste("Status:", overall_status),
    paste("Peak type:", peak_type),
    paste("All biological samples in consensus:", nrow(ss)),
    paste("All-sample consensus regions:", length(all_consensus)),
    paste("Eligible model conditions:", paste(eligible_order, collapse = ", ")),
    paste("Excluded model conditions:", paste(setdiff(observed_order, eligible_order), collapse = ", ")),
    paste("Planned pairwise comparisons:", nrow(plan)),
    paste("Successful comparisons:", sum(comparison_summary$status == "SUCCESS")),
    paste("Failed comparisons:", comparison_failures),
    paste("Comparison summary:", comparison_summary_file)
), module_summary_file)

# Preserve the established root-level two-condition files while the universal
# layout always keeps the authoritative result in comparisons/<comparison_id>/.
if (nrow(plan) == 1L && comparison_summary$status[[1]] == "SUCCESS") {
    legacy_source <- file.path(comparisons_dir, plan$comparison_id[[1]])
    file.copy(file.path(legacy_source, "diffbind_results.csv"),
              file.path(out_dir, "diffbind_results.csv"), overwrite = TRUE)
    file.copy(file.path(legacy_source, "diffbind_results_all_sites.csv"),
              file.path(out_dir, "diffbind_results_all_sites.csv"), overwrite = TRUE)
    for (plot_name in c("diffbind_pca.png", "diffbind_heatmap.png",
                        "diffbind_ma.png", "diffbind_volcano.png")) {
        file.copy(file.path(legacy_source, plot_name), file.path(out_dir, plot_name),
                  overwrite = TRUE)
    }
}

unlink(eligible_ss_file)
while (sink.number() > 0L) sink()
cat("DiffBind complete:", nrow(plan), "planned comparisons;",
    comparison_failures, "failed\n")
if (comparison_failures) quit(save = "no", status = 1L)
