# Shared, dependency-light peak annotation helpers for ATACseq2tracks.
# Requires GenomicRanges, IRanges and rtracklayer (already workflow dependencies).

annotation_canonical_names <- function(genome, with_chr = TRUE) {
    maximum <- switch(tolower(genome), hg38 = 22L, mm39 = 19L,
                      stop("Unsupported genome for annotation: ", genome))
    values <- c(as.character(seq_len(maximum)), "X", "Y")
    if (with_chr) paste0("chr", values) else values
}

annotation_canonicalize <- function(gr, genome) {
    allowed <- c(annotation_canonical_names(genome, TRUE),
                 annotation_canonical_names(genome, FALSE))
    gr[as.character(GenomicRanges::seqnames(gr)) %in% allowed]
}

annotation_harmonize_seqnames <- function(gr, reference_names) {
    if (!length(gr)) return(gr)
    current <- as.character(GenomicRanges::seqnames(gr))
    if (length(intersect(unique(current), unique(reference_names)))) return(gr)
    reference_has_chr <- any(grepl("^chr", reference_names))
    current_has_chr <- any(grepl("^chr", current))
    if (reference_has_chr && !current_has_chr) {
        GenomeInfoDb::seqlevels(gr) <- paste0("chr", GenomeInfoDb::seqlevels(gr))
    } else if (!reference_has_chr && current_has_chr) {
        GenomeInfoDb::seqlevels(gr) <- sub("^chr", "", GenomeInfoDb::seqlevels(gr))
    }
    gr
}

annotation_first_field <- function(gr, candidates, fallback = NA_character_) {
    metadata <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
    for (candidate in candidates) {
        if (candidate %in% names(metadata)) {
            value <- as.character(metadata[[candidate]])
            value[is.na(value) | value == ""] <- fallback
            return(value)
        }
    }
    rep(fallback, length(gr))
}

annotation_collapse_hits <- function(query_hits, values, n_query) {
    output <- rep(NA_character_, n_query)
    if (!length(query_hits)) return(output)
    split_values <- split(values, query_hits)
    collapsed <- vapply(split_values, function(x) {
        x <- sort(unique(x[!is.na(x) & x != ""]))
        if (length(x)) paste(x, collapse = ";") else NA_character_
    }, character(1))
    output[as.integer(names(collapsed))] <- collapsed
    output
}

annotation_read_ccre <- function(path, genome, reference_names) {
    compressed <- grepl("[.]gz$", path, ignore.case = TRUE)
    connection <- if (compressed) gzfile(path, "rt") else path
    if (compressed) on.exit(close(connection), add = TRUE)
    table <- utils::read.delim(
        connection, header = FALSE, sep = "\t", quote = "", comment.char = "",
        fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE
    )
    if (ncol(table) < 3L) stop("cCRE file has fewer than three BED columns: ", path)
    starts <- suppressWarnings(as.integer(table[[2]]))
    ends <- suppressWarnings(as.integer(table[[3]]))
    valid <- !is.na(starts) & !is.na(ends) & starts >= 0L & ends > starts
    table <- table[valid, , drop = FALSE]
    starts <- starts[valid]
    ends <- ends[valid]
    if (!nrow(table)) stop("cCRE file has no valid BED intervals: ", path)

    identifiers <- if (ncol(table) >= 4L) as.character(table[[4]]) else {
        paste0("ccre_", seq_len(nrow(table)))
    }
    identifiers[is.na(identifiers) | identifiers == ""] <-
        paste0("ccre_", which(is.na(identifiers) | identifiers == ""))

    known <- paste(
        c("PLS", "pELS", "dELS", "DNase-H3K4me3", "CA-CTCF",
          "CA-H3K4me3", "CA-TF", "CTCF-only", "CTCF-bound", "CA", "TF"),
        collapse = "|"
    )
    candidate_columns <- if (ncol(table) >= 4L) 4:ncol(table) else integer()
    class_column <- NA_integer_
    if (length(candidate_columns)) {
        scores <- vapply(candidate_columns, function(column) {
            sum(grepl(paste0("(^|[^A-Za-z0-9])(", known, ")([^A-Za-z0-9]|$)"),
                      as.character(table[[column]]), perl = TRUE), na.rm = TRUE)
        }, numeric(1))
        if (max(scores) > 0L) class_column <- candidate_columns[[which.max(scores)]]
    }
    raw_class <- if (!is.na(class_column)) as.character(table[[class_column]]) else {
        rep("unclassified", nrow(table))
    }

    normalize_class <- function(value) {
        if (is.na(value) || value == "") return("unclassified")
        hits <- regmatches(value, gregexpr(known, value, perl = TRUE))[[1]]
        hits <- unique(hits[hits != ""])
        if (length(hits)) paste(hits, collapse = ";") else trimws(value)
    }
    ccre_class <- vapply(raw_class, normalize_class, character(1))

    ranges <- GenomicRanges::GRanges(
        seqnames = as.character(table[[1]]),
        ranges = IRanges::IRanges(start = starts + 1L, end = ends)
    )
    S4Vectors::mcols(ranges)$ccre_id <- identifiers
    S4Vectors::mcols(ranges)$ccre_class <- ccre_class
    ranges <- annotation_canonicalize(ranges, genome)
    annotation_harmonize_seqnames(ranges, reference_names)
}

annotate_peak_ranges <- function(peaks, region_ids, genome, gtf_file,
                                 ccre_file = "", ccre_source = "",
                                 promoter_upstream = 2000L,
                                 promoter_downstream = 500L) {
    if (length(peaks) != length(region_ids)) stop("Annotation IDs do not match peak ranges")
    if (!file.exists(gtf_file)) stop("GTF annotation not found: ", gtf_file)
    promoter_upstream <- as.integer(promoter_upstream)
    promoter_downstream <- as.integer(promoter_downstream)
    if (is.na(promoter_upstream) || promoter_upstream < 0L ||
        is.na(promoter_downstream) || promoter_downstream < 0L) {
        stop("Promoter upstream/downstream values must be non-negative integers")
    }

    peaks <- annotation_canonicalize(peaks, genome)
    if (length(peaks) != length(region_ids)) {
        stop("Noncanonical intervals reached peak annotation")
    }
    reference_names <- unique(as.character(GenomicRanges::seqnames(peaks)))
    gtf <- annotation_canonicalize(rtracklayer::import(gtf_file), genome)
    gtf <- annotation_harmonize_seqnames(gtf, reference_names)
    feature_type <- if ("type" %in% names(S4Vectors::mcols(gtf))) {
        as.character(S4Vectors::mcols(gtf)$type)
    } else rep("", length(gtf))

    genes <- gtf[feature_type == "gene"]
    exons <- gtf[feature_type == "exon"]
    if (!length(genes)) {
        stop("GTF contains no gene features; a gene-level GTF is required: ", gtf_file)
    }
    if (!length(exons)) warning("GTF contains no exon features: ", gtf_file)

    gene_ids <- annotation_first_field(genes, c("gene_id", "ID", "gene"))
    gene_names <- annotation_first_field(genes, c("gene_name", "Name"))
    gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]

    tss <- GenomicRanges::resize(genes, width = 1L, fix = "start", ignore.strand = FALSE)
    S4Vectors::mcols(tss)$gene_id <- gene_ids
    S4Vectors::mcols(tss)$gene_name <- gene_names
    plus <- as.character(GenomicRanges::strand(tss)) != "-"
    promoter_start <- ifelse(
        plus,
        GenomicRanges::start(tss) - promoter_upstream,
        GenomicRanges::start(tss) - promoter_downstream
    )
    promoter_end <- ifelse(
        plus,
        GenomicRanges::start(tss) + promoter_downstream,
        GenomicRanges::start(tss) + promoter_upstream
    )
    promoters <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(tss),
        ranges = IRanges::IRanges(start = pmax(1L, promoter_start), end = promoter_end),
        strand = GenomicRanges::strand(tss)
    )

    promoter_hit <- IRanges::overlapsAny(peaks, promoters, ignore.strand = TRUE)
    exon_hit <- if (length(exons)) {
        IRanges::overlapsAny(peaks, exons, ignore.strand = TRUE)
    } else rep(FALSE, length(peaks))
    gene_hit <- IRanges::overlapsAny(peaks, genes, ignore.strand = TRUE)

    # In this simple GTF classifier, intron means gene-body sequence not covered
    # by any annotated exon. The separate gene_body_overlap column remains
    # non-exclusive and makes that definition transparent.
    intron_hit <- gene_hit & !exon_hit & !promoter_hit
    gene_context <- ifelse(
        promoter_hit, "promoter",
        ifelse(exon_hit, "exon",
               ifelse(intron_hit, "intron",
                      ifelse(gene_hit, "other_gene_body", "distal_intergenic")))
    )

    nearest <- GenomicRanges::distanceToNearest(peaks, tss, ignore.strand = TRUE)
    nearest_subject <- rep(NA_integer_, length(peaks))
    nearest_subject[S4Vectors::queryHits(nearest)] <- S4Vectors::subjectHits(nearest)
    nearest_gene_id <- nearest_gene_name <- rep(NA_character_, length(peaks))
    nearest_tss_distance <- rep(NA_integer_, length(peaks))
    has_nearest <- !is.na(nearest_subject)
    nearest_gene_id[has_nearest] <- gene_ids[nearest_subject[has_nearest]]
    nearest_gene_name[has_nearest] <- gene_names[nearest_subject[has_nearest]]
    peak_midpoint <- floor((GenomicRanges::start(peaks) + GenomicRanges::end(peaks)) / 2)
    tss_position <- GenomicRanges::start(tss)[nearest_subject[has_nearest]]
    signed_distance <- peak_midpoint[has_nearest] - tss_position
    nearest_strand <- as.character(GenomicRanges::strand(tss))[nearest_subject[has_nearest]]
    signed_distance[nearest_strand == "-"] <- -signed_distance[nearest_strand == "-"]
    nearest_tss_distance[has_nearest] <- as.integer(signed_distance)

    result <- data.frame(
        region_id = as.character(region_ids),
        chrom = as.character(GenomicRanges::seqnames(peaks)),
        start = GenomicRanges::start(peaks) - 1L,
        end = GenomicRanges::end(peaks),
        gene_context = gene_context,
        promoter_overlap = promoter_hit,
        exon_overlap = exon_hit,
        intron_overlap = intron_hit,
        gene_body_overlap = gene_hit,
        nearest_gene_id = nearest_gene_id,
        nearest_gene_name = nearest_gene_name,
        nearest_tss_distance_bp = nearest_tss_distance,
        promoter_definition = paste0("TSS-", promoter_upstream, "/+", promoter_downstream, "bp"),
        ccre_overlap = FALSE,
        ccre_primary_id = NA_character_,
        ccre_primary_class = NA_character_,
        ccre_all_ids = NA_character_,
        ccre_all_classes = NA_character_,
        enhancer_like = FALSE,
        ccre_source = ifelse(ccre_source == "", NA_character_, ccre_source),
        stringsAsFactors = FALSE
    )

    if (ccre_file != "") {
        if (!file.exists(ccre_file)) stop("Configured cCRE annotation not found: ", ccre_file)
        ccre <- annotation_read_ccre(ccre_file, genome, reference_names)
        hits <- IRanges::findOverlaps(peaks, ccre, ignore.strand = TRUE)
        if (length(hits)) {
            query <- S4Vectors::queryHits(hits)
            subject <- S4Vectors::subjectHits(hits)
            overlap_width <- GenomicRanges::width(GenomicRanges::pintersect(
                peaks[query], ccre[subject], ignore.strand = TRUE
            ))
            hit_ids <- as.character(S4Vectors::mcols(ccre)$ccre_id[subject])
            hit_classes <- as.character(S4Vectors::mcols(ccre)$ccre_class[subject])
            priority_values <- c(PLS = 1L, pELS = 2L, dELS = 3L, `CA-CTCF` = 4L,
                                 `CTCF-only` = 5L, `CTCF-bound` = 6L,
                                 `DNase-H3K4me3` = 7L, `CA-H3K4me3` = 8L,
                                 `CA-TF` = 9L, CA = 10L, TF = 11L)
            first_class <- sub(";.*$", "", hit_classes)
            priority <- unname(priority_values[first_class])
            priority[is.na(priority)] <- 999L
            ordered <- order(query, -overlap_width, priority, hit_ids)
            primary_rows <- ordered[!duplicated(query[ordered])]
            primary_query <- query[primary_rows]
            result$ccre_primary_id[primary_query] <- hit_ids[primary_rows]
            result$ccre_primary_class[primary_query] <- hit_classes[primary_rows]
            result$ccre_overlap[unique(query)] <- TRUE
            result$ccre_all_ids <- annotation_collapse_hits(query, hit_ids, length(peaks))
            result$ccre_all_classes <- annotation_collapse_hits(query, hit_classes, length(peaks))
            result$enhancer_like <- grepl("(^|;)(pELS|dELS)(;|$)",
                                          result$ccre_all_classes)
            result$enhancer_like[is.na(result$enhancer_like)] <- FALSE
        }
    }
    result
}

annotation_join <- function(table, annotations, id_column) {
    matches <- match(as.character(table[[id_column]]), annotations$region_id)
    if (anyNA(matches)) stop("Result rows could not be matched to shared peak annotations")
    annotation_columns <- setdiff(names(annotations), c("region_id", "chrom", "start", "end"))
    cbind(table, annotations[matches, annotation_columns, drop = FALSE])
}

write_peak_annotation_summary <- function(result_table, significant, path) {
    selected <- result_table[significant, , drop = FALSE]
    summarize_one <- function(column, annotation_type) {
        if (!column %in% names(selected) || !nrow(selected)) return(NULL)
        values <- as.character(selected[[column]])
        values[is.na(values) | values == ""] <- "unannotated"
        counts <- sort(table(values), decreasing = TRUE)
        data.frame(annotation_type = annotation_type, category = names(counts),
                   significant_sites = as.integer(counts), stringsAsFactors = FALSE)
    }
    summary <- rbind(
        summarize_one("gene_context", "gene_context"),
        summarize_one("ccre_primary_class", "ccre_primary_class")
    )
    if (is.null(summary)) {
        summary <- data.frame(annotation_type = character(), category = character(),
                              significant_sites = integer(), stringsAsFactors = FALSE)
    }
    utils::write.table(summary, path, sep = "\t", row.names = FALSE,
                       col.names = TRUE, quote = FALSE)
}
