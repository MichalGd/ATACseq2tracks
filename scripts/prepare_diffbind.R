#!/usr/bin/env Rscript
# =============================================================================
# ATACseq2tracks v3.0.4 — DiffBind sample sheet preparation (narrow + broad)
# Usage: Rscript scripts/prepare_diffbind.R <ss.csv> <bamDir> <peaksDir> <outDir> <genome>
# Produces: diffbind_samplesheet_<genome>_narrow.csv
#           diffbind_samplesheet_<genome>_broad.csv
#
# Robust mode: rows with empty/missing peak files are dropped with a warning
# rather than passed to DiffBind (which would crash on them).
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: prepare_diffbind.R <ss.csv> <bamDir> <peaksDir> <outDir> <genome>")

ss_file <- args[1]; bam_dir <- args[2]; peaks_dir <- args[3]
out_dir <- args[4]; genome  <- tolower(args[5])
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ss    <- read.csv(ss_file, stringsAsFactors = FALSE)
ss    <- ss[tolower(ss$genome) == genome, ]
ss_ip <- ss[tolower(ss$is_control) %in% c("false", "0", "no"), ]
ss_ip <- ss_ip[!duplicated(ss_ip$sample_id), ]   # collapse tech-rep rows

make_bam <- function(sid, rep, dir) {
    p <- file.path(dir, paste0(sid, "_bioR", rep, "_dedup_blFilt.bam"))
    if (!file.exists(p)) p <- file.path(dir, paste0(sid, "_bioR", rep, "_dedup.bam"))
    p
}

make_peak <- function(sid, rep, peaks_dir, ptype) {
    ext <- if (ptype == "broad") "broadPeak" else "narrowPeak"
    subdir <- ptype
    # Per-replicate layout
    p <- file.path(peaks_dir, "per_replicate", paste0(sid, "_bioR", rep), subdir,
                   paste0(sid, "_bioR", rep, "_peaks.", ext))
    if (file.exists(p)) return(p)
    # Flat layout fallback
    p2 <- file.path(peaks_dir, subdir, paste0(sid, "_peaks.", ext))
    if (file.exists(p2)) return(p2)
    return(NA_character_)
}

build_ss <- function(ptype) {
    df <- do.call(rbind, lapply(seq_len(nrow(ss_ip)), function(i) {
        sid     <- ss_ip$sample_id[i]
        rep     <- ss_ip$replicate[i]
        ctrl_id <- ss_ip$control_id[i]
        ctrl_rep <- if (!is.na(ctrl_id) && nchar(ctrl_id) > 0) {
            rows <- ss[ss$sample_id == ctrl_id, "replicate"]
            if (length(rows) > 0) rows[1] else rep
        } else NA_character_

        data.frame(
            SampleID   = sid,
            Tissue     = ss_ip$cell_type[i],
            Factor     = ss_ip$factor[i],
            Condition  = ss_ip$condition[i],
            Treatment  = ss_ip$treatment[i],
            Batch      = if ("batch" %in% names(ss_ip)) ss_ip$batch[i] else NA_character_,
            Replicate  = rep,
            bamReads   = make_bam(sid, rep, bam_dir),
            ControlID  = ifelse(is.na(ctrl_id) | nchar(ctrl_id) == 0, NA_character_, ctrl_id),
            bamControl = ifelse(is.na(ctrl_id) | nchar(ctrl_id) == 0, NA_character_,
                                make_bam(ctrl_id, ctrl_rep, bam_dir)),
            Peaks      = make_peak(sid, rep, peaks_dir, ptype),
            PeakCaller = ptype,
            stringsAsFactors = FALSE
        )
    }))

    n_before <- nrow(df)

    # Drop rows where BAM does not exist
    df <- df[file.exists(df$bamReads), ]
    n_no_bam <- n_before - nrow(df)
    if (n_no_bam > 0)
        message("WARN [", ptype, "]: dropped ", n_no_bam, " row(s) with missing BAM")

    # Drop rows with empty or missing peak files
    has_peaks <- !is.na(df$Peaks) &
                 file.exists(df$Peaks) &
                 (file.info(df$Peaks)$size > 0)
    n_no_peaks <- sum(!has_peaks)
    if (n_no_peaks > 0) {
        message("WARN [", ptype, "]: dropped ", n_no_peaks,
                " row(s) with empty/missing peaks: ",
                paste(df$SampleID[!has_peaks], collapse = ", "))
        df <- df[has_peaks, ]
    }

    if (nrow(df) == 0) {
        message("WARN [", ptype, "]: no valid samples remain -- output CSV will be empty")
    } else {
        message("[", ptype, "]: ", nrow(df), " / ", n_before, " samples retained for DiffBind")
    }

    df
}

for (ptype in c("narrow", "broad")) {
    result <- build_ss(ptype)
    out_csv <- file.path(out_dir, paste0("diffbind_samplesheet_", genome, "_", ptype, ".csv"))
    write.csv(result, out_csv, row.names = FALSE, quote = TRUE)
    message("Written: ", out_csv)
}
