#!/usr/bin/env Rscript
# =============================================================================
# fastq2tracks v3.0.4 — ChIPQC module
#
# Usage:
#   Rscript scripts/run_chipqc.R \
#       <samplesheet.csv> <bamDir> <peaksDir> <outDir> \
#       <genome: hg38|mm39> [workers:20] [peakType:narrow|broad] \
#       [/absolute/path/to/config.conf]
#
# Robust mode: samples with empty/missing peak files are silently dropped
# before ChIPQC is run. If ALL samples lack peaks the script exits 0 with
# a warning (does not abort the pipeline).
# =============================================================================
suppressPackageStartupMessages({
    library(ChIPQC); library(BiocParallel)
    library(dplyr);  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop(
    "Usage: Rscript run_chipqc.R <ss.csv> <bamDir> <peaksDir> <outDir> <genome> [workers] [peakType] [config.conf]")

samplesheet_csv <- args[1]; bam_dir  <- args[2]
peaks_dir       <- args[3]; out_dir  <- args[4]
genome_key      <- tolower(args[5])
n_workers       <- as.integer(ifelse(length(args) >= 6 && nchar(args[6]) > 0, args[6], 20))
peak_type_pref  <- ifelse(length(args) >= 7 && nchar(args[7]) > 0, tolower(args[7]), "narrow")

# ---- Resolve config.conf ----------------------------------------------------
config_sh_path <- if (length(args) >= 8 && nchar(args[8]) > 0) {
    args[8]
} else if (nchar(Sys.getenv("F2T_CONFIG")) > 0) {
    Sys.getenv("F2T_CONFIG")
} else {
    script_dir <- tryCatch(
        normalizePath(dirname(sys.frame(1)$ofile)),
        error = function(e) getwd()
    )
    file.path(dirname(script_dir), "config", "config.conf")
}
if (!file.exists(config_sh_path))
    stop("config.conf not found: ", config_sh_path,
         "\nPass as 8th argument or export F2T_CONFIG.")
cat("Config:", config_sh_path, "\n")

cfg_lines <- readLines(config_sh_path, warn = FALSE)
get_cfg <- function(key) {
    line <- grep(paste0("^\\s*", key, "\\s*="), cfg_lines, value = TRUE)[1]
    if (is.na(line)) return(NULL)
    val <- sub(paste0("^\\s*", key, "\\s*=\\s*['\"]?"), "", line)
    val <- sub("['\"]?\\s*(#.*)?$", "", val); trimws(val)
}

anno_rds <- blacklist_rds <- NULL
if (genome_key == "hg38") {
    anno_rds      <- get_cfg("CHIPQC_ANNOTATION_HG38")
    blacklist_rds <- get_cfg("CHIPQC_BLACKLIST_HG38_RDS")
} else {
    anno_rds      <- get_cfg("CHIPQC_ANNOTATION_MM39")
    blacklist_rds <- get_cfg("CHIPQC_BLACKLIST_MM39_RDS")
}

# ---- Build ChIPQC samplesheet -----------------------------------------------
ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE)
ss <- ss[tolower(ss$genome) == genome_key, ]
ss_ip <- ss[tolower(ss$is_control) %in% c("false", "0", "no"), ]
ss_ip <- ss_ip[!duplicated(paste0(ss_ip$sample_id, "_bioR", ss_ip$replicate)), ]

make_peak_path <- function(sid, rep, ptype) {
    subdir <- if (ptype == "narrow") "narrow" else "broad"
    ext    <- if (ptype == "narrow") "narrowPeak" else "broadPeak"
    # Per-replicate path (fastq2tracks layout)
    p <- file.path(peaks_dir, "per_replicate", paste0(sid, "_bioR", rep), subdir,
                   paste0(sid, "_bioR", rep, "_peaks.", ext))
    if (file.exists(p)) return(p)
    # Flat layout fallback
    p2 <- file.path(peaks_dir, subdir, paste0(sid, "_peaks.", ext))
    if (file.exists(p2)) return(p2)
    return(NA_character_)
}

chipqc_ss <- do.call(rbind, lapply(seq_len(nrow(ss_ip)), function(i) {
    sid <- ss_ip$sample_id[i]
    rep <- ss_ip$replicate[i]
    bam_canon <- file.path(bam_dir, paste0(sid, "_bioR", rep, "_canonical.bam"))
    bam <- if (file.exists(bam_canon)) bam_canon else file.path(bam_dir, paste0(sid, "_bioR", rep, "_dedup_blFilt.bam"))
    pk  <- make_peak_path(sid, rep, peak_type_pref)
    data.frame(SampleID   = paste0(sid, "_bioR", rep),
               Tissue     = ss_ip$cell_type[i],
               Factor     = ss_ip$factor[i],
               Condition  = ss_ip$condition[i],
               Treatment  = ss_ip$treatment[i],
               Replicate  = rep,
               bamReads   = bam,
               Peaks      = pk,
               PeakCaller = peak_type_pref,
               stringsAsFactors = FALSE)
}))

# ---- ROBUST FILTER: drop rows with empty/missing BAM or peak files ----------
n_before <- nrow(chipqc_ss)

# Drop missing BAM
chipqc_ss <- chipqc_ss[file.exists(chipqc_ss$bamReads), ]
n_no_bam  <- n_before - nrow(chipqc_ss)
if (n_no_bam > 0)
    message("WARN: dropped ", n_no_bam, " sample(s) with missing BAM files")

# Drop missing or empty peak files
has_peaks <- !is.na(chipqc_ss$Peaks) &
             file.exists(chipqc_ss$Peaks) &
             (file.info(chipqc_ss$Peaks)$size > 0)
n_no_peaks <- sum(!has_peaks)
if (n_no_peaks > 0) {
    message("WARN: dropped ", n_no_peaks, " sample(s) with empty/missing peak files: ",
            paste(chipqc_ss$SampleID[!has_peaks], collapse = ", "))
    chipqc_ss <- chipqc_ss[has_peaks, ]
}

# If nothing left, exit gracefully
if (nrow(chipqc_ss) == 0) {
    message("WARN: no samples with valid BAM + peaks for genome=", genome_key,
            " -- skipping ChIPQC (pipeline continues)")
    quit(save = "no", status = 0)
}

message("ChIPQC will run on ", nrow(chipqc_ss), " / ", n_before, " samples")

# ---- Write filtered samplesheet for reference -------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ss_out <- file.path(out_dir, paste0("chipqc_samplesheet.csv"))
write.csv(chipqc_ss, ss_out, row.names = FALSE, quote = TRUE)

# ---- Run ChIPQC -------------------------------------------------------------
register(MulticoreParam(workers = min(n_workers, nrow(chipqc_ss))))

anno_arg      <- if (!is.null(anno_rds)      && file.exists(anno_rds))      readRDS(anno_rds)      else genome_key
blacklist_arg <- if (!is.null(blacklist_rds) && file.exists(blacklist_rds)) readRDS(blacklist_rds) else NULL

experiment <- tryCatch(
    ChIPQC(chipqc_ss, annotation = anno_arg, blacklist = blacklist_arg,
        chromosomes = if (genome_key == "hg38") paste0("chr", c(1:22,"X","Y")) else paste0("chr", c(1:19,"X","Y"))),
    error = function(e) {
        message("ERROR in ChIPQC(): ", conditionMessage(e))
        quit(save = "no", status = 1)
    }
)

ChIPQCreport(experiment,
             reportName  = paste0("ChIPQC_", genome_key),
             reportFolder = out_dir,
             facetBy     = c("Condition", "Factor"))

message("ChIPQC complete. Report in: ", out_dir)
