#!/usr/bin/env Rscript
# =============================================================================
# fastq2tracks v3.0 — DiffBind sample sheet preparation
# Usage: Rscript scripts/prepare_diffbind.R <samplesheet.csv> <filteredBamDir> <peaksDir> <outDir> <genome>
# Produces: diffbind_samplesheet_<genome>.csv
# Does NOT run DiffBind — only prepares files for downstream analysis.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript prepare_diffbind.R <ss.csv> <bamDir> <peaksDir> <outDir> <genome>")

ss_file   <- args[1]
bam_dir   <- args[2]
peaks_dir <- args[3]
out_dir   <- args[4]
genome    <- args[5]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ss <- read.csv(ss_file, stringsAsFactors = FALSE)
ss <- ss[tolower(ss$genome) == genome, ]
ss_ip  <- ss[tolower(ss$is_control) %in% c("false","0","no"), ]
ss_ctr <- ss[tolower(ss$is_control) %in% c("true","1","yes"), ]

make_bam_path <- function(sid, dir) {
  p <- file.path(dir, paste0(sid, "_dedup_blFilt.bam"))
  if (!file.exists(p)) p <- file.path(dir, paste0(sid, "_dedup.bam"))
  p
}

make_peak_path <- function(sid, peaks_dir, mode) {
  narrow <- file.path(peaks_dir, "per_replicate", sid, paste0(sid, "_peaks.narrowPeak"))
  broad  <- file.path(peaks_dir, "per_replicate", sid, paste0(sid, "_peaks.broadPeak"))
  if (grepl("narrow", mode, ignore.case=TRUE) && file.exists(narrow)) return(narrow)
  if (grepl("broad",  mode, ignore.case=TRUE) && file.exists(broad))  return(broad)
  if (file.exists(narrow)) return(narrow)
  if (file.exists(broad))  return(broad)
  return(NA_character_)
}

# Build DiffBind samplesheet
diffbind_ss <- lapply(seq_len(nrow(ss_ip)), function(i) {
  sid      <- ss_ip$sample_id[i]
  ctrl_id  <- ss_ip$control_id[i]
  mode     <- ss_ip$macs2_mode[i]

  bam_ip   <- make_bam_path(sid, bam_dir)
  peak     <- make_peak_path(sid, peaks_dir, mode)
  bam_ctrl <- if (!is.na(ctrl_id) & nchar(ctrl_id) > 0)
                make_bam_path(ctrl_id, bam_dir) else NA_character_

  data.frame(
    SampleID       = sid,
    Tissue         = ss_ip$cell_type[i],
    Factor         = ss_ip$factor[i],
    Condition      = ss_ip$condition[i],
    Treatment      = ss_ip$treatment[i],
    Replicate      = ss_ip$replicate[i],
    bamReads       = bam_ip,
    ControlID      = ifelse(is.na(ctrl_id) | nchar(ctrl_id)==0, NA_character_, ctrl_id),
    bamControl     = bam_ctrl,
    Peaks          = peak,
    PeakCaller     = ifelse(grepl("narrow", mode, ignore.case=TRUE), "narrow", "broad"),
    stringsAsFactors = FALSE
  )
})
diffbind_ss <- do.call(rbind, diffbind_ss)

out_file <- file.path(out_dir, paste0("diffbind_samplesheet_", genome, ".csv"))
write.csv(diffbind_ss, out_file, row.names = FALSE)

cat("DiffBind samplesheet written:", out_file, "\n")
cat("Samples:", nrow(diffbind_ss), "\n")

# Warn about missing files
missing <- c(
  diffbind_ss$bamReads[!file.exists(diffbind_ss$bamReads)],
  diffbind_ss$Peaks[!is.na(diffbind_ss$Peaks) & !file.exists(diffbind_ss$Peaks)]
)
if (length(missing) > 0) {
  cat("WARNING: Files not yet present (expected after full run):\n")
  cat(missing, sep="\n")
}
