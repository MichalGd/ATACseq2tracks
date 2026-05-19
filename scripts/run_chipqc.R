#!/usr/bin/env Rscript
# =============================================================================
# fastq2tracks v3.0 — ChIPQC module
# Usage: Rscript scripts/run_chipqc.R <samplesheet.csv> <filteredBamDir> <peaksDir> <outDir> <genome> <workers>
#   genome  : hg38 or mm39
#   workers : number of BiocParallel workers (default 20)
# Uses pre-built annotation + blacklist RDS objects from config.sh
# =============================================================================
suppressPackageStartupMessages({
  library(ChIPQC)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_chipqc.R <samplesheet.csv> <filteredBamDir> <peaksDir> <outDir> <genome> [workers]")
}

samplesheet_csv <- args[1]
bam_dir         <- args[2]
peaks_dir       <- args[3]
out_dir         <- args[4]
genome_key      <- args[5]
n_workers       <- as.integer(ifelse(length(args) >= 6, args[6], 20))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load config paths ---
config_sh <- readLines(file.path(dirname(getwd()), "fastq2tracks", "config", "config.sh"),
                        warn = FALSE)
get_cfg <- function(key) {
  line <- grep(paste0("^", key, "="), config_sh, value = TRUE)[1]
  if (is.na(line)) return(NULL)
  gsub(paste0("^", key, '="?|"?$'), "", line)
}

if (genome_key == "hg38") {
  anno_rds      <- get_cfg("CHIPQC_ANNOTATION_HG38")
  blacklist_rds <- get_cfg("CHIPQC_BLACKLIST_HG38_RDS")
} else if (genome_key == "mm39") {
  anno_rds      <- get_cfg("CHIPQC_ANNOTATION_MM39")
  blacklist_rds <- get_cfg("CHIPQC_BLACKLIST_MM39_RDS")
} else {
  stop("Unsupported genome: ", genome_key)
}

cat("Loading annotation:", anno_rds, "\n")
cat("Loading blacklist: ", blacklist_rds, "\n")
annotation <- readRDS(anno_rds)
blacklist   <- readRDS(blacklist_rds)

# --- Parse samplesheet ---
ss <- read.csv(samplesheet_csv, stringsAsFactors = FALSE)
ss <- ss[tolower(ss$genome) == genome_key, ]
ss <- ss[tolower(ss$is_control) %in% c("false", "0", "no"), ]

if (nrow(ss) == 0) stop("No IP samples found for genome: ", genome_key)

# Build ChIPQC samplesheet
chipqc_ss <- lapply(seq_len(nrow(ss)), function(i) {
  sid  <- ss$sample_id[i]
  
  # BAM: prefer blacklist-filtered
  bam <- file.path(bam_dir, paste0(sid, "_dedup_blFilt.bam"))
  if (!file.exists(bam)) bam <- file.path(bam_dir, paste0(sid, "_dedup.bam"))
  
  # Peaks: prefer narrowPeak, fall back to broadPeak
  peak <- file.path(peaks_dir, "per_replicate", sid,
                    paste0(sid, "_peaks.narrowPeak"))
  if (!file.exists(peak))
    peak <- file.path(peaks_dir, "per_replicate", sid,
                      paste0(sid, "_peaks.broadPeak"))
  if (!file.exists(peak)) peak <- NA_character_

  data.frame(
    SampleID   = sid,
    Tissue     = ss$cell_type[i],
    Factor     = ss$factor[i],
    Condition  = ss$condition[i],
    Treatment  = ss$treatment[i],
    Replicate  = ss$replicate[i],
    bamReads   = bam,
    Peaks      = peak,
    PeakCaller = ifelse(grepl("narrow", ss$macs2_mode[i], ignore.case=TRUE),
                        "narrow", "broad"),
    stringsAsFactors = FALSE
  )
})
chipqc_ss <- do.call(rbind, chipqc_ss)

# Remove samples with missing BAM
missing_bam <- !file.exists(chipqc_ss$bamReads)
if (any(missing_bam)) {
  cat("WARNING: Dropping samples with missing BAMs:\n")
  cat(chipqc_ss$SampleID[missing_bam], sep="\n")
  chipqc_ss <- chipqc_ss[!missing_bam, ]
}

cat("ChIPQC samples:", nrow(chipqc_ss), "\n")
write.csv(chipqc_ss, file.path(out_dir, "chipqc_samplesheet.csv"), row.names = FALSE)

# --- Configure BiocParallel ---
register(MulticoreParam(workers = n_workers, progressbar = TRUE))
cat("BiocParallel workers:", n_workers, "\n")

# --- Run ChIPQC ---
cat("Running ChIPQC...\n")
experiment <- ChIPQC(
  experiment  = chipqc_ss,
  annotation  = annotation,
  blacklist   = blacklist,
  chromosomes = if (genome_key == "hg38") {
    paste0("chr", c(1:22, "X", "Y"))
  } else {
    paste0("chr", c(1:19, "X", "Y"))
  },
  BPPARAM = bpparam()
)

# --- Reports ---
cat("Generating ChIPQC report...\n")
ChIPQCreport(
  object    = experiment,
  reportName= paste0("ChIPQC_report_", genome_key),
  reportFolder = out_dir,
  facet     = TRUE,
  colourBy  = "Condition"
)

# --- Summary tables ---
qc_metrics <- QCmetrics(experiment)
write.csv(as.data.frame(qc_metrics),
          file.path(out_dir, "chipqc_metrics_summary.csv"))

frip <- data.frame(
  SampleID = chipqc_ss$SampleID,
  FRiP     = round(frip(experiment), 4)
)
write.csv(frip, file.path(out_dir, "chipqc_frip.csv"), row.names = FALSE)

cat("ChIPQC complete. Outputs in:", out_dir, "\n")
