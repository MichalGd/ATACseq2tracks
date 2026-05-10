#!/usr/bin/env bash
set -euo pipefail

# generate_multiqc_unified_report_v4.sh
#
# Fixes issues seen in v3b:
# - Unquoted heredoc caused bash to execute backticks inside the Rmd (command substitution),
#   leading to errors like:
#     ", asset_dir_rel, : command not found"
#     "*/plots/png/: No such file or directory"
# - Ensures the report is portable:
#     * mode=selfcontained embeds PNGs as data: URIs (no pandoc needed)
#     * mode=assets copies PNGs into *_assets/ and references them relatively
#
# Usage:
#   ./generate_multiqc_unified_report_v4.sh <InFolder> <OutFolder> [mode]
#
# mode:
#   selfcontained (default)  -> single HTML file, images embedded as base64 data URIs
#   assets                   -> HTML + *_assets folder (copy both together)

if [[ $# -lt 2 ]]; then
  cat <<'USAGE'
Usage: generate_multiqc_unified_report_v4.sh <InFolder> <OutFolder> [mode]

mode:
  selfcontained (default)  -> single HTML (images embedded)
  assets                   -> HTML + *_assets folder (copy both)

Example:
  ./generate_multiqc_unified_report_v4.sh /dysk2/groupFolders/micgdu/bioinformatics/run3_20aug2025_test \
    /dysk2/groupFolders/micgdu/bioinformatics/run3_20aug2025_test/reports/multiqc_summary_20260111 \
    selfcontained
USAGE
  exit 1
fi

InFolder="$1"
OutFolder="$2"
MODE="${3:-selfcontained}"

if [[ ! -d "$InFolder" ]]; then
  echo "ERROR: Input folder does not exist: $InFolder" >&2
  exit 1
fi
mkdir -p "$OutFolder"

if ! command -v R >/dev/null 2>&1; then
  echo "ERROR: R is not installed or not in PATH" >&2
  exit 1
fi

case "$MODE" in
  selfcontained|assets) ;;
  *)
    echo "ERROR: mode must be 'selfcontained' or 'assets' (got: $MODE)" >&2
    exit 1
    ;;
esac

REPORT_NAME="multiqc_unified_report_$(date +%Y%m%d_%H%M%S)"
REPORT_RMD="${OutFolder}/${REPORT_NAME}.Rmd"
HTML_OUT="${OutFolder}/${REPORT_NAME}.html"
ASSET_REL="${REPORT_NAME}_assets"
ASSET_DIR="${OutFolder}/${ASSET_REL}"

echo "[INFO] Input : $InFolder"
echo "[INFO] Output: $OutFolder"
echo "[INFO] Mode  : $MODE"
echo "[INFO] Assets: $ASSET_DIR"

# -------------------------------------------------------------------
# Copy plots into an assets folder next to the report.
# We include:
#   - MultiQC plots:   */multiQC/*_plots/png/*.png
#   - Any pipeline:    */plots/png/*.png  or */plots/PNG/*.png
# -------------------------------------------------------------------
mkdir -p "$ASSET_DIR"

mapfile -t PNGS < <(
  find "$InFolder" -type f -name "*.png" \
    \( \
      -path "*/multiQC/*_plots/png/*.png" -o \
      -path "*/multiQC/*_plots/PNG/*.png" -o \
      -path "*/plots/png/*.png"          -o \
      -path "*/plots/PNG/*.png" \
    \) \
    | sort
)

if [[ ${#PNGS[@]} -eq 0 ]]; then
  echo "[WARN] No PNGs found in expected plot folders. The report will render but without plots."
else
  echo "[INFO] Copying ${#PNGS[@]} PNGs into assets folder..."
  for f in "${PNGS[@]}"; do
    rel="${f#$InFolder/}"                 # path relative to InFolder
    dest="${ASSET_DIR}/${rel}"
    mkdir -p "$(dirname "$dest")"
    if [[ ! -f "$dest" ]] || ! cmp -s "$f" "$dest"; then
      cp -f "$f" "$dest"
    fi
  done
fi

# -------------------------------------------------------------------
# Write R Markdown (IMPORTANT: quoted heredoc to prevent bash expanding backticks / $)
# -------------------------------------------------------------------
cat > "$REPORT_RMD" <<'RMD'
---
title: "Unified ChIP-seq QC & Alignment Summary"
subtitle: "MultiQC + pipeline control plots (portable report)"
author: "Automated report"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    theme: bootstrap
    highlight: tango
    code_folding: hide
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center",
  out.width = "100%"
)

suppressPackageStartupMessages({
  library(readr)
  library(DT)
})

# Passed from bash via environment variables
base_folder <- normalizePath(Sys.getenv("UAM_INPUT_FOLDER"), winslash = "/", mustWork = FALSE)
out_folder  <- normalizePath(Sys.getenv("UAM_OUTPUT_FOLDER"), winslash = "/", mustWork = FALSE)
asset_rel   <- Sys.getenv("UAM_ASSET_REL")
mode        <- Sys.getenv("UAM_MODE", "assets")

# Ensure knitr evaluates with out_folder as working root (relative src works)
knitr::opts_knit$set(root.dir = out_folder)

clean_names <- function(x) {
  x <- gsub("^#\\s*", "", x)
  x <- gsub("[._]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

as_sample_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  colnames(df) <- clean_names(colnames(df))
  if (!("Sample" %in% colnames(df))) colnames(df)[1] <- "Sample"
  df
}

read_tsv_safe <- function(rel_path, label) {
  full_path <- file.path(base_folder, rel_path)
  if (!file.exists(full_path)) {
    cat("❌ ", label, " not found: ", full_path, "\n", sep = "")
    return(NULL)
  }
  tryCatch({
    df <- readr::read_tsv(full_path, show_col_types = FALSE, progress = FALSE)
    df <- as_sample_df(df)
    cat("✅ Loaded ", label, ": ", nrow(df), " rows\n", sep = "")
    df
  }, error = function(e) {
    cat("❌ Error reading ", label, ": ", e$message, "\n", sep = "")
    NULL
  })
}

pick_cols_regex <- function(df, patterns) {
  if (is.null(df)) return(NULL)
  keep <- rep(FALSE, ncol(df))
  for (p in patterns) keep <- keep | grepl(p, colnames(df), ignore.case = TRUE)
  keep[colnames(df) == "Sample"] <- TRUE
  df[, keep, drop = FALSE]
}

prefix_cols <- function(df, prefix) {
  if (is.null(df)) return(NULL)
  out <- df
  names(out) <- ifelse(names(out) == "Sample", "Sample", paste0(prefix, names(out)))
  out
}

merge_by_sample <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) return(NULL)
  Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE, sort = FALSE), dfs)
}

datatable_scroll <- function(df, caption) {
  if (is.null(df) || nrow(df) == 0) {
    cat("*No data available for ", caption, "*\n", sep = "")
    return(invisible(NULL))
  }
  DT::datatable(
    df,
    caption = caption,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      scrollY = "600px",
      fixedHeader = TRUE,
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel")
    ),
    extensions = c("Buttons", "FixedHeader"),
    rownames = FALSE,
    class = "cell-border stripe"
  )
}

# ---------- Portable plot embedding ----------
# assets mode: <img src="report_assets/.../plot.png">
# selfcontained: embed each PNG as data URI (no pandoc required)
base64_png <- function(p_full) {
  b64_bin <- Sys.which("base64")
  if (nzchar(b64_bin)) {
    # GNU coreutils: base64 -w 0; if -w unsupported, try without it
    out <- tryCatch(system2(b64_bin, c("-w", "0", shQuote(p_full)), stdout = TRUE, stderr = TRUE),
                    error = function(e) NULL)
    if (is.null(out) || any(grepl("invalid option", out, ignore.case = TRUE))) {
      out <- system2(b64_bin, c(shQuote(p_full)), stdout = TRUE)
    }
    return(paste(out, collapse = ""))
  }
  if (requireNamespace("base64enc", quietly = TRUE)) {
    return(base64enc::base64encode(p_full))
  }
  return(NA_character_)
}

img_src_for <- function(p_full, p_rel) {
  if (identical(mode, "selfcontained")) {
    b64 <- base64_png(p_full)
    if (!is.na(b64) && nzchar(b64)) {
      return(paste0("data:image/png;base64,", b64))
    }
  }
  p_rel
}

embed_pngs_from_asset_dir <- function(rel_dir, title, max_plots = 80) {
  asset_dir_rel <- file.path(asset_rel, rel_dir)
  asset_dir_full <- file.path(out_folder, asset_dir_rel)

  if (!dir.exists(asset_dir_full)) {
    cat("*Missing plot directory in assets:* `", asset_dir_rel, "`\n\n", sep = "")
    return(invisible(NULL))
  }

  pngs_full <- list.files(asset_dir_full, pattern = "\\.png$", full.names = TRUE)
  if (!length(pngs_full)) {
    cat("*No PNG plots in:* `", asset_dir_rel, "`\n\n", sep = "")
    return(invisible(NULL))
  }

  pngs_full <- sort(pngs_full)
  if (length(pngs_full) > max_plots) {
    cat("*Note:* showing first ", max_plots, " of ", length(pngs_full), " plots from `", asset_dir_rel, "`.\n\n", sep = "")
    pngs_full <- pngs_full[seq_len(max_plots)]
  }

  cat("### ", title, "\n\n", sep = "")
  for (p_full in pngs_full) {
    p_rel <- sub(paste0("^", out_folder, "/?"), "", p_full)
    cap <- tools::file_path_sans_ext(basename(p_full))
    cap <- gsub("[-_]+", " ", cap)
    cap <- gsub("\\s+", " ", cap)
    cap <- trimws(cap)

    src <- img_src_for(p_full, p_rel)

    cat("<div style=\"margin: 0 0 18px 0;\">\n")
    cat("  <div style=\"font-weight: 600; margin: 6px 0;\">", cap, "</div>\n", sep = "")
    cat("  <img src=\"", src, "\" style=\"max-width: 100%; height: auto; border: 1px solid #eee; border-radius: 6px;\">\n", sep = "")
    cat("</div>\n\n")
  }
}

discover_extra_plot_dirs_in_assets <- function() {
  dirs <- list.dirs(file.path(out_folder, asset_rel), recursive = TRUE, full.names = TRUE)
  dirs <- dirs[grepl("/plots/(png|PNG)$", dirs)]
  dirs <- dirs[!grepl("/multiQC/", dirs)]
  has_png <- vapply(dirs, function(d) length(list.files(d, pattern = "\\.png$", full.names = TRUE)) > 0, logical(1))
  dirs <- dirs[has_png]
  rels <- sub(paste0("^", file.path(out_folder, asset_rel), "/?"), "", dirs)
  sort(unique(rels))
}
```

# Overview

**Base folder (input):** `r base_folder`  
**Report folder (output):** `r out_folder`  
**Assets folder:** `r file.path(out_folder, asset_rel)`  
**Mode:** `r mode`

This report:
- builds **one per-sample table** combining key metrics from pre-trimming, post-trimming, alignment and deduplication,
- embeds **all PNG control plots** from MultiQC and pipeline plot folders.

---

# 1. Per-sample summary table (single table)

```{r tables}
pretrim_stats  <- read_tsv_safe("multiQC/multiQC_unTrimmed/multiQC_unTrimmed_data/multiqc_general_stats.txt", "Pre-trimming general stats")
posttrim_stats <- read_tsv_safe("multiQC/multiQC_trimmed/multiQC_trimmed_data/multiqc_general_stats.txt", "Post-trimming general stats")
align_stats    <- read_tsv_safe("multiQC/multiQC_alignments/multiQC_aligments_data/multiqc_general_stats.txt", "Alignment general stats")
dedup_stats    <- read_tsv_safe("multiQC/multiQC_deduplication/multiQC_deduplication_data/multiqc_general_stats.txt", "Deduplication general stats")

pre_key  <- pick_cols_regex(pretrim_stats,  c("sample$", "total", "gc", "percent.*duplicate|duplicates", "avg.*length|average.*length|sequence length"))
post_key <- pick_cols_regex(posttrim_stats, c("sample$", "total", "gc", "percent.*duplicate|duplicates", "avg.*length|average.*length|sequence length"))

align_key <- pick_cols_regex(align_stats, c(
  "sample$",
  "overall.*align|overall alignment",
  "mapped|mapping",
  "properly paired",
  "concordant",
  "total.*reads|reads"
))

dedup_key <- pick_cols_regex(dedup_stats, c(
  "sample$",
  "percent.*dup|percent duplication|duplicate",
  "estimated.*library",
  "read.*pairs.*examined|pairs examined",
  "optical"
))

samples_master <- merge_by_sample(list(
  prefix_cols(pre_key,  "PRE "),
  prefix_cols(post_key, "POST "),
  prefix_cols(align_key,"ALIGN "),
  prefix_cols(dedup_key,"DEDUP ")
))

datatable_scroll(samples_master, "Key QC + Alignment Metrics (one table for all samples)")
```

<details>
<summary><b>Full exported table (all columns from MultiQC general stats, very wide)</b></summary>

```{r full_table}
pre_full   <- prefix_cols(pretrim_stats,  "PRE FULL ")
post_full  <- prefix_cols(posttrim_stats, "POST FULL ")
align_full <- prefix_cols(align_stats,    "ALIGN FULL ")
dedup_full <- prefix_cols(dedup_stats,    "DEDUP FULL ")

samples_full <- merge_by_sample(list(pre_full, post_full, align_full, dedup_full))
datatable_scroll(samples_full, "Full MultiQC general stats across all stages (wide)")
```

</details>

---

# 2. Control plots (embedded / portable)

## 2.1 MultiQC - Pre-trimming FastQC plots

```{r pretrim_plots, results='asis'}
embed_pngs_from_asset_dir("multiQC/multiQC_unTrimmed/multiQC_unTrimmed_plots/png", "Pre-trimming FastQC (MultiQC)")
```

## 2.2 MultiQC - Post-trimming FastQC plots

```{r posttrim_plots, results='asis'}
embed_pngs_from_asset_dir("multiQC/multiQC_trimmed/multiQC_trimmed_plots/png", "Post-trimming FastQC (MultiQC)")
```

## 2.3 MultiQC - Alignment plots

```{r align_plots, results='asis'}
embed_pngs_from_asset_dir("multiQC/multiQC_alignments/multiQC_aligments_plots/png", "Alignment QC (MultiQC)")
```

## 2.4 MultiQC - Deduplication plots (if present)

```{r dedup_plots, results='asis'}
embed_pngs_from_asset_dir("multiQC/multiQC_deduplication/multiQC_deduplication_plots/png", "Deduplication QC (MultiQC)")
```

## 2.5 Extra pipeline plots (auto-discovered)

```{r extra_plots, results='asis'}
extra_dirs <- discover_extra_plot_dirs_in_assets()
if (!length(extra_dirs)) {
  cat("*No extra `*/plots/png/` directories were found (outside `multiQC/`).*\n")
} else {
  cat("*Discovered* ", length(extra_dirs), " extra plot directories.\n\n", sep = "")
  for (rel in extra_dirs) {
    embed_pngs_from_asset_dir(rel, paste0("Extra plots: ", rel))
  }
}
```

---

# 3. Copying this report

- **If mode = selfcontained:** copy just the HTML file anywhere (images are embedded).
- **If mode = assets:** copy the HTML **and** the folder `r asset_rel` (keep it next to the HTML).
RMD

echo "[INFO] Rmd created: $REPORT_RMD"
echo "[INFO] Rendering HTML report..."

export UAM_INPUT_FOLDER="$InFolder"
export UAM_OUTPUT_FOLDER="$OutFolder"
export UAM_ASSET_REL="$ASSET_REL"
export UAM_MODE="$MODE"

R --slave --no-restore --no-save -e "
tryCatch({
  req <- c('rmarkdown','knitr','readr','DT')
  missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    cat('ERROR: Missing R packages: ', paste(missing, collapse=', '), '\n', sep='')
    cat('Install them first (e.g. in R): install.packages(c(\"', paste(missing, collapse='\",\"'), '\"), repos=\"https://cran.r-project.org/\")\n', sep='')
    quit(status = 2)
  }
  library(rmarkdown)
  rmarkdown::render(
    input = Sys.getenv('UAM_REPORT_RMD', '$REPORT_RMD'),
    output_format = 'html_document',
    output_file = '$REPORT_NAME.html',
    output_dir = Sys.getenv('UAM_OUTPUT_FOLDER'),
    quiet = FALSE
  )
  cat('SUCCESS: $HTML_OUT\n')
}, error = function(e) {
  cat('ERROR: ', conditionMessage(e), '\n', sep='')
  quit(status = 1)
})
" 2>&1

if [[ "$MODE" == "selfcontained" ]]; then
  echo "[INFO] selfcontained mode: HTML includes images via data URIs."
  echo "[INFO] Assets folder kept for debugging: $ASSET_DIR"
else
  echo "[INFO] assets mode: copy HTML + assets folder together:"
  echo "       $HTML_OUT"
  echo "       $ASSET_DIR"
fi

echo "[DONE] HTML report: $HTML_OUT"
