#!/bin/bash

# Enhanced R Markdown Report Generator for DNAfastqBigWig Pipeline with PDF Generation

# This script generates a comprehensive report documenting the execution of DNAfastqBigWig_human_main_v5_31aug2025.sh

# and automatically renders it to both HTML and PDF formats

# Usage: ./generate_pipeline_report_pdf.sh [format]

# Example: ./generate_pipeline_report_pdf.sh /path/to/results pipeline_report_$(date +%Y%m%d) both

# Formats: html, pdf, both (default: both)

if [[ $# -lt 2 ]]; then

echo "Usage: $0 [format]"

echo "Example: $0 /dysk2/results pipeline_report_$(date +%Y%m%d) both"

echo "Formats: html, pdf, both (default: both)"

exit 1

fi

OUTPUT_FOLDER="$1"

REPORT_NAME="$2"

FORMAT="${3:-both}"

REPORT_FILE="${OUTPUT_FOLDER}/${REPORT_NAME}.Rmd"

# Validate format parameter

if [[ ! "$FORMAT" =~ ^(html|pdf|both)$ ]]; then

echo "ERROR: Format must be 'html', 'pdf', or 'both'"

exit 1

fi

# Get system information

HOSTNAME=$(hostname)

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

USER=$(whoami)

echo "=== DNAfastqBigWig Pipeline Report Generator with PDF Support ==="

echo "Output folder: $OUTPUT_FOLDER"

echo "Report name: $REPORT_NAME"

echo "Format: $FORMAT"

echo "Timestamp: $TIMESTAMP"

echo ""

# Check if output folder exists

if [[ ! -d "$OUTPUT_FOLDER" ]]; then

echo "ERROR: Output folder does not exist: $OUTPUT_FOLDER"

exit 1

fi

# Function to check R and required packages

check_r_dependencies() {

echo "Checking R and required packages..."

# Check if R is installed

if ! command -v R &> /dev/null; then

echo "ERROR: R is not installed or not in PATH"

return 1

fi

# Check if required R packages are available

R --slave --no-restore --no-save -e "

required_packages <- c('rmarkdown', 'knitr', 'dplyr', 'ggplot2', 'readr', 'lubridate', 'stringr', 'DT')

missing_packages <- required_packages[!sapply(required_packages, require, character.only=TRUE, quietly=TRUE)]

if(length(missing_packages) > 0) {

cat('ERROR: Missing R packages:', paste(missing_packages, collapse=', '), '\n')

cat('Install with: install.packages(c(', paste(paste0('\"', missing_packages, '\"'), collapse=', '), '))\n')

quit(status=1)

} else {

cat('All required R packages are available\n')

}

" 2>/dev/null

return $?

}

# Function to check LaTeX dependencies for PDF generation

check_latex_dependencies() {

echo "Checking LaTeX dependencies for PDF generation..."

# Check for common LaTeX installations

latex_found=false

if command -v pdflatex &> /dev/null; then

latex_found=true

echo "Found pdflatex: $(which pdflatex)"

elif command -v xelatex &> /dev/null; then

latex_found=true

echo "Found xelatex: $(which xelatex)"

elif command -v lualatex &> /dev/null; then

latex_found=true

echo "Found lualatex: $(which lualatex)"

fi

# Check for TinyTeX (common R installation)

if [[ ! "$latex_found" == true ]]; then

R --slave --no-restore --no-save -e "

if(tinytex:::is_tinytex()) {

cat('Found TinyTeX installation\n')

quit(status=0)

} else {

cat('No LaTeX installation found\n')

quit(status=1)

}

" 2>/dev/null

if [[ $? -eq 0 ]]; then

latex_found=true

fi

fi

if [[ "$latex_found" == true ]]; then

echo "LaTeX installation found - PDF generation supported"

return 0

else

echo "WARNING: No LaTeX installation found"

echo "PDF generation may fail. Consider installing:"

echo " - TinyTeX: R -e \"tinytex::install_tinytex()\""

echo " - Full LaTeX: sudo apt-get install texlive-latex-recommended texlive-latex-extra"

return 1

fi

}

# Run dependency checks

if ! check_r_dependencies; then

echo "ERROR: R dependency check failed"

exit 1

fi

latex_available=true

if [[ "$FORMAT" == "pdf" ]] || [[ "$FORMAT" == "both" ]]; then

if ! check_latex_dependencies; then

if [[ "$FORMAT" == "pdf" ]]; then

echo "ERROR: PDF format requested but LaTeX not available"

exit 1

else

echo "WARNING: Will skip PDF generation and create HTML only"

latex_available=false

FORMAT="html"

fi

fi

fi

echo ""

echo "Creating R Markdown report..."

# Create enhanced R Markdown report with PDF-optimized formatting

cat > "$REPORT_FILE" << 'EOF'

---

title: "DNA FASTQ to BigWig Pipeline Report"

subtitle: "DNAfastqBigWig_human_main_v5_31aug2025.sh Execution Summary"

author: "Automated Pipeline Report"

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

  pdf_document:

    toc: true

    toc_depth: 3

    number_sections: true

    highlight: tango

    df_print: kable

    fig_caption: true

    keep_tex: false

    latex_engine: pdflatex

    includes:

      in_header:

        - \usepackage{booktabs}

        - \usepackage{longtable}

        - \usepackage{array}

        - \usepackage{multirow}

        - \usepackage{wrapfig}

        - \usepackage{float}

        - \usepackage{colortbl}

        - \usepackage{pdflscape}

        - \usepackage{tabu}

        - \usepackage{threeparttable}

        - \usepackage{threeparttablex}

        - \usepackage[normalem]{ulem}

        - \usepackage{makecell}

        - \usepackage{xcolor}

    geometry: "margin=1in"

    fontsize: 11pt

    linestretch: 1.2

---

```{r setup, include=FALSE}

knitr::opts_chunk$set(

  echo = TRUE,

  warning = FALSE,

  message = FALSE,

  fig.width = 10,

  fig.height = 6,

  fig.align = 'center',

  out.width = '100%'

)

# Load required libraries

suppressPackageStartupMessages({

  library(knitr)

  library(dplyr)

  library(ggplot2)

  library(readr)

  library(lubridate)

  library(stringr)

  library(DT)

})

# Set output format specific options

output_format <- knitr::opts_knit$get("rmarkdown.pandoc.to")

if(is.null(output_format)) output_format <- "html"

# Configure table output based on format

if(output_format == "latex") {

  options(knitr.table.format = "latex")

  knitr::opts_chunk$set(fig.pos = 'H')

} else {

  options(knitr.table.format = "html")

}

# Set working directory to output folder

EOF

echo "setwd(\"$OUTPUT_FOLDER\")" >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'EOF'

```

\newpage

# Executive Summary

This report documents the execution of the DNA FASTQ to BigWig processing pipeline (`DNAfastqBigWig_human_main_v5_31aug2025.sh`) which performs comprehensive analysis of paired-end DNA sequencing data from raw FASTQ files to normalized genome coverage tracks.

## Pipeline Overview

The pipeline executes the following major steps:

1. **Quality Control (Initial)** - FastQC analysis of raw FASTQ files

2. **Adapter Trimming** - TrimGalore processing for quality and adapter removal

3. **Quality Control (Post-trimming)** - FastQC analysis of trimmed files

4. **Genome Alignment** - Bowtie2 mapping to hg38 reference genome

5. **Duplicate Removal** - Picard deduplication of aligned reads

6. **Coverage Generation** - Creation of normalized bedGraph and bigWig files

7. **Quality Reporting** - MultiQC summaries at each major step

\newpage

# System Information

## Computing Environment

```{r system_info}

# System information

system_info <- data.frame(

  Parameter = c("Hostname", "User", "Operating System", "R Version", "Report Generated", "Pipeline Script"),

  Value = c(

EOF

echo "    \"$HOSTNAME\"," >> "$REPORT_FILE"

echo "    \"$USER\"," >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'EOF'

    system("uname -a", intern = TRUE),

    R.version.string,

EOF

echo "    \"$TIMESTAMP\"," >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'EOF'

    "DNAfastqBigWig_human_main_v5_31aug2025.sh"

  ),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(system_info, caption = "System and Environment Information",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

} else {

  kable(system_info, caption = "System and Environment Information")

}

```

## Hardware Resources

```{r hardware_info}

# Get hardware information

cpu_info <- system("nproc", intern = TRUE)

memory_info <- system("free -h | grep '^Mem:' | awk '{print $2}'", intern = TRUE)

disk_info <- system("df -h . | tail -1 | awk '{print $2\" (\"$4\" available)\"}'", intern = TRUE)

hardware_info <- data.frame(

  Resource = c("CPU Cores", "Total Memory", "Disk Space (Available)"),

  Specification = c(cpu_info, memory_info, disk_info),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(hardware_info, caption = "Hardware Resources",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

} else {

  kable(hardware_info, caption = "Hardware Resources")

}

```

## Software Versions

```{r software_versions}

# Check software versions

get_version <- function(cmd) {

  tryCatch({

    result <- system(cmd, intern = TRUE)[1]

    if(is.na(result) || result == "") return("Not available")

    return(result)

  }, error = function(e) "Not available")

}

software_versions <- data.frame(

  Software = c("FastQC", "TrimGalore", "Bowtie2", "Samtools", "Picard", "Bedtools", "MultiQC"),

  Version_Command = c(

    "fastqc --version 2>&1 | head -1",

    "trim_galore --version 2>&1 | head -1",

    "bowtie2 --version 2>&1 | head -1",

    "samtools --version 2>&1 | head -1",

    "java -jar /home/micgdu/software/picard.jar MarkDuplicates --version 2>&1 | head -1",

    "bedtools --version 2>&1 | head -1",

    "multiqc --version 2>&1 | head -1"

  ),

  stringsAsFactors = FALSE

)

software_versions$Version <- sapply(software_versions$Version_Command, get_version)

software_versions$Version_Command <- NULL

# Format table based on output

if(output_format == "latex") {

  kable(software_versions, caption = "Software Versions Used",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

} else {

  kable(software_versions, caption = "Software Versions Used")

}

```

\newpage

# Pipeline Execution Analysis

## Script-Level Resource Usage and Timing

```{r script_analysis}

# Define the scripts used in the pipeline with their resource profiles

pipeline_scripts <- data.frame(

  Step = c(

    "1. Initial FastQC",

    "2. TrimGalore",

    "3. Post-trim FastQC",

    "4. Bowtie2 Alignment",

    "5. Picard Deduplication",

    "6. Coverage Generation"

  ),

  Script_Name = c(

    "fastqc_batch_v1_30aug2025.sh",

    "trimgalore_batch_final_v2_30aug2025.sh",

    "fastqc_batch_v1_30aug2025.sh",

    "bowtie2_human_batch_v1_31aug.sh",

    "picard_deduplication_batch_31aug2025_v8.sh",

    "genomecoverage_batch_v1_31aug2025.sh"

  ),

  Lower_Level_Script = c(

    "Built-in FastQC",

    "Built-in TrimGalore",

    "Built-in FastQC",

    "bowtie2_dovetail_pairedEnd_Hsapiens_31Aug25.sh",

    "picard_deduplication_28aug2025.sh",

    "genomeCoverage_DNA_human_26aug2025.sh"

  ),

  Concurrent_Jobs = c("2x parameter", "1x parameter", "2x parameter", "Parameter", "Parameter", "Parameter"),

  Threads_Per_Job = c("10 (fixed)", "8 (fixed)", "10 (fixed)", "15 (Bowtie2)", "1 (Picard)", "10 (samtools)"),

  Memory_Per_Job = c("~2GB", "~4GB", "~2GB", "~8GB", "128GB (Java heap)", "~16GB"),

  Primary_Tool = c("FastQC", "TrimGalore + cutadapt", "FastQC", "Bowtie2 + samtools", "Picard MarkDuplicates", "bedtools + samtools"),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(pipeline_scripts, caption = "Pipeline Scripts and Resource Requirements",

        booktabs = TRUE, longtable = TRUE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "repeat_header"),

                            font_size = 9) %>%

    kableExtra::landscape()

} else {

  kable(pipeline_scripts, caption = "Pipeline Scripts and Resource Requirements")

}

```

## Theoretical Resource Consumption

```{r resource_calculation}

# Calculate theoretical resource usage based on script parameters

calculate_resources <- function(max_jobs_param = 8) {

  resource_calc <- data.frame(

    Step = pipeline_scripts$Step,

    Jobs = c(

      2 * max_jobs_param, # FastQC: 2x parameter

      1 * max_jobs_param, # TrimGalore: 1x parameter

      2 * max_jobs_param, # FastQC: 2x parameter

      1 * max_jobs_param, # Bowtie2: 1x parameter

      1 * max_jobs_param, # Picard: 1x parameter

      1 * max_jobs_param  # Coverage: 1x parameter

    ),

    Threads_Per_Job = c(10, 8, 10, 15, 1, 10),

    Memory_GB_Per_Job = c(2, 4, 2, 8, 128, 16),

    stringsAsFactors = FALSE

  )

  resource_calc$Total_Threads <- resource_calc$Jobs * resource_calc$Threads_Per_Job

  resource_calc$Total_Memory_GB <- resource_calc$Jobs * resource_calc$Memory_GB_Per_Job

  return(resource_calc)

}

# Calculate for default parameter of 8

resource_usage <- calculate_resources(8)

# Format table based on output

if(output_format == "latex") {

  kable(resource_usage, caption = "Theoretical Resource Usage (max\\_jobs parameter = 8)",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"),

                            font_size = 10)

} else {

  kable(resource_usage, caption = "Theoretical Resource Usage (max_jobs parameter = 8)")

}

# Summary statistics

peak_threads <- max(resource_usage$Total_Threads)

peak_memory <- max(resource_usage$Total_Memory_GB)

total_thread_hours <- sum(resource_usage$Total_Threads) * 2 # Assuming 2 hours average per step

```

**Peak Resource Requirements:**

- Maximum concurrent threads: `r peak_threads`

- Maximum memory usage: `r peak_memory` GB

- Total computational thread-hours: ~`r total_thread_hours` (estimated for medium dataset)

\newpage

# File System Analysis

## Directory Structure

```{r directory_structure}

# List the directory structure created by the pipeline

directories <- c(

  "fastQC/",

  "fastQC/fastQC_unTrimmed/",

  "fastQC/fastQC_trimmed/",

  "multiQC/",

  "multiQC/multiQC_unTrimmed/",

  "multiQC/multiQC_trimmed/",

  "multiQC/multiQC_alignments/",

  "multiQC/multiQC_deduplication/",

  "trimmedFastq/",

  "bams/",

  "dedupBams/",

  "bedGraph/",

  "NormBedGraph/",

  "bigwig/"

)

dir_info <- data.frame(

  Directory = directories,

  Purpose = c(

    "FastQC reports root directory",

    "FastQC reports for raw FASTQ files",

    "FastQC reports for trimmed FASTQ files",

    "MultiQC reports root directory",

    "MultiQC summary for raw FASTQ analysis",

    "MultiQC summary for trimmed FASTQ analysis",

    "MultiQC summary for alignment statistics",

    "MultiQC summary for deduplication statistics",

    "Trimmed FASTQ files (*_val_*.fq.gz)",

    "Aligned BAM files (*_sorted_stChr.bam)",

    "Deduplicated BAM files (*_dedup.bam)",

    "Raw bedGraph coverage files",

    "Normalized bedGraph coverage files",

    "BigWig coverage files (*_Snorm.bw)"

  ),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(dir_info, caption = "Pipeline Output Directory Structure",

        booktabs = TRUE, longtable = TRUE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "repeat_header"),

                            font_size = 9)

} else {

  kable(dir_info, caption = "Pipeline Output Directory Structure")

}

```

## File Size Analysis

```{r file_analysis}

# Function to get file sizes and counts

get_directory_info <- function(dir_path) {

  if (!dir.exists(dir_path)) {

    return(data.frame(

      Directory = basename(dir_path),

      Files = 0,

      Total_Size_GB = 0,

      Avg_Size_MB = 0,

      stringsAsFactors = FALSE

    ))

  }

  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)

  if (length(files) == 0) {

    return(data.frame(

      Directory = basename(dir_path),

      Files = 0,

      Total_Size_GB = 0,

      Avg_Size_MB = 0,

      stringsAsFactors = FALSE

    ))

  }

  file_info <- file.info(files)

  total_size_bytes <- sum(file_info$size, na.rm = TRUE)

  total_size_gb <- total_size_bytes / (1024^3)

  avg_size_mb <- (total_size_bytes / length(files)) / (1024^2)

  data.frame(

    Directory = basename(dir_path),

    Files = length(files),

    Total_Size_GB = round(total_size_gb, 2),

    Avg_Size_MB = round(avg_size_mb, 2),

    stringsAsFactors = FALSE

  )

}

# Analyze each directory

dir_analysis <- do.call(rbind, lapply(directories, function(d) {

  get_directory_info(file.path(getwd(), d))

}))

# Add totals row

totals <- data.frame(

  Directory = "**TOTAL**",

  Files = sum(dir_analysis$Files),

  Total_Size_GB = sum(dir_analysis$Total_Size_GB),

  Avg_Size_MB = round(mean(dir_analysis$Avg_Size_MB), 2),

  stringsAsFactors = FALSE

)

dir_analysis_with_totals <- rbind(dir_analysis, totals)

# Format table based on output

if(output_format == "latex") {

  kable(dir_analysis_with_totals, caption = "File System Usage Analysis",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

} else {

  kable(dir_analysis_with_totals, caption = "File System Usage Analysis")

}

```

\newpage

# Log and Report Analysis

## Batch Processing Logs

```{r log_analysis}

# Function to find and analyze log directories

find_log_dirs <- function(pattern) {

  log_dirs <- list.dirs(".", recursive = TRUE, full.names = FALSE)

  log_dirs[grepl(pattern, log_dirs)]

}

# Find all batch processing log directories

log_types <- c(

  "fastqc_batch_logs",

  "trimgalore_batch_logs",

  "bowtie2_batch_logs",

  "picard_batch_logs",

  "genomecov_batch_logs"

)

log_summary <- data.frame(

  Log_Type = c(

    "FastQC Batch Logs",

    "TrimGalore Batch Logs",

    "Bowtie2 Batch Logs",

    "Picard Batch Logs",

    "GenomeCov Batch Logs"

  ),

  Pattern = paste0(log_types, "_*"),

  Purpose = c(

    "Individual FastQC job logs, timing, and errors",

    "Individual TrimGalore job logs and statistics",

    "Individual Bowtie2 alignment logs and metrics",

    "Individual Picard deduplication logs and metrics",

    "Individual genome coverage generation logs"

  ),

  Contains = c(

    "Job logs, error logs, main batch log, PID file",

    "Job logs, error logs, main batch log, PID file",

    "Job logs, error logs, main batch log, PID file",

    "Job logs, error logs, main batch log, PID file",

    "Job logs, error logs, main batch log, PID file"

  ),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(log_summary, caption = "Batch Processing Log Structure",

        booktabs = TRUE, longtable = TRUE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "repeat_header"),

                            font_size = 8)

} else {

  kable(log_summary, caption = "Batch Processing Log Structure")

}

# Try to find actual log directories

actual_logs <- character(0)

for (pattern in log_types) {

  found <- find_log_dirs(pattern)

  if (length(found) > 0) {

    actual_logs <- c(actual_logs, found)

  }

}

```

`r if(length(actual_logs) > 0) paste("**Found Log Directories:**", paste("- ", actual_logs, collapse = "\n"), collapse = "\n") else "**Note:** No batch processing log directories found in current location. Log directories are created in the respective output folders during pipeline execution."`

\newpage

# Performance Metrics and Summary

## Pipeline Efficiency Features

```{r efficiency_features}

efficiency_features <- data.frame(

  Feature = c(

    "Continuous Job Replacement",

    "Signal Immunity",

    "Comprehensive Logging",

    "Smart Skip Logic",

    "PID Tracking",

    "Progress Monitoring",

    "Resource Optimization"

  ),

  Description = c(

    "New jobs start immediately when others finish, maximizing CPU utilization",

    "Scripts immune to SIGHUP, SIGINT, SIGTERM - safe for remote execution",

    "Detailed logs for each job, main process, and error tracking",

    "Automatically skips files that have already been processed",

    "Tracks process IDs for job management and cleanup",

    "Regular status updates every 5 minutes during execution",

    "Parameterized job counts and thread usage for different system capacities"

  ),

  Benefit = c(

    "Reduced total processing time",

    "Prevents data loss from disconnections",

    "Easy troubleshooting and monitoring",

    "Enables pipeline restart/resume",

    "Clean termination and resource cleanup",

    "Real-time execution monitoring",

    "Optimal resource utilization"

  ),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(efficiency_features, caption = "Pipeline Efficiency and Reliability Features",

        booktabs = TRUE, longtable = TRUE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "repeat_header"),

                            font_size = 8)

} else {

  kable(efficiency_features, caption = "Pipeline Efficiency and Reliability Features")

}

```

## Final Resource Requirements Summary

```{r final_summary}

# Create final summary table

final_summary <- data.frame(

  Metric = c(

    "Total Pipeline Steps",

    "Core Processing Scripts",

    "Maximum Concurrent Jobs (default)",

    "Peak Thread Usage (default)",

    "Peak Memory Usage (default)",

    "Primary Output File Types"

  ),

  Value = c(

    "7 major steps",

    "8 specialized scripts",

    "16 jobs (2x FastQC steps)",

    paste0(peak_threads, " threads (", peak_threads/10, " FastQC jobs × 10 threads)"),

    paste0(peak_memory, " GB (", peak_memory/128, " Picard jobs × 128GB)"),

    "BigWig, BAM, FastQC reports, MultiQC summaries"

  ),

  stringsAsFactors = FALSE

)

# Format table based on output

if(output_format == "latex") {

  kable(final_summary, caption = "Pipeline Resource Requirements Summary",

        booktabs = TRUE, longtable = FALSE) %>%

    kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

} else {

  kable(final_summary, caption = "Pipeline Resource Requirements Summary")

}

```

\newpage

# Complete File Inventory

This section provides a comprehensive inventory of all files generated during the pipeline execution, organized by directory and file type.

```{r complete_file_inventory, results='asis'}

# Function to format file sizes

format_size <- function(bytes) {

  if (is.na(bytes) || bytes == 0) return("0 B")

  units <- c("B", "KB", "MB", "GB", "TB")

  unit_index <- min(floor(log(bytes, 1024)) + 1, length(units))

  size_value <- bytes / (1024^(unit_index - 1))

  return(paste0(round(size_value, 2), " ", units[unit_index]))

}

# Function to get detailed file information with sizes (non-recursive for main files only)

get_detailed_file_info <- function(dir_path, base_path = ".") {

  full_path <- file.path(base_path, dir_path)

  

  if (!dir.exists(full_path)) {

    return(list(

      path = dir_path,

      folder_size = "0 B",

      files = data.frame(

        Filename = character(0),

        Size = character(0),

        Type = character(0),

        stringsAsFactors = FALSE

      )

    ))

  }

  # Get folder size (including subdirectories)

  folder_size_cmd <- paste0("du -sb '", full_path, "' 2>/dev/null | cut -f1")

  folder_size_result <- system(folder_size_cmd, intern = TRUE)

  folder_size_bytes <- if(length(folder_size_result) > 0) as.numeric(folder_size_result) else 0

  folder_size <- format_size(folder_size_bytes)

  

  # Get files in directory (NON-recursive to avoid log contamination)

  files <- list.files(full_path, full.names = TRUE, recursive = FALSE)

  files <- files[!dir.exists(files)] # Only files, not subdirectories

  

  if (length(files) == 0) {

    return(list(

      path = dir_path,

      folder_size = folder_size,

      files = data.frame(

        Filename = character(0),

        Size = character(0),

        Type = character(0),

        stringsAsFactors = FALSE

      )

    ))

  }

  # Get file info

  file_info <- file.info(files)

  filenames <- basename(files)

  sizes_bytes <- file_info$size

  sizes_formatted <- sapply(sizes_bytes, format_size)

  

  # Determine file types based on extensions

  get_file_type <- function(filename) {

    ext <- tools::file_ext(tolower(filename))

    if (ext == "") return("No extension")

    switch(ext,

      "fq.gz" = "FASTQ (compressed)",

      "fastq.gz" = "FASTQ (compressed)",

      "bam" = "BAM alignment",

      "bai" = "BAM index",

      "sam" = "SAM alignment",

      "bw" = "BigWig coverage",

      "bigwig" = "BigWig coverage",

      "bedgraph" = "BedGraph coverage",

      "gz" = "Compressed file",

      "html" = "HTML report",

      "zip" = "ZIP archive",

      "txt" = "Text file",

      "log" = "Log file",

      "pid" = "Process ID file",

      "json" = "JSON data",

      "csv" = "CSV data",

      "tsv" = "TSV data",

      "sizes" = "Chromosome sizes",

      paste0(toupper(ext), " file")

    )

  }

  file_types <- sapply(filenames, get_file_type)

  

  # Create data frame

  files_df <- data.frame(

    Filename = filenames,

    Size = sizes_formatted,

    Type = file_types,

    stringsAsFactors = FALSE

  )

  

  # Sort by size (convert back to numeric for sorting)

  size_numeric <- sapply(sizes_bytes, function(x) if(is.na(x)) 0 else x)

  files_df <- files_df[order(size_numeric, decreasing = TRUE), ]

  

  return(list(

    path = dir_path,

    folder_size = folder_size,

    files = files_df

  ))

}

# Define main pipeline directories to analyze (exclude log directories)

main_directories <- c(

  "fastQC/fastQC_unTrimmed",

  "fastQC/fastQC_trimmed", 

  "multiQC/multiQC_unTrimmed",

  "multiQC/multiQC_trimmed",

  "multiQC/multiQC_alignments",

  "multiQC/multiQC_deduplication",

  "trimmedFastq",

  "bams",

  "dedupBams",

  "bigwig"

)

# Get file information for all directories

all_file_info <- lapply(main_directories, get_detailed_file_info)

names(all_file_info) <- main_directories

# Generate output for each directory

for (i in seq_along(all_file_info)) {

  dir_name <- names(all_file_info)[i]

  dir_info <- all_file_info[[i]]

  

  cat("\n## Directory:", dir_name, "\n\n")

  cat("**Folder Path:** `", file.path(getwd(), dir_name), "`  \n")

  cat("**Folder Size:** ", dir_info$folder_size, "  \n\n")

  

  if (nrow(dir_info$files) > 0) {

    if(output_format == "latex") {

      print(kable(dir_info$files, caption = paste("Files in", dir_name),

                  booktabs = TRUE, longtable = FALSE) %>%

        kableExtra::kable_styling(latex_options = c("striped", "hold_position"),

                                font_size = 8))

    } else {

      print(kable(dir_info$files, caption = paste("Files in", dir_name)))

    }

  } else {

    cat("*No files found in this directory.*\n")

  }

  

  cat("\n")

}

```

## File Type Summary

```{r file_type_summary}

# Aggregate file information by type (only from main pipeline directories)

type_summary <- list()

for (dir_name in names(all_file_info)) {

  dir_info <- all_file_info[[dir_name]]

  if (nrow(dir_info$files) > 0) {

    for (i in 1:nrow(dir_info$files)) {

      file_type <- dir_info$files$Type[i]

      file_size_str <- dir_info$files$Size[i]

      

      # Convert size string back to bytes for aggregation

      size_parts <- strsplit(file_size_str, " ")[[1]]

      size_value <- as.numeric(size_parts[1])

      size_unit <- size_parts[2]

      size_bytes <- switch(size_unit,

        "B" = size_value,

        "KB" = size_value * 1024,

        "MB" = size_value * 1024^2,

        "GB" = size_value * 1024^3,

        "TB" = size_value * 1024^4,

        size_value

      )

      

      if (file_type %in% names(type_summary)) {

        type_summary[[file_type]]$count <- type_summary[[file_type]]$count + 1

        type_summary[[file_type]]$total_bytes <- type_summary[[file_type]]$total_bytes + size_bytes

        type_summary[[file_type]]$locations <- unique(c(type_summary[[file_type]]$locations, dir_name))

      } else {

        type_summary[[file_type]] <- list(

          count = 1,

          total_bytes = size_bytes,

          locations = dir_name

        )

      }

    }

  }

}

# Convert to data frame

if (length(type_summary) > 0) {

  type_summary_df <- data.frame(

    File_Type = names(type_summary),

    Count = sapply(type_summary, function(x) x$count),

    Total_Size = sapply(type_summary, function(x) format_size(x$total_bytes)),

    Locations = sapply(type_summary, function(x) paste(x$locations, collapse = ", ")),

    stringsAsFactors = FALSE

  )

  

  # Sort by total size (descending)

  size_bytes_for_sort <- sapply(type_summary, function(x) x$total_bytes)

  type_summary_df <- type_summary_df[order(size_bytes_for_sort, decreasing = TRUE), ]

  

  if(output_format == "latex") {

    kable(type_summary_df, caption = "File Type Summary - Pipeline Files Only",

          booktabs = TRUE, longtable = TRUE) %>%

      kableExtra::kable_styling(latex_options = c("striped", "repeat_header"),

                              font_size = 8)

  } else {

    kable(type_summary_df, caption = "File Type Summary - Pipeline Files Only")

  }

} else {

  cat("*No files found for type summary.*\n")

}

```

## Storage Summary

```{r storage_summary}

# Calculate total storage used by main pipeline directories only

total_bytes <- 0

directory_sizes <- data.frame(

  Directory = character(0),

  Size_GB = numeric(0),

  Percentage = numeric(0),

  stringsAsFactors = FALSE

)

for (dir_name in names(all_file_info)) {

  dir_info <- all_file_info[[dir_name]]

  

  # Convert folder size to bytes

  size_str <- dir_info$folder_size

  if (size_str != "0 B") {

    size_parts <- strsplit(size_str, " ")[[1]]

    size_value <- as.numeric(size_parts[1])

    size_unit <- size_parts[2]

    size_bytes <- switch(size_unit,

      "B" = size_value,

      "KB" = size_value * 1024,

      "MB" = size_value * 1024^2,

      "GB" = size_value * 1024^3,

      "TB" = size_value * 1024^4,

      size_value

    )

    

    total_bytes <- total_bytes + size_bytes

    directory_sizes <- rbind(directory_sizes, data.frame(

      Directory = dir_name,

      Size_GB = round(size_bytes / (1024^3), 3),

      Percentage = 0, # Will calculate after total

      stringsAsFactors = FALSE

    ))

  }

}

# Calculate percentages

if (nrow(directory_sizes) > 0) {

  directory_sizes$Percentage <- round((directory_sizes$Size_GB * 1024^3 / total_bytes) * 100, 1)

  directory_sizes <- directory_sizes[order(directory_sizes$Size_GB, decreasing = TRUE), ]

  

  if(output_format == "latex") {

    kable(directory_sizes, caption = "Storage Usage by Directory - Pipeline Files Only",

          booktabs = TRUE, longtable = FALSE) %>%

      kableExtra::kable_styling(latex_options = c("striped", "hold_position"))

  } else {

    kable(directory_sizes, caption = "Storage Usage by Directory - Pipeline Files Only")

  }

  

  cat("\n**Total Pipeline Storage Usage:** ", format_size(total_bytes), "\n")

} else {

  cat("*No directory size information available.*\n")

}

```

## Conclusion

This automated DNA sequencing analysis pipeline provides a comprehensive workflow from raw FASTQ files to normalized genome coverage tracks. Key achievements include:

- **Automated Quality Control:** Multi-stage QC with FastQC and MultiQC reporting

- **Robust Processing:** Signal-immune batch processing with comprehensive logging

- **Scalable Architecture:** Parameterized resource allocation for different system capacities

- **Production Ready:** Skip logic enables restart/resume capabilities

### Recommendations

- **Minimum system:** 32 cores, 256GB RAM, 1TB+ storage

- **Optimal system:** 64+ cores, 512GB+ RAM, fast SSD storage

- **Monitor disk space:** Intermediate files can be 3-5x input size

- **Use parameter tuning** based on available resources

### Generated Files and Locations

All output files are organized in the specified output directory with the following structure:

- **Analysis Results:** BAM files, coverage tracks, quality reports

- **Logs and Metrics:** Comprehensive logging for troubleshooting and monitoring

- **Intermediate Files:** Trimmed FASTQ files, alignment statistics, duplication metrics

---

*Report generated automatically by the DNAfastqBigWig pipeline reporting system.*

*For questions about this pipeline, consult the individual script documentation and log files.*

EOF

echo "R Markdown report created: $REPORT_FILE"

# Function to render the report

render_report() {

local output_format="$1"

local output_file="${OUTPUT_FOLDER}/${REPORT_NAME}.${output_format}"

echo "Rendering $output_format report..."

R --slave --no-restore --no-save -e "

tryCatch({

rmarkdown::render(

input = '$REPORT_FILE',

output_format = '${output_format}_document',

output_file = '$output_file',

quiet = FALSE

)

cat('SUCCESS: $output_format report generated: $output_file\n')

}, error = function(e) {

cat('ERROR generating $output_format report:\n')

cat(conditionMessage(e), '\n')

quit(status=1)

})

" 2>&1

return $?

}

# Generate reports based on requested format

echo ""

echo "Generating report(s)..."

html_success=false

pdf_success=false

# Generate HTML report

if [[ "$FORMAT" == "html" ]] || [[ "$FORMAT" == "both" ]]; then

if render_report "html"; then

html_success=true

echo "✅ HTML report generated successfully"

else

echo "❌ HTML report generation failed"

fi

fi

# Generate PDF report

if [[ "$FORMAT" == "pdf" ]] || ([[ "$FORMAT" == "both" ]] && [[ "$latex_available" == true ]]); then

if render_report "pdf"; then

pdf_success=true

echo "✅ PDF report generated successfully"

else

echo "❌ PDF report generation failed"

if [[ "$FORMAT" == "pdf" ]]; then

echo "ERROR: PDF generation failed and no fallback requested"

exit 1

else

echo "WARNING: PDF generation failed, but HTML version was created"

fi

fi

fi

echo ""

echo "=== Report Generation Summary ==="

echo "R Markdown source: $REPORT_FILE"

if [[ "$html_success" == true ]]; then

echo "HTML report: ${OUTPUT_FOLDER}/${REPORT_NAME}.html"

fi

if [[ "$pdf_success" == true ]]; then

echo "PDF report: ${OUTPUT_FOLDER}/${REPORT_NAME}.pdf"

fi

if [[ "$html_success" == false ]] && [[ "$pdf_success" == false ]]; then

echo "ERROR: No reports were generated successfully"

exit 1

fi

echo ""

echo "✅ Report generation completed successfully!"

# Clean up if requested (optional)

if [[ "$CLEANUP_RMD" == "true" ]]; then

rm "$REPORT_FILE"

echo "Cleaned up R Markdown source file"

fi