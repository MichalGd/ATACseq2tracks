#!/bin/bash
# This script mapps reads stored in fastq files with Bowtie2 and generates 1.) bam alignment files 2.) coverage files 
# stored as read-depth normalised and unnormalised bedGrpaph files and bigWig files. It anlyses property and quality of 
# fastq files with FastQC (REF) and summarises this information with multiQC (REF). AFter first FastQC check adapter and
# low quqlity sequences are removed with TrimmGalore (REF) from fastq files and FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
# and multiQC (https://github.com/MultiQC/MultiQC) are used again to asses
# changes. Trimmed fatq files are mapped with Bowtie2 with "dovetail" setting enabling work with overlapping paired-end 
# reads. Mapped reads are stored in bam files. Duplicated alignments are removed with picard (REF). Deduplicated bam files 
# are used to generate raw and read-depth normalised  genoem coverages stored as bedGraph files with bedtools. Lastly 
# normalised bedGraph files serve to generate bigWig files used later for data browsing. 
#
# User has to provide path to the folder with fastq files with reads belonging to single species - currently either mouse 
# or human [input folder] the folder where all results, los and reports will be stored in teh set of subfolder. 
# Additionally user has to define number of samples jobs processed in  parallel - suggested number for our system: 8     
#
# script assumes that one sample is represented by 1 pair of fastq files stored in seperate files in one folder anmed:
# path/to/folder/*_1.fq.gz (or fastq}.gz)
# path/to/folder/*_2.fq.gz (or fastq}.gz)
# It has been tested initially on 150 bp paire-end reads geberated by Illumina NovaSeqX Plus
# usage ./DNAfastqBigWig_human_main_v5_31aug2025.sh [inputFolder] [outputFolder] [max jobs] [species]
# [inputFolder] - folder with fastq files; seperately human and mouse; distinct scripts for each
# [outputFolder] - for all script outputs
# [max jobs] - maximal number of jobs to be run in parallel
# [species] - mouse or human
  
# ### Scripts used in the workflow
# /home/micgdu/workflows/RNAseq/scripts/DNAfastqBigWig_human_main_v5_31aug2025.sh  # master script
# /home/micgdu/workflows/RNAseq/scripts/fastqc_batch_v1_30aug2025.sh
# /home/micgdu/workflows/RNAseq/scripts/trimgalore_batch_final_v2_30aug2025.sh
# /home/micgdu/workflows/RNAseq/scripts/bowtie2_human_batch_v1_31aug.sh
# /home/micgdu/workflows/RNAseq/scripts/bowtie2_dovetail_pairedEnd_Hsapiens_31Aug25.sh
# /home/micgdu/workflows/RNAseq/scripts/picard_deduplication_batch_31aug2025_v8.sh
# /home/micgdu/workflows/RNAseq/scripts/picard_deduplication_28aug2025.sh
# /home/micgdu/workflows/RNAseq/scripts/genomecoverage_batch_v1_31aug2025.sh
# /home/micgdu/workflows/RNAseq/scripts/genomeCoverage_DNA_human_26aug2025.sh
# /home/micgdu/workflows/RNAseq/scripts/generate_multiqc_unified_report_v2.sh

# validate of $4 variable:
species="${4:-}"

if [[ "$species" != "human" && "$species" != "mouse" ]]; then
  echo "wrong species only mouse and human available" >&2
  exit 1
fi

### 1. Craeting folders:
cd $2  
# $2/run3_20aug2025/scripts/ #in the future lower-level scripts allong the amjor one should be copied into this folder and mede executable with "chmod a+x"  
mkdir $2/reports/
#mkdir $2/fastQC/
#mkdir $2/fastQC/fastQC_unTrimmed/
#mkdir $2/fastQC/fastQC_trimmed/
#mkdir $2/multiQC/
#mkdir $2/multiQC/multiQC_unTrimmed
#mkdir $2/multiQC/multiQC_trimmed
mkdir $2/multiQC/multiQC_alignments
mkdir $2/multiQC/multiQC_deduplication
mkdir $2/trimmedFastq
mkdir $2/bams
mkdir $2/dedupBams
mkdir $2/bedGraph
mkdir $2/NormBedGraph
mkdir $2/bigwig
  
### 2. Quality control of unprocessed fastq files :
# Analysis with fastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# Enhanced FastQC Batch Processing Script - Continuous Job Replacement for Maximum Efficiency
# Usage: ./fastqc_batch_continuous.sh <input_folder> <output_folder> [max_jobs]
# Example: ./fastqc_batch_continuous.sh /path/to/fastqs /path/to/output 8
  
# /home/micgdu/workflows/RNAseq/scripts/fastqc_batch.1.0.sh $1 $2/fastQC/fastQC_unTrimmed/ $3
 
# Summarizing fastQC results with multiQC
# source /home/micgdu/myenv/bin/activate
# multiqc   $2/fastQC/ -n multiQC_unTrimmed -o $2/multiQC/multiQC_unTrimmed/ --data-format tsv --export
 
### 3. trimming fastq files with TrmmGalore based on quality & sequence - adpater removal
# Enhanced TrimGalore Batch Processing Script - Continuous Job Replacement for Maximum Efficiency
# Usage: ./trimgalore_batch_continuous.sh <input_folder> <output_folder> [max_jobs]
# Example: ./trimgalore_batch_continuous.sh /path/to/fastqs /path/to/output 8
 
# /home/micgdu/workflows/RNAseq/scripts/trimgalore_batch.1.0.sh $1 $2/trimmedFastq $3
 
### 4. Quality control of trimmed fastq files :
 
# /home/micgdu/workflows/RNAseq/scripts/fastqc_batch.1.0.sh $2/trimmedFastq $2/fastQC/fastQC_trimmed/ $3
 
# Summarizing fastQC results with multiQC
source /home/micgdu/myenv/bin/activate
#multiqc   $2/fastQC/fastQC_trimmed/ -n multiQC_trimmed -o $2/multiQC/multiQC_trimmed --data-format tsv --export
 
### 5. Alignment of trimmed fastq files with Bowtie2 to hg38 genome
# Enhanced Bowtie2 Batch Processing Script - Continuous Job Replacement for Maximum Efficiency
# Usage: ./bowtie2_batch_continuous.sh <index> <input_folder> <output_folder> [max_jobs]
# Example: ./bowtie2_batch_continuous.sh /path/to/index /path/to/fastqs /path/to/output 8
 
if [[ "$species" == "human" ]]; then
  /home/micgdu/workflows/RNAseq/scripts/bowtie2_human_batch.2.1.sh \
    /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38 \
    "$2/trimmedFastq" \
    "$2/bams" \
    "$3"
else  # mouse
  /home/micgdu/workflows/RNAseq/scripts/bowtie2_mouse_batch.2.1.sh \
    /home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39 \
    "$2/trimmedFastq" \
    "$2/bams" \
    "$3"
fi

# Summarizing fastQC results with multiQC
source /home/micgdu/myenv/bin/activate
multiqc   $2/bams/ -n multiQC_aligments -o $2/multiQC/multiQC_alignments --data-format tsv --export
 
### 6. Removal of read duplicates with picard

# modifying bam file header - adding @RG tag
cd "$2/bams/" || exit 1

samples=( *stChr.bam )
if [[ ! -e "${samples[0]}" ]]; then
  echo "ERROR: No *stChr.bam files found in $2/bams/" >&2
  exit 1
fi

rg_pids=()
for sample in "${samples[@]}"; do
  samtools addreplacerg \
    -r "@RG\tID:RG1\tSM:third_run_NovaSeqXplus_aug2025\tPL:Illumina\tLB:Library.fa" \
    -o "${sample/.bam/H.bam}" \
    "$sample" &
  rg_pids+=($!)
done

rg_failed=0
for pid in "${rg_pids[@]}"; do
  if ! wait "$pid"; then
    rg_failed=1
  fi
done

if [[ "$rg_failed" -ne 0 ]]; then
  echo "ERROR: samtools addreplacerg failed for at least one BAM; stopping before Picard." >&2
  exit 1
fi

for sample in "${samples[@]}"; do
  rg_bam="${sample/.bam/H.bam}"
  if [[ ! -s "$rg_bam" ]] || ! samtools quickcheck "$rg_bam" 2>/dev/null; then
    echo "ERROR: Invalid or truncated BAM after addreplacerg: $rg_bam" >&2
    exit 1
  fi
done

# removing the duplicates with picard
# Usage: ./picard_deduplication_batch_28aug2025_v7.sh <input_folder> <output_folder> [max_jobs]
# Example: ./picard_deduplication_batch_28aug2025_v7.sh /path/to/bams /path/to/output 8
/home/micgdu/workflows/RNAseq/scripts/picard_deduplication_batch.2.1.sh "$2/bams/" "$2/dedupBams" "$3"

# Summarizing fastQC results with multiQC
source /home/micgdu/myenv/bin/activate
multiqc   $2/dedupBams -n multiQC_deduplication -o $2/multiQC/multiQC_deduplication --data-format tsv --export

### 7. Generation of coverages, rpm normalization and generation of bedGraph & bigWig files
# Usage: ./genomecoverage_batch_v1_31aug2025.sh <input_folder> <genome> [max_jobs]
# Example: ./genomecoverage_batch_continuous.sh /path/to/bams hg38 8

if [[ "$species" == "human" ]]; then
  /home/micgdu/workflows/RNAseq/scripts/genomecoverage_batch.1.0.sh $2/dedupBams hg38 $3 $4
else  # mouse
  /home/micgdu/workflows/RNAseq/scripts/genomecoverage_batch.1.0.sh $2/dedupBams mm39 $3 $4
fi

# move outputs first
mv "$2"/dedupBams/*.bw "$2"/bigwig
mv "$2"/dedupBams/*.bedGraph.gz "$2"/bedGraph
mv "$2"/dedupBams/*Snorm*.bedGraph.gz "$2"/NormBedGraph

### 8. Creating generic UCSC bigwig track lines
# UCSC BigWig Track Generator Script
# Creates UCSC track lines for bigWig files in the specified folder
# Usage: ./create_ucsc_tracks.sh [output_folder] [remote_url_base]
# Example: ./create_ucsc_tracks.sh /dysk2/results https://myserver.com/data

/home/micgdu/workflows/RNAseq/scripts/create_ucsc_tracks.1.0.sh $2  http://your-server.com/data

# removing unnecessary large files:
rm "$2"/bams/*.bam "$2"/bams/*.bai
rm $2/trimmedFastq/*.gz

### 8. Generating analysis report
## Usage with format options
#./generate_pipeline_report_pdf.sh <output_folder> <report_name> [format]
## Examples:
#./generate_pipeline_report_pdf.sh /dysk2/results report_$(date +%Y%m%d) both    # HTML + PDF
#./generate_pipeline_report_pdf.sh /dysk2/results report_$(date +%Y%m%d) pdf     # PDF only  
#./generate_pipeline_report_pdf.sh /dysk2/results report_$(date +%Y%m%d) html    # HTML only

/home/micgdu/workflows/RNAseq/scripts/generate_pipeline_report.1.0.sh $2 "pipeline_report_$(date +%Y%m%d)" html

### 9. Generating report with summary of multiqc analysis of fastq files, trimming, aligmnet and deduplication
# Enhanced MultiQC Report Aggregator for DNAfastqBigWig Pipeline
# This script generates a comprehensive R Markdown HTML report that extracts and displays
# key plots and tables from all MultiQC reports generated by DNAfastqBigWig_human_main_v7_1sept2025.sh
# Usage: ./generate_multiqc_summary_report.sh <output_folder> <report_name>
# Example: ./generate_multiqc_summary_report.sh /dysk2/results multiqc_summary_$(date +%Y%m%d)
# Expected folder structure
# /your/output/folder/
# ├── multiQC/multiQC_unTrimmed/
# ├── multiQC/multiQC_trimmed/
# ├── multiQC/multiQC_alignments/
# └── multiQC/multiQC_deduplication/

###################################
# Enhanced MultiQC Report Aggregator for DNAfastqBigWig Pipeline v8
# This script generates a comprehensive R Markdown HTML report that extracts and displays
# key plots and tables from all MultiQC reports generated by DNAfastqBigWig_human_main_v8_1sept2025_full.sh

# Usage: ./generate_multiqc_summary_report_v8_complete.sh <input_folder> <output_folder>
# Example: ./generate_multiqc_summary_report_v8_complete.sh /path/to/pipeline/results /path/to/reports

/home/micgdu/workflows/RNAseq/scripts/generate_multiqc_unified_report.1.0.sh $2 $2/reports/multiqc_summary_$(date +%Y%m%d) selfcontained

# folder structure:

# inputFolder/
# ├── trimmedFastq/
# ├── bams/
# ├── dedupBams/
# ├── bedGraph/
# ├── NormBedGraph/
# ├── bigwig/
# ├── fastQC/
# │   ├── fastQC_unTrimmed/
# │   └── fastQC_trimmed/
# ├──── multiQC/
# │   ├── multiQC_unTrimmed/
# │   ├── multiQC_trimmed/
# │   ├── multiQC_alignments/
# │   └── multiQC_deduplication/
# └── reports/