#!/bin/bash

# Enhanced FastQC Batch Processing Script - Continuous Job Replacement for Maximum Efficiency
# Usage: ./fastqc_batch_continuous.sh <input_folder> <output_folder> [max_jobs]
# Example: ./fastqc_batch_continuous.sh /path/to/fastqs /path/to/output 8

# Function to display usage information
usage() {
    echo "Usage: $0 <input_folder> <output_folder> [max_jobs]"
    echo ""
    echo "Parameters:"
    echo "  input_folder   - Directory containing FASTQ files (*.fq.gz)"
    echo "  output_folder  - Directory where FastQC reports will be saved"
    echo "  max_jobs       - Jobs parameter (actual concurrent jobs will be 2x this, default: 8, optional)"
    echo ""
    echo "Examples:"
    echo "  $0 /dysk2/fastqs /dysk2/fastqc_output 8     # Will run 16 concurrent jobs with 10 threads each"
    echo "  $0 /dysk2/fastqs /dysk2/fastqc_output     # Uses default 8 (16 actual jobs)"
    echo ""
    echo "Note: Jobs are started continuously - as soon as one finishes, the next begins!"
    echo ""
    exit 1
}

# Check if required parameters are provided
if [[ $# -lt 2 ]]; then
    echo "ERROR: Missing required parameters"
    echo ""
    usage
fi

# Parse command line arguments
inFolder="$1"
outFolder="$2"
max_jobs_param="${3:-8}"  # Default to 8 if not provided
max_jobs=$((2 * max_jobs_param))  # Double the provided parameter

# Validate input parameters
if [[ ! -d "$inFolder" ]]; then
    echo "ERROR: Input folder does not exist: $inFolder"
    exit 1
fi

# Create output folder if it doesn't exist
if [[ ! -d "$outFolder" ]]; then
    echo "Creating output folder: $outFolder"
    mkdir -p "$outFolder"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to create output folder: $outFolder"
        exit 1
    fi
fi

# Validate max_jobs parameter
if ! [[ "$max_jobs_param" =~ ^[0-9]+$ ]] || [[ "$max_jobs_param" -lt 1 ]]; then
    echo "ERROR: max_jobs parameter must be a positive integer, got: $max_jobs_param"
    exit 1
fi

# Set fixed number of threads for FastQC
threads=10

# Create FastQC output directory
mkdir -p "$outFolder/fastQC"

# Set up logging in output folder
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_SUBFOLDER="$outFolder/fastqc_batch_logs_${TIMESTAMP}"
MAIN_LOG="$LOG_SUBFOLDER/fastqc_batch_${TIMESTAMP}.log"
ERROR_LOG="$LOG_SUBFOLDER/fastqc_batch_errors_${TIMESTAMP}.log"
JOBS_LOG_DIR="$LOG_SUBFOLDER/individual_jobs"
MASTER_PID_FILE="$LOG_SUBFOLDER/fastqc_batch_${TIMESTAMP}.pid"

# Create log directories
mkdir -p "$LOG_SUBFOLDER"
mkdir -p "$JOBS_LOG_DIR"

# Make script immune to hangup signals
exec > >(tee -a "$MAIN_LOG") 2>&1

# Create PID file for monitoring
echo $$ > "$MASTER_PID_FILE"

# Logging functions
log_main() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [PID:$$] $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [PID:$$] ERROR: $1" | tee -a "$ERROR_LOG"
}

# Set up signal handling for graceful shutdown
cleanup() {
    exit_code=$?
    log_main "Script terminated (exit code: $exit_code)"
    # Kill remaining jobs
    if [[ ${#job_pids[@]} -gt 0 ]]; then
        log_main "Terminating ${#job_pids[@]} remaining jobs..."
        for pid in "${job_pids[@]}"; do
            [[ -n "$pid" ]] && kill -TERM "$pid" 2>/dev/null
        done
        sleep 2
        for pid in "${job_pids[@]}"; do
            [[ -n "$pid" ]] && kill -KILL "$pid" 2>/dev/null
        done
    fi
    rm -f "$MASTER_PID_FILE" 2>/dev/null
    exit $exit_code
}

trap cleanup SIGHUP SIGINT SIGTERM

# Change to input directory
cd "$inFolder" || {
    echo "ERROR: Cannot access input folder: $inFolder"
    exit 1
}

log_main "=== Starting FastQC Batch Processing (Continuous Job Replacement) ==="
log_main "Input: $inFolder"
log_main "Output: $outFolder"
log_main "Jobs parameter: $max_jobs_param (actual concurrent jobs: $max_jobs)"
log_main "Threads per job: $threads (fixed)"
log_main "Mode: Continuous replacement (new job starts immediately when one finishes)"
log_main "Log folder: $LOG_SUBFOLDER"

# Find FASTQ files that need processing
all_samples=($(ls -f *.{fq,fastq}.gz 2>/dev/null))
samples_to_process=()

log_main "Scanning ${#all_samples[@]} FASTQ files..."

for sample in "${all_samples[@]}"; do
    if [[ ! -f "$sample" ]]; then
        log_error "File not found: $sample"
        continue
    fi
    
    sample_name=$(basename "$sample" .fq.gz)
    if [[ -d "$outFolder/fastQC/${sample_name}_fastqc" ]]; then
        log_main "SKIP: $sample (output exists)"
        continue
    fi
    
    samples_to_process+=("$sample")
done

log_main "Found ${#samples_to_process[@]} files to process (${#all_samples[@]} total, $((${#all_samples[@]} - ${#samples_to_process[@]})) skipped)"

if [[ ${#samples_to_process[@]} -eq 0 ]]; then
    log_main "No files need processing. All outputs already exist."
    exit 0
fi

# Job management arrays
declare -a job_pids=()
declare -a job_samples=()
declare -A sample_start_times=()

# Function to start a job
start_job() {
    local sample="$1"
    local sample_name=$(basename "$sample" .fq.gz)
    local job_log="$JOBS_LOG_DIR/${sample_name}_job.log"
    local job_error_log="$JOBS_LOG_DIR/${sample_name}_error.log"
    
    # Pre-create log files
    touch "$job_log"
    touch "$job_error_log"
    
    # Start job
    (
        start_time=$(date)
        
        # Initialize log file
        cat > "$job_log" << EOF
=== FastQC Analysis Job ===
Sample: $sample
Started: $start_time
PID: $$
Input file: $inFolder/$sample
Output directory: $outFolder/fastQC/
Threads: $threads
Working directory: $(pwd)
=================================
EOF
        
        # Execute FastQC command
        echo "Executing: fastqc --outdir \"$outFolder/fastQC/\" --format fastq --threads $threads \"$inFolder/$sample\"" >> "$job_log"
        echo "==================== FASTQC OUTPUT ====================" >> "$job_log"
        if fastqc --outdir "$outFolder/fastQC/" --format fastq --threads $threads "$inFolder/$sample" >> "$job_log" 2>> "$job_error_log"; then
            end_time=$(date)
            echo "==================== JOB COMPLETED ====================" >> "$job_log"
            echo "SUCCESS: $end_time" >> "$job_log"
            # Verify output directory was created
            if [[ -d "$outFolder/fastQC/${sample_name}_fastqc" ]]; then
                echo "Output directory verified: $outFolder/fastQC/${sample_name}_fastqc" >> "$job_log"
                exit 0
            else
                echo "ERROR: FastQC output directory not created despite success exit code" >> "$job_log"
                echo "FAILED: $end_time (no output directory)" >> "$job_log"
                exit 1
            fi
        else
            exit_code=$?
            end_time=$(date)
            echo "==================== JOB FAILED ====================" >> "$job_log"
            echo "FAILED: $end_time (exit code: $exit_code)" >> "$job_log"
            exit $exit_code
        fi
    ) &
    
    local job_pid=$!
    
    # Validate and add PID
    if [[ -n "$job_pid" ]] && kill -0 "$job_pid" 2>/dev/null; then
        job_pids+=("$job_pid")
        job_samples+=("$sample")
        sample_start_times["$sample"]=$(date)
        return 0
    else
        log_error "Failed to start job for: $sample"
        rm -f "$job_log" "$job_error_log"
        return 1
    fi
}

# Function to clean up completed jobs and return completed count
cleanup_completed_jobs() {
    local completed_count=0
    local new_pids=()
    local new_samples=()
    
    for i in "${!job_pids[@]}"; do
        local pid="${job_pids[i]}"
        local sample="${job_samples[i]}"
        
        if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
            # Job still running
            new_pids+=("$pid")
            new_samples+=("$sample")
        else
            # Job completed
            local start_time="${sample_start_times[$sample]}"
            log_main "COMPLETED: $sample [PID:$pid] (started: $start_time)"
            unset sample_start_times["$sample"]
            ((completed_count++))
        fi
    done
    
    job_pids=("${new_pids[@]}")
    job_samples=("${new_samples[@]}")
    
    return $completed_count
}

# Main processing loop with continuous job replacement
current_sample_index=0
completed_jobs=0
failed_jobs=0

log_main "Starting continuous job processing..."

while [[ $current_sample_index -lt ${#samples_to_process[@]} ]] || [[ ${#job_pids[@]} -gt 0 ]]; do
    # Clean up completed jobs
    cleanup_completed_jobs
    completed_this_cycle=$?
    
    if [[ $completed_this_cycle -gt 0 ]]; then
        completed_jobs=$((completed_jobs + completed_this_cycle))
        log_main "Progress: $completed_jobs completed, ${#job_pids[@]} active, $((${#samples_to_process[@]} - current_sample_index)) queued"
    fi
    
    # Start new jobs to fill available slots
    while [[ ${#job_pids[@]} -lt $max_jobs ]] && [[ $current_sample_index -lt ${#samples_to_process[@]} ]]; do
        sample="${samples_to_process[$current_sample_index]}"
        
        if start_job "$sample"; then
            log_main "STARTED: $sample [PID:${job_pids[-1]}] (${#job_pids[@]}/$max_jobs slots, job $((current_sample_index + 1))/${#samples_to_process[@]})"
        else
            log_error "Failed to start job for: $sample"
            ((failed_jobs++))
        fi
        
        ((current_sample_index++))
    done
    
    # Short sleep to prevent busy waiting
    sleep 2
    
    # Progress report every 5 minutes
    current_time=$(date +%s)
    if [[ -z "$last_progress_time" ]]; then
        last_progress_time=$current_time
    elif [[ $((current_time - last_progress_time)) -ge 300 ]]; then
        if [[ ${#job_pids[@]} -gt 0 ]]; then
            log_main "Status update: $completed_jobs completed, ${#job_pids[@]} active jobs, $((${#samples_to_process[@]} - current_sample_index)) remaining"
            log_main "Active jobs: ${job_samples[*]}"
        fi
        last_progress_time=$current_time
    fi
done

log_main "=== All Jobs Submitted and Completed ==="

# Generate final summary
success_count=0
fail_count=0

for sample in "${samples_to_process[@]}"; do
    sample_name=$(basename "$sample" .fq.gz)
    output_dir="$outFolder/fastQC/${sample_name}_fastqc"
    sample_name=$(basename "$sample" .fq.gz)
    job_log="$JOBS_LOG_DIR/${sample_name}_job.log"
    
    if [[ -d "$output_dir" ]]; then
        ((success_count++))
    elif [[ -f "$job_log" ]] && grep -q "SUCCESS:" "$job_log"; then
        ((success_count++))
    else
        ((fail_count++))
    fi
done

skipped_count=$((${#all_samples[@]} - ${#samples_to_process[@]}))

log_main "FINAL SUMMARY:"
log_main "  Processed successfully: $success_count"
log_main "  Failed: $fail_count"
log_main "  Skipped (output exists): $skipped_count"
log_main "  Total samples: ${#all_samples[@]}"
log_main ""
log_main "Efficiency: Continuous job replacement maintained $max_jobs concurrent jobs"
log_main "Log files: $LOG_SUBFOLDER"

if [[ $fail_count -gt 0 ]]; then
    log_error "Some jobs failed. Check individual job logs in: $JOBS_LOG_DIR"
fi

# Cleanup
rm -f "$MASTER_PID_FILE" 2>/dev/null

log_main "=== FastQC batch processing completed ==="

exit 0