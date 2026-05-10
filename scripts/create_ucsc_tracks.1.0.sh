#!/bin/bash

# UCSC BigWig Track Generator Script
# Creates UCSC track lines for bigWig files in the specified folder

# Usage: ./create_ucsc_tracks.sh [output_folder] [remote_url_base]
# Example: ./create_ucsc_tracks.sh /dysk2/results https://myserver.com/data

# Function to display usage information
usage() {
    echo "Usage: $0  [remote_url_base]"
    echo ""
    echo "Parameters:"
    echo "  output_folder - Directory containing bigWig files (*.bw)"
    echo "  remote_url_base - Base URL where bigWig files will be accessible remotely"
    echo ""
    echo "Examples:"
    echo "  $0 /dysk2/results https://myserver.com/data"
    echo "  $0 /dysk2/results https://genome.ucsc.edu/goldenpath/help/examples"
    echo ""
    echo "Output:"
    echo "  Creates 'ucsc_tracks.txt' in the output folder with track definitions"
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
output_folder="$1"
remote_url_base="$2"

# Validate output folder exists
if [[ ! -d "$output_folder" ]]; then
    echo "ERROR: Output folder does not exist: $output_folder"
    exit 1
fi

# Check if bigwig subfolder exists
bigwig_folder="$output_folder/bigwig"
if [[ ! -d "$bigwig_folder" ]]; then
    echo "ERROR: BigWig folder does not exist: $bigwig_folder"
    echo "Expected folder structure: $output_folder/bigwig/"
    exit 1
fi

# Clean up remote URL base (remove trailing slash if present)
remote_url_base="${remote_url_base%/}"

# Output file for track definitions
track_file="$output_folder/ucsc_tracks.txt"

# Initialize the track file
echo "# UCSC BigWig Track Definitions" > "$track_file"
echo "# Generated on: $(date)" >> "$track_file"
echo "# Base URL: $remote_url_base" >> "$track_file"
echo "" >> "$track_file"

# Find all bigWig files and sort them
cd "$bigwig_folder" || {
    echo "ERROR: Cannot access bigwig folder: $bigwig_folder"
    exit 1
}

# Get list of bigWig files (*.bw and *.bigwig)
bigwig_files=($(ls -1 *.bw *.bigwig 2>/dev/null | sort))

# Check if any bigWig files were found
if [[ ${#bigwig_files[@]} -eq 0 ]]; then
    echo "WARNING: No bigWig files found in $bigwig_folder"
    echo "Looking for files with extensions: .bw, .bigwig"
    echo ""
    echo "Available files:"
    ls -la "$bigwig_folder"
    exit 1
fi

echo "Found ${#bigwig_files[@]} bigWig files in $bigwig_folder"

# Generate track definitions
n=1
for filename in "${bigwig_files[@]}"; do
    
    # Extract sample name from filename (remove extension)
    sample_name=$(basename "$filename" .bw)
    sample_name=$(basename "$sample_name" .bigwig)
    
    # Create track line according to the specified scheme
    track_line="track type=bigWig name=\"track_$n\" description=\"$filename\" bigDataUrl=$remote_url_base/$filename visibility=full color=0,0,0 priority=6 autoScale=off alwaysZero=on gridDefault=on graphType=bar windowingFunction=mean viewLimits=0:0.5 maxHeightPixels=100:50:8"
    
    # Write track line to file
    echo "$track_line" >> "$track_file"
    
    # Add blank line between tracks (except after the last track)
    if [[ $n -lt ${#bigwig_files[@]} ]]; then
        echo "" >> "$track_file"
    fi
    
    echo "Track $n: $filename -> track_$n"
    
    ((n++))
done

echo ""
echo "=== Track Generation Complete ==="
echo "Generated ${#bigwig_files[@]} track definitions"
echo "Output file: $track_file"
echo ""
echo "To use these tracks:"
echo "1. Upload your bigWig files to: $remote_url_base/"
echo "2. Copy the contents of $track_file"
echo "3. Paste into UCSC Genome Browser custom track input"
echo ""
echo "Track file preview:"
echo "===================="
head -15 "$track_file"

# Also create a summary file with file information
summary_file="$output_folder/bigwig_summary.txt"
echo "# BigWig Files Summary" > "$summary_file"
echo "# Generated on: $(date)" >> "$summary_file"
echo "# Total files: ${#bigwig_files[@]}" >> "$summary_file"
echo "" >> "$summary_file"

echo "Track_Number|Filename|Sample_Name|File_Size" >> "$summary_file"
echo "-------------|---------|-----------|----------" >> "$summary_file"

n=1
for filename in "${bigwig_files[@]}"; do
    sample_name=$(basename "$filename" .bw)
    sample_name=$(basename "$sample_name" .bigwig)
    
    # Get file size in human readable format
    file_size=$(du -h "$bigwig_folder/$filename" | cut -f1)
    
    echo "$n|$filename|$sample_name|$file_size" >> "$summary_file"
    ((n++))
done

echo ""
echo "Summary file created: $summary_file"

exit 0