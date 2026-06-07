#!/usr/bin/env bash
# Example only. Edit paths before running.
set -euo pipefail

RESULTS=/data/project/ATACseq2tracks_human

../scripts/generate_pipeline_report.1.0.sh   "$RESULTS"   "pipeline_report_$(date +%Y%m%d)"   html

../scripts/generate_multiqc_unified_report.1.0.sh   "$RESULTS"   "$RESULTS/reports/multiqc_summary_$(date +%Y%m%d)"   selfcontained
