#!/bin/bash
# =============================================================================
# fastq2tracks v3.0 — Pipeline report wrapper
# Usage: bash scripts/generate_pipeline_report.sh <outDir> <reportName> <format>
# Wraps the existing generate_pipeline_report.1.0.sh logic (now in scripts/)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
exec "$(dirname "$0")/generate_multiqc_unified_report.sh" "$@"
