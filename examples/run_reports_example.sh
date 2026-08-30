#!/usr/bin/env bash
# Example only. Edit paths before running.
set -euo pipefail

CONFIG="/data/project/human/config/config.conf"

/opt/bioinformatics/workflows/ATACseq2tracks/current/utilities/regenerate_reports.sh \
    --config "$CONFIG"
