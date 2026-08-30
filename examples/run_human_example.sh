#!/usr/bin/env bash
# Example only. Edit paths before running.
set -euo pipefail

CONFIG="/data/project/human/config/config.conf"
LOG="/data/project/human/atacseq2tracks.log"

atacseq2tracks --config "$CONFIG" --preflight-only
nohup atacseq2tracks --config "$CONFIG" > "$LOG" 2>&1 &
echo "$!" > /data/project/human/atacseq2tracks.pid
