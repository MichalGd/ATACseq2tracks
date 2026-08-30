#!/usr/bin/env bash
set -euo pipefail

INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG=""
OUTPUT_OVERRIDE=""
usage() {
    cat <<'EOF'
Usage: regenerate_reports.sh --config /absolute/path/config.conf [--output-dir DIR]
Regenerates reports only; no alignment, peak calling, tracks or differential models run.
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="${2:?missing value for --config}"; shift 2 ;;
        --output-dir) OUTPUT_OVERRIDE="${2:?missing value for --output-dir}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
tmp_dir="$(mktemp -d)"
trap 'rm -rf -- "$tmp_dir"' EXIT
python3 "${INSTALL_DIR}/scripts/resolve_config.py" \
    --config "$CONFIG" --template "${INSTALL_DIR}/config/config.conf" \
    --shell-output "$tmp_dir/resolved.conf" --tsv-output "$tmp_dir/resolved.tsv"
# shellcheck disable=SC1090
source "$tmp_dir/resolved.conf"
[[ -n "$OUTPUT_OVERRIDE" ]] && OUTPUT_DIR="$OUTPUT_OVERRIDE"
[[ -d "$OUTPUT_DIR" ]] || { echo "ERROR: output directory not found: $OUTPUT_DIR" >&2; exit 1; }
export F2T_CONFIG="$tmp_dir/resolved.conf"
bash "${INSTALL_DIR}/scripts/generate_pipeline_report.sh" "$OUTPUT_DIR" \
    "${OUTPUT_DIR}/reports/pipeline_report_$(date +%Y%m%d)" html
