#!/usr/bin/env bash
# Download, validate and install the native-GRCh38 ENCODE4 expanded cCRE registry.
# No liftOver is required: the source file is already on GRCh38/hg38.
set -euo pipefail

SOURCE_URL="https://users.moore-lab.org/ENCODE-cCREs/Supplementary-Data/Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
SOURCE_NAME="Supplementary-Data-1.GRCh38-cCREs-V4.bed.gz"
EXPECTED_RECORDS=2348854
DEFAULT_OUTPUT_DIR="/home/micgdu/Analysis/utilities/UCSC/CREs/human/hg38"

usage() {
    cat <<EOF
Usage:
  bash utilities/prepare_encode4_hg38_ccre.sh [--output-dir DIR] [--force]
  bash utilities/prepare_encode4_hg38_ccre.sh --validate-only FILE

The default installed file is:
  ${DEFAULT_OUTPUT_DIR}/${SOURCE_NAME}

Set CCRE_BED_HG38 to that file. Use --output-dir on another server.
EOF
}

sha256_file() {
    local path="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$path" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$path" | awk '{print $1}'
    else
        echo "ERROR: sha256sum or shasum is required" >&2
        return 1
    fi
}

validate_bed() {
    local path="$1" expected_records="${2:-0}"
    [[ -s "$path" ]] || { echo "ERROR: cCRE BED is missing or empty: $path" >&2; return 1; }
    gzip -t "$path" || { echo "ERROR: invalid gzip stream: $path" >&2; return 1; }
    gzip -dc "$path" | awk -F '\t' -v expected="$expected_records" '
        BEGIN {
            valid_class["PLS"]=1; valid_class["pELS"]=1; valid_class["dELS"]=1
            valid_class["CA-CTCF"]=1; valid_class["CA-H3K4me3"]=1
            valid_class["CA-TF"]=1; valid_class["CA"]=1; valid_class["TF"]=1
            valid_class["CTCF-bound"]=1; valid_class["CTCF-only"]=1
            valid_class["DNase-H3K4me3"]=1
        }
        /^#/ || /^[[:space:]]*$/ { next }
        {
            row++
            canonical = ($1 == "chrX" || $1 == "chrY" ||
                         $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/)
            if (!canonical) {
                print "ERROR: noncanonical chromosome at record " row ": " $1 > "/dev/stderr"
                exit 20
            }
            if (NF < 6 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 ||
                $4 == "" || $5 == "" || $6 == "") {
                print "ERROR: invalid six-column cCRE BED record " row > "/dev/stderr"
                exit 21
            }
            count = split($6, classes, ",")
            for (class_index = 1; class_index <= count; class_index++) {
                if (!(classes[class_index] in valid_class)) {
                    print "ERROR: unsupported cCRE class at record " row ": " classes[class_index] > "/dev/stderr"
                    exit 22
                }
            }
        }
        END {
            if (row == 0) {
                print "ERROR: no cCRE records found" > "/dev/stderr"
                exit 23
            }
            if (expected > 0 && row != expected) {
                print "ERROR: expected " expected " cCRE records but found " row > "/dev/stderr"
                exit 24
            }
            print row
        }
    '
}

OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
FORCE=false
VALIDATE_ONLY=""
while (( $# > 0 )); do
    case "$1" in
        --output-dir)
            (( $# >= 2 )) || { echo "ERROR: --output-dir requires a value" >&2; exit 2; }
            OUTPUT_DIR="$2"; shift 2 ;;
        --force) FORCE=true; shift ;;
        --validate-only)
            (( $# >= 2 )) || { echo "ERROR: --validate-only requires a file" >&2; exit 2; }
            VALIDATE_ONLY="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -n "$VALIDATE_ONLY" ]]; then
    records="$(validate_bed "$VALIDATE_ONLY" 0)"
    echo "OK: compatible hg38 cCRE BED ($records records): $VALIDATE_ONLY"
    exit 0
fi

for tool in gzip awk date; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: required tool not found: $tool" >&2; exit 1; }
done
if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: curl or wget is required" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
TARGET="${OUTPUT_DIR}/${SOURCE_NAME}"
PROVENANCE="${OUTPUT_DIR}/${SOURCE_NAME%.bed.gz}.provenance.tsv"
if [[ -e "$TARGET" && "$FORCE" != "true" ]]; then
    records="$(validate_bed "$TARGET" "$EXPECTED_RECORDS")"
    echo "OK: existing compatible ENCODE4 hg38 cCRE reference ($records records): $TARGET"
    echo "Use --force to download and replace it."
    exit 0
fi

TMP_DIR="$(mktemp -d "${OUTPUT_DIR}/.encode4_hg38_ccre.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT
DOWNLOAD="${TMP_DIR}/${SOURCE_NAME}"
if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 2 --output "$DOWNLOAD" "$SOURCE_URL"
else
    wget --tries=3 --output-document="$DOWNLOAD" "$SOURCE_URL"
fi

records="$(validate_bed "$DOWNLOAD" "$EXPECTED_RECORDS")"
source_sha256="$(sha256_file "$DOWNLOAD")"
chmod 0644 "$DOWNLOAD"
mv -f "$DOWNLOAD" "$TARGET"
installed_sha256="$(sha256_file "$TARGET")"
[[ "$source_sha256" == "$installed_sha256" ]] || {
    echo "ERROR: checksum changed while installing $TARGET" >&2
    exit 1
}

{
    printf 'field\tvalue\n'
    printf 'source_url\t%s\n' "$SOURCE_URL"
    printf 'source_filename\t%s\n' "$SOURCE_NAME"
    printf 'downloaded_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'assembly\tGRCh38/hg38\n'
    printf 'registry\tENCODE4 expanded Registry of cCREs\n'
    printf 'records\t%s\n' "$records"
    printf 'bed_columns\tchrom,chromStart,chromEnd,accession1,accession2,class\n'
    printf 'sha256\t%s\n' "$installed_sha256"
} > "$PROVENANCE"
chmod 0644 "$PROVENANCE"

echo "OK: installed ENCODE4 hg38 cCRE reference"
echo "BED: $TARGET"
echo "Records: $records"
echo "SHA256: $installed_sha256"
echo "Provenance: $PROVENANCE"
echo "Set: CCRE_BED_HG38=\"$TARGET\""
