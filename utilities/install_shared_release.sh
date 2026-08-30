#!/usr/bin/env bash
set -euo pipefail

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RELEASE_ROOT="${ATACSEQ2TRACKS_RELEASE_ROOT:-/opt/bioinformatics/workflows/ATACseq2tracks/releases}"
STABLE_LINK="${ATACSEQ2TRACKS_STABLE_LINK:-/opt/bioinformatics/workflows/ATACseq2tracks/current}"
COMMAND_LINK="${ATACSEQ2TRACKS_COMMAND_LINK:-/usr/local/bin/atacseq2tracks}"
MAIN_ENV="${ATACSEQ2TRACKS_MAIN_ENV:-/opt/miniconda/envs/ATACseq2tracks}"
ATAQV_ENV="${ATACSEQ2TRACKS_ATAQV_ENV:-/opt/miniconda/envs/ataqv-tools}"
VERSION="$(tr -d '[:space:]' < "${SOURCE_DIR}/VERSION")"
DESTINATION="${RELEASE_ROOT}/${VERSION}"
STAGING="${RELEASE_ROOT}/.${VERSION}.stage.$$"

[[ $EUID -eq 0 ]] || { echo "ERROR: run with sudo" >&2; exit 1; }
[[ ! -e "$DESTINATION" ]] || { echo "ERROR: release already exists: $DESTINATION" >&2; exit 1; }
[[ ! -e "$STAGING" ]] || { echo "ERROR: staging path already exists: $STAGING" >&2; exit 1; }
for tool in python3 Rscript bowtie2 samtools bedtools trim_galore fastqc macs3 multiqc bamCoverage; do
    [[ -x "${MAIN_ENV}/bin/${tool}" ]] || { echo "ERROR: missing managed tool: ${MAIN_ENV}/bin/${tool}" >&2; exit 1; }
done
for tool in ataqv mkarv; do
    [[ -x "${ATAQV_ENV}/bin/${tool}" ]] || { echo "ERROR: missing ataqv-sidecar tool: ${ATAQV_ENV}/bin/${tool}" >&2; exit 1; }
done
unset PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER
export PATH="${MAIN_ENV}/bin:${ATAQV_ENV}/bin:/usr/local/bin:/usr/bin:/bin"
install -d -m 0755 "$RELEASE_ROOT"
cp -a "$SOURCE_DIR" "$STAGING"
find "$STAGING" -type d -exec chmod 0755 {} +
find "$STAGING" -type f -exec chmod 0644 {} +
find "$STAGING" -type f \( -name '*.sh' -o -name '*.py' -o -path '*/bin/*' \) -exec chmod 0755 {} +
chown -R root:root "$STAGING"

(cd "$STAGING" && bash tests/check_bash_syntax.sh)
mv "$STAGING" "$DESTINATION"
ln -sfn "$DESTINATION" "${STABLE_LINK}.${VERSION}.new"
mv -Tf "${STABLE_LINK}.${VERSION}.new" "$STABLE_LINK"
ln -sfn "${STABLE_LINK}/bin/atacseq2tracks" "${COMMAND_LINK}.${VERSION}.new"
mv -Tf "${COMMAND_LINK}.${VERSION}.new" "$COMMAND_LINK"
echo "Installed ATACseq2tracks ${VERSION}: $DESTINATION"
echo "Active release: $(readlink -f "$STABLE_LINK")"
echo "Command: $COMMAND_LINK"
