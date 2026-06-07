#!/bin/bash
# ATACseq2tracks v3.0.4 — MACS2 batch (per-replicate + pooled; tech-replicate aware)
# Usage: bash scripts/macs2_batch.sh
set -euo pipefail
_load_config() {
  if [[ -n "${F2T_CONFIG:-}" && -f "${F2T_CONFIG}" ]]; then
    source "${F2T_CONFIG}"
  else
    local _d; _d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    local _c="${_d}/../config/config.conf"
    [[ -f "$_c" ]] && source "$_c" || {
      echo "ERROR: config.sh not found. Export F2T_CONFIG or pass --config to atacseq2tracks.sh." >&2
      exit 1
    }
  fi
}
_load_config

SAMPLESHEET="$1"; BAM_DIR="$2"; PEAKS_DIR="$3"
mkdir -p "$PEAKS_DIR"
LOG="${PEAKS_DIR}/macs2_batch_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"; }
declare -A ip_bams ip_ctrls ip_modes ip_genomes
declare -A grp_bams grp_ctrl grp_mode grp_genome
declare -A seen_keys=()
while IFS=',' read -r sid fq1 fq2 layout genome assay factor condition treatment cell_type rep tech_rep is_ctr ctrl_id macs2_mode blacklist rest; do
  [[ "$sid" == "sample_id" ]] && continue
  sid="${sid//\"/}"; is_ctr="${is_ctr//\"/}"; ctrl_id="${ctrl_id//\"/}"
  genome="${genome//\"/}"; macs2_mode="${macs2_mode//\"/}"
  factor="${factor//\"/}"; condition="${condition//\"/}"
  treatment="${treatment//\"/}"; cell_type="${cell_type//\"/}"; rep="${rep//\"/}"
  [[ "${is_ctr,,}" == "true" || "$is_ctr" == "1" ]] && continue
  KEY="${sid}_bioR${rep}"
  [[ -n "${seen_keys[$KEY]+x}" ]] && continue
  seen_keys["$KEY"]=1
  BAM="${BAM_DIR}/${KEY}_dedup_blFilt.bam"
  [[ ! -f "$BAM" ]] && BAM="${BAM_DIR}/${KEY}_dedup.bam"
  [[ ! -f "$BAM" ]] && { log "WARN: BAM not found for $KEY"; continue; }
  CTRL_BAM="none"
  if [[ -n "$ctrl_id" ]]; then
    CB="${BAM_DIR}/${ctrl_id}_bioR${rep}_dedup_blFilt.bam"
    [[ ! -f "$CB" ]] && CB="${BAM_DIR}/${ctrl_id}_bioR1_dedup_blFilt.bam"
    [[ ! -f "$CB" ]] && CB="${BAM_DIR}/${ctrl_id}_dedup_blFilt.bam"
    [[ ! -f "$CB" ]] && CB="${BAM_DIR}/${ctrl_id}_dedup.bam"
    [[ -f "$CB" ]] && CTRL_BAM="$CB"
  fi
  ip_bams["$KEY"]="$BAM"; ip_ctrls["$KEY"]="$CTRL_BAM"
  ip_modes["$KEY"]="$macs2_mode"; ip_genomes["$KEY"]="$genome"
  GRP="${factor}__${condition}__${treatment}__${cell_type}__${genome}"
  grp_bams["$GRP"]+="$BAM "; grp_ctrl["$GRP"]="$CTRL_BAM"
  grp_mode["$GRP"]="$macs2_mode"; grp_genome["$GRP"]="$genome"
done < "$SAMPLESHEET"
log "=== Per-replicate MACS2 ==="
for KEY in "${!ip_bams[@]}"; do
  [[ "${ip_modes[$KEY],,}" == "none" ]] && { log "SKIP $KEY (macs2_mode=none)"; continue; }
  OUT="${PEAKS_DIR}/per_replicate/${KEY}"; mkdir -p "$OUT"
  (bash "$(dirname "$0")/macs2_peaks.sh" \
    "${ip_bams[$KEY]}" "${ip_ctrls[$KEY]}" "$OUT" \
    "both" "${ip_genomes[$KEY]}" "$KEY" >> "$LOG" 2>&1) &
done
wait
log "=== Pooled-replicate MACS2 ==="
for GRP in "${!grp_bams[@]}"; do
  BAMS=(${grp_bams[$GRP]})
  [[ ${#BAMS[@]} -eq 0 ]] && continue
  [[ "${grp_mode[$GRP],,}" == "none" ]] && continue
  OUT="${PEAKS_DIR}/pooled/${GRP}"; mkdir -p "$OUT"
  GRP_NAME="${GRP//__/_}_pooled"
  if [[ ${#BAMS[@]} -eq 1 ]]; then MERGED_BAM="${BAMS[0]}"
  else
    MERGED_BAM="${PEAKS_DIR}/pooled/${GRP_NAME}_tmp.bam"
    samtools merge -@ "${THREADS_SAMTOOLS}" -f "$MERGED_BAM" "${BAMS[@]}"
    samtools index -@ "${THREADS_SAMTOOLS}" "$MERGED_BAM"
  fi
  (bash "$(dirname "$0")/macs2_peaks.sh" \
    "$MERGED_BAM" "${grp_ctrl[$GRP]}" "$OUT" \
    "both" "${grp_genome[$GRP]}" "$GRP_NAME" >> "$LOG" 2>&1) &
done
wait
# Verify every expected peak file was produced
n_warn=0
for sid in "${!ip_bams[@]}"; do
    [[ "${ip_modes[$sid],,}" == "none" ]] && continue
    OUT="${PEAKS_DIR}/per_replicate/${sid}"
    if [[ ! -f "${OUT}/narrow/${sid}_peaks.narrowPeak" && ! -f "${OUT}/broad/${sid}_peaks.broadPeak" ]]; then
        log "ERROR: no peak file produced for ${sid}"
        (( n_warn++ )) || true
    fi
done
[[ $n_warn -gt 0 ]] && { log "FATAL: $n_warn sample(s) produced no peaks — aborting step"; exit 1; }
log "=== MACS2 batch complete ==="
