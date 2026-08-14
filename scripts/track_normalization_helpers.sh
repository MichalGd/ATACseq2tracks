#!/usr/bin/env bash
# Shared fragment/read and canonical-chromosome semantics for coverage normalization.

canonical_autosome_max() {
    case "${1,,}" in
        hg38) printf '22\n' ;;
        mm39) printf '19\n' ;;
        *) return 1 ;;
    esac
}

canonical_contigs_from_bam() {
    local bam="$1" genome="$2" maximum
    maximum="$(canonical_autosome_max "$genome")" || return 1
    samtools idxstats "$bam" | awk -v maximum="$maximum" '
        $1 != "*" {
            original=$1
            name=$1
            sub(/^chr/, "", name)
            if ((name ~ /^[0-9]+$/ && name+0 >= 1 && name+0 <= maximum) || name == "X" || name == "Y") {
                print original
            }
        }
    '
}

canonicalize_peak_file() {
    local input="$1" output="$2" genome="$3" maximum
    maximum="$(canonical_autosome_max "$genome")" || return 1
    awk -v maximum="$maximum" 'BEGIN{OFS="\t"}
        NF >= 3 {
            name=$1
            sub(/^chr/, "", name)
            if ((name ~ /^[0-9]+$/ && name+0 >= 1 && name+0 <= maximum) || name == "X" || name == "Y") {
                print $1,$2,$3
            }
        }
    ' "$input" > "$output"
}

resolve_signal_layout() {
    local bam="$1" requested="${2:-}"
    requested="${requested^^}"
    if [[ "$requested" == "PE" || "$requested" == "SE" ]]; then
        printf '%s\n' "$requested"
        return 0
    fi
    if (( $(samtools view -c -f 1 -F 4 "$bam" 2>/dev/null || echo 0) > 0 )); then
        printf 'PE\n'
    else
        printf 'SE\n'
    fi
}

configured_se_signal_mode() {
    local mode="${SE_SIGNAL_MODE:-read}"
    mode="${mode,,}"
    case "$mode" in
        read) printf 'read\n' ;;
        *) return 1 ;;
    esac
}

deeptools_signal_args() {
    local layout="${1^^}" destination="$2"
    # Bash nameref keeps the arguments safely separated when passed to deepTools.
    local -n result="$destination"
    result=()
    if [[ "$layout" == "PE" ]]; then
        # The filtered BAM already contains proper pairs. Requiring 0x2 + 0x40
        # selects one first-mate record per fragment, and --extendReads spans
        # the paired-end insert defined by TLEN.
        result=(--extendReads --samFlagInclude 66)
    elif [[ "$layout" == "SE" ]]; then
        # A single-end alignment does not identify a physical fragment. In
        # the supported read mode, one retained alignment is one observation
        # and no artificial extension is applied.
        configured_se_signal_mode >/dev/null || return 1
    else
        return 1
    fi
}

signal_unit_for_layout() {
    if [[ "${1^^}" == "PE" ]]; then
        printf 'fragment\n'
    elif [[ "${1^^}" == "SE" ]] && configured_se_signal_mode >/dev/null; then
        printf 'read\n'
    else
        return 1
    fi
}

signal_count_for_bam() {
    local bam="$1" layout="${2^^}"
    if [[ "$layout" == "PE" ]]; then
        samtools view -c -f 66 -F 4 "$bam"
    elif [[ "$layout" == "SE" ]] && configured_se_signal_mode >/dev/null; then
        samtools view -c -F 4 "$bam"
    else
        return 1
    fi
}

write_scaled_coverage_track() {
    local bam="$1" output="$2" format="$3" scale="$4" layout="$5"
    local bin_size="${6:-10}" threads="${7:-2}"
    local -a signal_args=()
    deeptools_signal_args "$layout" signal_args || return 1
    bamCoverage --bam "$bam" --outFileName "$output" --outFileFormat "$format" \
        --normalizeUsing None --scaleFactor "$scale" --binSize "$bin_size" \
        --numberOfProcessors "$threads" "${signal_args[@]}" \
        ${BAMCOVERAGE_COMMON_ARGS:-} || return 1
    [[ -s "$output" ]]
}
