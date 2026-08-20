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

# Scale raw host coverage to a fixed retained dm6 observation target.  The
# optional ratio is the declared dm6:host amount added before tagmentation;
# multiplying by it prevents samples with deliberately different spike-in
# additions from being treated as if they received the same reference amount.
spikein_scale_factor() {
    local target="$1" spikein_to_host_ratio="$2" dm6_signal_count="$3"
    awk -v target="$target" -v ratio="$spikein_to_host_ratio" -v fly="$dm6_signal_count" '
        BEGIN {
            if (target <= 0 || ratio <= 0 || fly <= 0) exit 1
            printf "%.12g\n", target * ratio / fly
        }
    '
}

# Print five tab-separated fragment/read diagnostics using the same PE/SE
# observation definition as coverage: MAPQ 0, 1-9, 10-29, >=30, and XS-tagged.
signal_mapq_bin_counts() {
    local bam="$1" layout="${2^^}" total q1 q10 q30 xs
    local -a args=()
    if [[ "$layout" == "PE" ]]; then
        args=(-f 66 -F 4)
    elif [[ "$layout" == "SE" ]] && configured_se_signal_mode >/dev/null; then
        args=(-F 4)
    else
        return 1
    fi
    total="$(samtools view -c "${args[@]}" "$bam")" || return 1
    q1="$(samtools view -c -q 1 "${args[@]}" "$bam")" || return 1
    q10="$(samtools view -c -q 10 "${args[@]}" "$bam")" || return 1
    q30="$(samtools view -c -q 30 "${args[@]}" "$bam")" || return 1
    xs="$(samtools view "${args[@]}" "$bam" | awk '$0 ~ /\tXS:i:/ {n++} END{print n+0}')" || return 1
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$((total - q1))" "$((q1 - q10))" "$((q10 - q30))" "$q30" "$xs"
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
