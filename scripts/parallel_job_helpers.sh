#!/usr/bin/env bash
# Shared bounded-process-pool helpers for ATACseq2tracks batch stages.
# The caller remains responsible for giving each worker unique output files.

parallel_require_positive_integer() {
    local name="${1:?setting name required}" value="${2:-}"
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: $name must be a positive integer; found: ${value:-<empty>}" >&2
        return 1
    }
}

parallel_pool_init() {
    local maximum="${1:?maximum jobs required}"
    parallel_require_positive_integer "parallel job limit" "$maximum" || return 1
    PARALLEL_POOL_MAX="$maximum"
    PARALLEL_POOL_FAILURES=0
    PARALLEL_POOL_PIDS=()
    PARALLEL_POOL_LABELS=()
    PARALLEL_POOL_FAILED_LABELS=()
}

parallel_pool_wait_one() {
    local pid label
    (( ${#PARALLEL_POOL_PIDS[@]} > 0 )) || return 0
    pid="${PARALLEL_POOL_PIDS[0]}"
    label="${PARALLEL_POOL_LABELS[0]}"
    if ! wait "$pid"; then
        PARALLEL_POOL_FAILURES=$((PARALLEL_POOL_FAILURES + 1))
        PARALLEL_POOL_FAILED_LABELS+=("$label")
    fi
    PARALLEL_POOL_PIDS=("${PARALLEL_POOL_PIDS[@]:1}")
    PARALLEL_POOL_LABELS=("${PARALLEL_POOL_LABELS[@]:1}")
}

parallel_pool_submit() {
    local label="${1:?job label required}"
    shift
    (( $# > 0 )) || { echo "ERROR: no worker command supplied for $label" >&2; return 1; }
    while (( ${#PARALLEL_POOL_PIDS[@]} >= PARALLEL_POOL_MAX )); do
        parallel_pool_wait_one
    done
    "$@" &
    PARALLEL_POOL_PIDS+=("$!")
    PARALLEL_POOL_LABELS+=("$label")
}

parallel_pool_wait_all() {
    while (( ${#PARALLEL_POOL_PIDS[@]} > 0 )); do
        parallel_pool_wait_one
    done
    (( PARALLEL_POOL_FAILURES == 0 ))
}

parallel_failed_labels_csv() {
    local IFS=,
    printf '%s' "${PARALLEL_POOL_FAILED_LABELS[*]:-}"
}

parallel_write_timing_row() {
    local file="${1:?timing file required}" scope="${2:?scope required}"
    local label="${3:?label required}" start="${4:?start epoch required}"
    local end="${5:?end epoch required}" status="${6:?status required}"
    local jobs="${7:-NA}" threads="${8:-NA}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$scope" "$label" "$start" "$end" "$((end - start))" "$jobs" "$threads" "$status" \
        > "$file"
}
