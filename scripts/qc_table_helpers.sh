#!/usr/bin/env bash
# Strict table parsers used by the v3.2.0 test-run hotfix.

picard_duplication_pct() {
    local metrics_file="$1"
    [[ -s "$metrics_file" ]] || return 1
    awk -F '\t' '
        /^LIBRARY(\t|$)/ {
            for (i = 1; i <= NF; i++) {
                header = $i
                gsub(/\r/, "", header)
                if (header == "PERCENT_DUPLICATION") column = i
            }
            next
        }
        column && NF > 1 {
            value = $column
            gsub(/^[[:space:]]+|[[:space:]\r]+$/, "", value)
            if (value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value >= 0 && value <= 1) {
                printf "%.2f\n", value * 100
                found = 1
                exit
            }
        }
        END { if (!found) exit 1 }
    ' "$metrics_file"
}

consensus_table_value_for_key() {
    local table_file="$1" wanted_key="$2" wanted_column="$3"
    [[ -s "$table_file" ]] || return 1
    awk -F '\t' -v wanted="$wanted_key" -v wanted_column="$wanted_column" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                header = $i
                gsub(/\r/, "", header)
                if (header == "key") key_column = i
                if (header == wanted_column) value_column = i
            }
            if (!key_column || !value_column) exit 2
            next
        }
        {
            key = $key_column
            gsub(/\r/, "", key)
            if (key == wanted) {
                value = $value_column
                gsub(/^[[:space:]]+|[[:space:]\r]+$/, "", value)
                if (value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value > 0) {
                    print value
                    found = 1
                    exit
                }
                exit 3
            }
        }
        END { if (!found) exit 1 }
    ' "$table_file"
}

consensus_size_factor_for_key() {
    consensus_table_value_for_key "$1" "$2" size_factor
}
