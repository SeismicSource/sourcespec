#!/usr/bin/env bash
# Copyright (c) 2021-2026 Claudio Satriano <satriano@ipgp.fr>
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/cleanup_whitespace.sh [check|fix] [PATH ...]

Removes or reports trailing whitespace on any line.

Modes:
  check   Report offending lines and exit with code 1 if any are found
  fix     Clean offending lines in place (default)

PATH:
  Optional file or directory paths. Default is current directory.

Examples:
  scripts/cleanup_whitespace.sh check
  scripts/cleanup_whitespace.sh fix tests
  scripts/cleanup_whitespace.sh check requake tests docs
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

mode="fix"
case "${1:-}" in
    check|fix)
        mode="$1"
        shift
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    "")
        ;;
    *)
        die "Unknown mode '$1'. Use 'check' or 'fix'."
        ;;
esac

require_cmd rg
require_cmd perl

if [[ $# -eq 0 ]]; then
    targets=(.)
else
    targets=("$@")
fi

# Text-oriented files where whitespace cleanup is generally safe.
ext_pattern='\.(py|sh|md|rst|txt|toml|ya?ml|json|ini|cfg|conf|csv)$'

collect_files() {
    local root
    for root in "$@"; do
        if [[ -d "$root" ]]; then
            rg --files "$root" | rg "$ext_pattern" || true
        elif [[ -f "$root" ]]; then
            if [[ "$root" =~ $ext_pattern ]]; then
                echo "$root"
            fi
        fi
    done | sort -u
}

files=()
while IFS= read -r line; do
    files+=("$line")
done < <(collect_files "${targets[@]}")

if [[ ${#files[@]} -eq 0 ]]; then
    echo "No matching text files found."
    exit 0
fi

if [[ "$mode" == "check" ]]; then
    found=0
    for file in "${files[@]}"; do
        if rg -n "[ \t]+$" "$file" >/dev/null; then
            found=1
            echo "$file"
            rg -n "[ \t]+$" "$file"
            echo
        fi
    done
    if [[ $found -eq 1 ]]; then
        echo "Trailing whitespace found."
        exit 1
    fi
    echo "No trailing whitespace found."
    exit 0
fi

for file in "${files[@]}"; do
    perl -pi -e 's/[ \t]+$//' "$file"
done

echo "Cleanup complete."
