#!/bin/bash
# Collect artifacts from sourcespec testruns
# Usage: ./collect_testrun_artifacts.sh [path_to_sourcespec_testruns] [destination_dir]
# Defaults: path="sourcespec_testruns", destination="testrun_artifacts"

set -euo pipefail

print_usage() {
    echo "Usage: $0 [path_to_sourcespec_testruns] [destination_dir]"
    echo "Defaults: path='sourcespec_testruns', destination='testrun_artifacts'"
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    print_usage
    exit 0
fi

TESTRUNS_PATH="${1:-sourcespec_testruns}"
DEST_DIR="${2:-testrun_artifacts}"

if [[ ! -d "$TESTRUNS_PATH" ]]; then
    echo "Directory $TESTRUNS_PATH not found!"
    echo "Please provide the correct path to sourcespec_testruns."
    exit 1
fi

# Resolve to absolute paths
TESTRUNS_PATH=$(cd "$TESTRUNS_PATH" && pwd)

REQUIRED_DIRS=(test_CDSA test_CRL test_IPOC test_ISNet)
for testdir in "${REQUIRED_DIRS[@]}"; do
    if [[ ! -d "$TESTRUNS_PATH/$testdir" ]]; then
        echo "Error: $testdir directory not found in $TESTRUNS_PATH"
        echo "This does not appear to be a valid sourcespec_testruns directory."
        exit 1
    fi
done

mkdir -p "$DEST_DIR"
# Resolve destination directory to an absolute path to avoid issues
# when archiving from subdirectories
DEST_DIR_ABS=$(cd "$DEST_DIR" && pwd)

created=0
for testdir in "${REQUIRED_DIRS[@]}"; do
    src_dir="$TESTRUNS_PATH/$testdir"
    sspec_dir="$src_dir/sspec_out"
    archive_path="$DEST_DIR_ABS/${testdir}_sspec_out.tar.gz"

    if [[ -d "$sspec_dir" ]]; then
        echo "Archiving sspec_out for $testdir -> $archive_path"
        # Exclude nested input_files directories (which may contain symlinks)
        # Works across GNU tar (Linux) and bsdtar (macOS/Windows Git Bash)
        tar --exclude='sspec_out/*/input_files' \
            --exclude='sspec_out/*/input_files/*' \
            -czf "$archive_path" -C "$src_dir" sspec_out
        created=$((created + 1))
    else
        echo "Warning: $sspec_dir not found; skipping $testdir"
    fi
done

echo "Created $created archive(s) in $DEST_DIR"
