#!/bin/bash
# Test runner for sourcespec
# This script runs test examples from the sourcespec_testruns repository
#
# Usage: ./do_testruns.sh [path_to_sourcespec_testruns]
# If no path is provided, defaults to "sourcespec_testruns"
# in the current directory

# Exit on error, undefined variable, or pipe failure
set -euo pipefail

# Print usage information
print_usage() {
    echo "Usage: $0 [path_to_sourcespec_testruns]"
    echo "Defaults: path='sourcespec_testruns' in the current directory"
}

# Show usage on --help or -h
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    print_usage
    exit 0
fi

# Use provided path or default to sourcespec_testruns
TESTRUNS_PATH="${1:-sourcespec_testruns}"

cd "$TESTRUNS_PATH" 2> /dev/null || {
    echo "Directory $TESTRUNS_PATH not found!"
    echo "Please clone the repository or provide the correct path."
    exit 1
}

# Verify that this is a valid testruns directory
REQUIRED_DIRS="test_CDSA test_CRL test_IPOC test_ISNet"
for testdir in $REQUIRED_DIRS; do
    if [ ! -d "$testdir" ]; then
        echo "Error: $testdir directory not found in $TESTRUNS_PATH"
        echo "This does not appear to be a valid sourcespec_testruns directory."
        exit 1
    fi
done

for testdir in $REQUIRED_DIRS; do
    echo ""
    echo "=========================================="
    echo "Running $testdir..."
    echo "=========================================="
    cd "$testdir"
    echo "y" | source_spec -U source_spec.conf
    ./run.sh
    cd ..
    echo ""
    echo "=========================================="
    echo "✓ $testdir completed"
    echo "=========================================="
    echo ""
done

echo "All testruns passed!"