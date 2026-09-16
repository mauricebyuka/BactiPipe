#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 ACCESSIONS DATASETS_EXECUTABLE OUTPUT_DIRECTORY" >&2
    exit 2
fi

accessions=$1
datasets_executable=$2
output_dir=$3
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
status_file="$output_dir/status.tsv"

mkdir -p "$output_dir"

record_status() {
    local state=$1
    local details=${2:-}
    printf '%s\tncbi_assembly_reports\t%s\t%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$state" "$details" \
        > "$status_file"
}

on_error() {
    local exit_code=$?
    trap - ERR
    record_status failed "exit_code=$exit_code"
    exit "$exit_code"
}
trap on_error ERR

record_status running
PYTHONPATH="$repo_root" python -m bactipipe.scripts.fetch_species_reports \
    --accessions "$accessions" \
    --datasets-executable "$datasets_executable" \
    --output-root "$output_dir"
record_status success "fetch_summary=fetch_summary.json"
