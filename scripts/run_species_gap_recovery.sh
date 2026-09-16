#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 CONTROLLED_GAPS_TSV DATASETS_EXECUTABLE OUTPUT_DIRECTORY" >&2
    exit 2
fi

controlled_gaps=$1
datasets_executable=$2
output_dir=$3
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
status_file="$output_dir/status.tsv"
reports="$output_dir/assembly_data_report.jsonl"

mkdir -p "$output_dir"

record_status() {
    local state=$1
    local details=${2:-}
    printf '%s\tcontrolled_gap_recovery\t%s\t%s\n' \
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

if [[ ! -s "$controlled_gaps" || ! -x "$datasets_executable" ]]; then
    echo "Controlled-gap TSV and executable NCBI Datasets client are required" >&2
    exit 1
fi

record_status running
awk -F '\t' 'NR > 1 && $2 ~ /^[0-9]+$/ {print $2}' "$controlled_gaps" \
    | sort -u > "$output_dir/taxids.txt.tmp"
mv "$output_dir/taxids.txt.tmp" "$output_dir/taxids.txt"
if [[ ! -s "$output_dir/taxids.txt" ]]; then
    echo "No controlled-gap taxids were found" >&2
    exit 1
fi

if [[ ! -s "$reports" ]]; then
    "$datasets_executable" summary genome taxon \
        --inputfile "$output_dir/taxids.txt" \
        --tax-exact-match \
        --assembly-source all \
        --assembly-version current \
        --mag exclude \
        --exclude-atypical \
        --exclude-multi-isolate \
        --as-json-lines > "$reports.tmp"
    if [[ ! -s "$reports.tmp" ]]; then
        echo "NCBI returned no assembly reports for the controlled gaps" >&2
        exit 1
    fi
    mv "$reports.tmp" "$reports"
fi

PYTHONPATH="$repo_root" python -m bactipipe.scripts.plan_species_references \
    --ncbi-reports "$reports" \
    --output "$output_dir/candidate_assessment.tsv"

awk -F '\t' 'NR == 1 || $17 == "selected_canonical"' \
    "$output_dir/candidate_assessment.tsv" > "$output_dir/selected_canonical.tsv.tmp"
mv "$output_dir/selected_canonical.tsv.tmp" "$output_dir/selected_canonical.tsv"

sha256sum \
    "$controlled_gaps" \
    "$output_dir/taxids.txt" \
    "$reports" \
    "$output_dir/candidate_assessment.tsv" \
    "$output_dir/selected_canonical.tsv" \
    > "$output_dir/SHA256SUMS"
record_status success "candidate_assessment=candidate_assessment.tsv"
