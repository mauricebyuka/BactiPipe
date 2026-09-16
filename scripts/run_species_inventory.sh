#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 GTDB_METADATA CONTROLLED_CATALOG NCBI_TAXONOMY_DIR OUTPUT_DIRECTORY" >&2
    exit 2
fi

gtdb_metadata=$1
controlled_catalog=$2
ncbi_taxonomy_dir=$3
output_dir=$4
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/.." && pwd)
status_file="$output_dir/status.tsv"

mkdir -p "$output_dir"

record_status() {
    local state=$1
    local details=${2:-}
    printf '%s\tspecies_inventory\t%s\t%s\n' \
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

if [[ ! -s "$gtdb_metadata" || ! -s "$controlled_catalog" || ! -s "$ncbi_taxonomy_dir/merged.dmp" ]]; then
    echo "GTDB metadata, controlled catalog, and NCBI merged.dmp must exist and be non-empty" >&2
    exit 1
fi

record_status running
PYTHONPATH="$repo_root" python -m bactipipe.scripts.plan_species_inventory \
    --gtdb-metadata "$gtdb_metadata" \
    --controlled-catalog "$controlled_catalog" \
    --ncbi-taxonomy-dir "$ncbi_taxonomy_dir" \
    --output "$output_dir/reference_inventory.tsv" \
    --accessions-output "$output_dir/assembly_accessions.txt" \
    --summary-output "$output_dir/inventory_summary.json" \
    --missing-controlled-output "$output_dir/missing_controlled_species.tsv" \
    --controlled-aliases-output "$output_dir/controlled_taxid_aliases.tsv"

sha256sum \
    "$gtdb_metadata" \
    "$controlled_catalog" \
    "$ncbi_taxonomy_dir/merged.dmp" \
    "$output_dir/reference_inventory.tsv" \
    "$output_dir/assembly_accessions.txt" \
    "$output_dir/inventory_summary.json" \
    "$output_dir/missing_controlled_species.tsv" \
    "$output_dir/controlled_taxid_aliases.tsv" \
    > "$output_dir/SHA256SUMS"
record_status success "checksums=written"
