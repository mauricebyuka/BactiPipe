#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 OUTPUT_DIRECTORY" >&2
    exit 2
fi

output_dir=$1
base_url=https://data.gtdb.ecogenomic.org/releases/release232/232.0
status_file="$output_dir/status.tsv"

mkdir -p "$output_dir"

record_status() {
    local state=$1
    local details=${2:-}
    printf '%s\tgtdb_metadata\t%s\t%s\n' \
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

record_status running "release=R232"
curl -fL --retry 8 --retry-all-errors -C - \
    -o "$output_dir/MD5SUM.txt" "$base_url/MD5SUM.txt"
curl -fL --retry 8 --retry-all-errors -C - \
    -o "$output_dir/VERSION.txt" "$base_url/VERSION.txt"
curl -fL --retry 8 --retry-all-errors -C - \
    -o "$output_dir/bac120_metadata_r232.tsv.gz" \
    "$base_url/bac120_metadata_r232.tsv.gz"
curl -fL --retry 8 --retry-all-errors -C - \
    -o "$output_dir/bac120_taxonomy_r232.tsv.gz" \
    "$base_url/bac120_taxonomy_r232.tsv.gz"

(
    cd "$output_dir"
    rg 'bac120_(metadata|taxonomy)_r232.tsv.gz$' MD5SUM.txt | md5sum -c -
)

record_status success "release=R232; checksums=verified"
