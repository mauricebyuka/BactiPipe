#!/usr/bin/env bash

set -Eeuo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: $0 OUTPUT_DIRECTORY [SOURCE_URL] [EXPECTED_BYTES]" >&2
    exit 2
fi

output_dir=$1
source_url=${2:-https://download.genepi.dk/speciesfinder_db.tar.gz}
expected_bytes=${3:-54858324043}
download_dir="$output_dir/downloads"
archive="$download_dir/speciesfinder_db.tar.gz"
partial="$archive.part"
extracted="$output_dir/extracted"
extracting="$output_dir/.extracted.building"
status_file="$output_dir/status.tsv"

mkdir -p "$download_dir"
touch "$status_file"

record_status() {
    printf '%s\tspeciesfinder_db\t%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" | tee -a "$status_file"
}

on_error() {
    local exit_code=$?
    local line_number=$1
    trap - ERR
    record_status "failed; line=$line_number; exit_code=$exit_code; resumable_partial_retained=true"
    exit "$exit_code"
}
trap 'on_error "$LINENO"' ERR

record_status "download_started_or_resumed; url=$source_url; expected_bytes=$expected_bytes"
if [[ ! -f "$archive" ]]; then
    curl --fail --location --continue-at - \
        --retry 20 --retry-all-errors --retry-delay 30 \
        --connect-timeout 30 --speed-limit 1024 --speed-time 300 \
        --output "$partial" "$source_url"
    [[ "$(stat -c '%s' "$partial")" == "$expected_bytes" ]]
    mv "$partial" "$archive"
fi
[[ "$(stat -c '%s' "$archive")" == "$expected_bytes" ]]
record_status "download_size_verified; bytes=$expected_bytes"

sha256sum "$archive" > "$archive.sha256"
record_status "download_sha256_complete; sha256=$(cut -d' ' -f1 "$archive.sha256")"

tar -tzf "$archive" > "$download_dir/archive_inventory.txt"
grep -E '(^|/)bacteria/.*\.(comp\.b|seq\.b|length\.b|name|tax)$' \
    "$download_dir/archive_inventory.txt" >/dev/null
record_status "archive_inventory_verified"

if [[ ! -d "$extracted" ]]; then
    [[ ! -e "$extracting" ]]
    mkdir -p "$extracting"
    tar -xzf "$archive" -C "$extracting" --wildcards --no-anchored 'bacteria/*'
    mv "$extracting" "$extracted"
fi

database_prefix=$(find "$extracted" -type f -path '*/bacteria/*.name' -print -quit)
[[ -n "$database_prefix" ]]
database_prefix=${database_prefix%.name}
for suffix in .comp.b .seq.b .length.b .name; do
    [[ -s "${database_prefix}${suffix}" ]]
done
taxonomy_file=$(find "$(dirname "$database_prefix")" -maxdepth 1 -type f -name '*.tax' -print -quit)
[[ -n "$taxonomy_file" && -s "$taxonomy_file" ]]
printf '%s\n' "$database_prefix" > "$output_dir/bacteria_database_prefix.txt"
printf '%s\n' "$taxonomy_file" > "$output_dir/bacteria_taxonomy_file.txt"
record_status "ready_for_preparation; database_prefix=$database_prefix; taxonomy_file=$taxonomy_file"
