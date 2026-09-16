#!/usr/bin/env bash
set -euo pipefail
# Installed in genepid/bin; use the isolated KMA runtime in the sibling env.
kma_dispatch_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${kma_dispatch_dir}/../../bactipipe-kma/bin/kma" "$@"
