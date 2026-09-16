#!/usr/bin/env bash
set -euo pipefail

environment_name="${1:-genepid}"

# These exact legacy wrappers declare obsolete dependency pins (notably
# Biopython 1.73), although their assembly-mode code works with the compatible
# modern dependencies declared in genepid.yml. Install the wrappers without
# re-solving those stale metadata constraints.
conda install -y -n "$environment_name" --no-deps \
  https://conda.anaconda.org/bioconda/noarch/mlst-cge-2.0.9-hdfd78af_0.tar.bz2 \
  https://conda.anaconda.org/bioconda/noarch/seqsero2-1.3.2-pyhdfd78af_0.conda \
  https://conda.anaconda.org/bioconda/noarch/serotypefinder-2.0.1-py39hdfd78af_0.tar.bz2

genepid_prefix="$(conda run -n "$environment_name" python -c 'import sys; print(sys.prefix)')"
genepid_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
test -x "${genepid_prefix}/../bactipipe-kma/bin/kma"
install -m 0755 "${genepid_script_dir}/kma_dispatch.sh" "${genepid_prefix}/bin/kma"
