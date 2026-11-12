#!/usr/bin/env bash
set -euo pipefail
# Build the validated correction in a separate build stage. Never mutate a DB.
kma_source="${1:?Usage: install_kma_reassign_fix.sh PINNED_SOURCE INSTALL_PREFIX}"
kma_prefix="${2:?Missing installation prefix}"
kma_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
test "$(git -C "$kma_source" rev-parse HEAD)" = 5cfd66b310437035d2f866d97dfba8b85dc7e408
test "$(sha256sum "$kma_source/reassign.c" | cut -d ' ' -f 1)" = 515547a1ffb49d2db9dd6f8f9afa611ea5e41e6e591aa65f18e90c1a6ec6454d
patch --batch --fuzz=0 -d "$kma_source" -p1 -i "$kma_script_dir/patches/kma-1.6.13-reassign-progress.patch"
patch --batch --fuzz=0 -d "$kma_source" -p1 -i "$kma_script_dir/patches/kma-1.6.13-build-label.patch"
make -C "$kma_source" -j2 CFLAGS='-O2 -g -Wall -std=c99' kma
cc -O2 -std=c99 -I "$kma_source" -o "$kma_source/reassign_probe" \
    "$kma_script_dir/patches/kma_reassign_probe.c" "$kma_source/libkma.a" -lm -lpthread -lz
python3 - "$kma_source/reassign_probe" <<'PY'
import subprocess
import sys
for consensus, candidate in [('A' * 1290, 'C' * 168), ('A' * 400, 'A' * 96),
                             ('C' * 90 + 'A' * 168, 'A' * 168)]:
    result = int(subprocess.check_output([sys.argv[1], consensus, candidate], text=True, timeout=3))
    assert result == consensus.find(candidate), result
PY
install -D -m 0755 "$kma_source/kma" "$kma_prefix/bin/kma"
test "$("$kma_prefix/bin/kma" -v 2>&1)" = KMA-1.6.13+adl-reassign-progress-1
python3 - "$kma_prefix" "$kma_script_dir" <<'PY'
import hashlib
import json
import sys
from pathlib import Path
prefix, scripts = map(Path, sys.argv[1:])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
record = {
    'upstream_version': '1.6.13',
    'upstream_commit': '5cfd66b310437035d2f866d97dfba8b85dc7e408',
    'patch_id': 'adl-reassign-progress-1',
    'version': '1.6.13+adl-reassign-progress-1',
    'algorithm_patch_sha256': digest(scripts / 'patches/kma-1.6.13-reassign-progress.patch'),
    'label_patch_sha256': digest(scripts / 'patches/kma-1.6.13-build-label.patch'),
    'binary_sha256': digest(prefix / 'bin/kma'),
    'reassignment_enabled': True,
}
output = prefix / 'share/bactipipe/kma-build.json'
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(record, indent=2) + '\n')
PY
