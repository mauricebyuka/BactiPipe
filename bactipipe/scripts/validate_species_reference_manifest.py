#!/usr/bin/env python3
"""Validate every reference, checksum, role, and species in a candidate manifest."""

from __future__ import annotations

import argparse
import sys

from bactipipe.species_identification import ReferenceCatalog


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args(argv)
    try:
        catalog = ReferenceCatalog.from_tsv(args.manifest)
        catalog.validate_files(require_checksums=True)
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe manifest validation failed: {exc}\n")
        return 1
    print(
        f"references={len(catalog.records)}\t"
        f"species={len({record.species_taxid for record in catalog.records})}\t"
        "status=valid"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
