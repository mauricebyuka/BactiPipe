#!/usr/bin/env python3
"""Build, but never activate, a reviewed species-identification database."""

from __future__ import annotations

import argparse
import sys

from bactipipe.species_database import build_identification_database


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build an immutable inactive BactiPipe species database candidate."
    )
    parser.add_argument("--reviewed-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument(
        "--ncbi-taxonomy-dir",
        required=True,
        help="Pinned NCBI taxonomy directory containing nodes.dmp and names.dmp",
    )
    parser.add_argument("--version")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--kraken-max-db-size",
        type=int,
        help=(
            "Optional Kraken2 hash-size ceiling in bytes. Use only after "
            "sensitivity validation."
        ),
    )
    parser.add_argument(
        "--kraken-library-batch-size",
        type=int,
        default=100,
        help="References per temporary Kraken import batch",
    )
    args = parser.parse_args(argv)
    try:
        output = build_identification_database(
            reviewed_manifest=args.reviewed_manifest,
            output_root=args.output_root,
            ncbi_taxonomy_dir=args.ncbi_taxonomy_dir,
            version=args.version,
            threads=args.threads,
            kraken_max_db_size=args.kraken_max_db_size,
            kraken_library_batch_size=args.kraken_library_batch_size,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe species database build failed: {exc}\n")
        return 1
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
