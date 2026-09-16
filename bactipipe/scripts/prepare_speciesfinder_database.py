#!/usr/bin/env python3
"""Prepare an immutable SpeciesFinder candidate database from a verified archive."""

from __future__ import annotations

import argparse
import sys

from bactipipe.speciesfinder_database import prepare_speciesfinder_database


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the official SpeciesFinder bacterial KMA index for bounded "
            "candidate nomination and exact-template ANI validation."
        )
    )
    parser.add_argument("--database-prefix", required=True)
    parser.add_argument("--taxonomy-file", required=True)
    parser.add_argument("--archive-sha256-file", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--version")
    parser.add_argument(
        "--production",
        action="store_true",
        help="Mark active only after the scientific validation gate has passed.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = prepare_speciesfinder_database(
            database_prefix=args.database_prefix,
            taxonomy_file=args.taxonomy_file,
            archive_sha256_file=args.archive_sha256_file,
            output_root=args.output_root,
            version=args.version,
            status="production" if args.production else "candidate",
        )
    except (OSError, RuntimeError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe SpeciesFinder database preparation failed: {exc}\n")
        return 1
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
