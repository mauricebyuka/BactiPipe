#!/usr/bin/env python3
"""Reconcile proposed species references against current NCBI reports."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.reference_inventory import reconcile_reference_inventory


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Resolve current accession versions and promote available canonicals."
    )
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--accession-resolution", required=True)
    parser.add_argument("--ncbi-reports", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--pre-inventory-missing")
    parser.add_argument(
        "--ncbi-taxonomy-dir",
        required=True,
        help="Pinned NCBI taxdump directory containing merged.dmp",
    )
    args = parser.parse_args(argv)
    try:
        summary = reconcile_reference_inventory(
            args.inventory,
            args.accession_resolution,
            args.ncbi_reports,
            args.output_root,
            pre_inventory_missing_tsv=args.pre_inventory_missing,
            ncbi_taxonomy_dir=args.ncbi_taxonomy_dir,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe inventory reconciliation failed: {exc}\n")
        return 1
    print(
        f"current_references={summary['current_reference_count']}\t"
        f"current_species={summary['current_species_count']}\t"
        f"missing_species={summary['missing_species_count']}\t"
        f"taxonomy_review={summary['taxonomy_review_reference_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
