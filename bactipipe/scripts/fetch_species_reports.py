#!/usr/bin/env python3
"""Fetch current NCBI assembly reports for a proposed reference inventory."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.ncbi_reference_metadata import fetch_assembly_reports


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Fetch current NCBI assembly metadata in resumable accession batches."
    )
    parser.add_argument("--accessions", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--datasets-executable", default="datasets")
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        summary = fetch_assembly_reports(
            args.accessions,
            args.output_root,
            datasets_executable=args.datasets_executable,
            batch_size=args.batch_size,
            retries=args.retries,
        )
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe NCBI report retrieval failed: {exc}\n")
        return 1
    print(
        f"requested={summary['requested_count']}\t"
        f"retrieved={summary['retrieved_count']}\t"
        f"missing={summary['missing_count']}\t"
        f"status={summary['status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
