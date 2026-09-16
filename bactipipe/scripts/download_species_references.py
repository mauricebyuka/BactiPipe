#!/usr/bin/env python3
"""Download reconciled bacterial reference assemblies with checkpoints."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.reference_sequences import (
    download_reference_genomes,
    write_reviewed_manifest,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Download reconciled NCBI assemblies and write a reviewed manifest."
    )
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--datasets-executable", default="datasets")
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    args = parser.parse_args(argv)
    try:
        summary = download_reference_genomes(
            args.inventory,
            args.output_root,
            datasets_executable=args.datasets_executable,
            batch_size=args.batch_size,
            retries=args.retries,
            timeout_seconds=args.timeout_seconds,
        )
        manifest = write_reviewed_manifest(
            args.inventory,
            f"{args.output_root}/sequence_download_audit.tsv",
            f"{args.output_root}/reviewed_reference_manifest.tsv",
        )
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe reference download failed: {exc}\n")
        return 1
    print(
        f"downloaded={summary['downloaded_count']}\t"
        f"missing={summary['missing_count']}\t"
        f"manifest_references={manifest['manifest_references']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
