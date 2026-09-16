#!/usr/bin/env python3
"""Write the reviewed species-reference manifest from verified downloads."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.reference_sequences import write_reviewed_manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--sequence-audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    try:
        summary = write_reviewed_manifest(
            args.inventory, args.sequence_audit, args.output
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe reviewed-manifest generation failed: {exc}\n")
        return 1
    print(
        f"references={summary['manifest_references']}\t"
        f"species={summary['manifest_species']}\t"
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
