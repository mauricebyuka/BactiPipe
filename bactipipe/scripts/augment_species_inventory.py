#!/usr/bin/env python3
"""Add audited controlled-gap assemblies to a reconciled inventory."""

from __future__ import annotations

import argparse
import json
import sys

from bactipipe.reference_supplements import augment_inventory_with_gap_candidates


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--candidate-assessment", required=True)
    parser.add_argument("--controlled-catalog", required=True)
    parser.add_argument("--controlled-gaps", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args(argv)
    try:
        summary = augment_inventory_with_gap_candidates(
            args.inventory,
            args.candidate_assessment,
            args.controlled_catalog,
            args.controlled_gaps,
            args.output_root,
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"BactiPipe inventory augmentation failed: {exc}\n")
        return 1
    print(
        f"supplements={summary['supplement_reference_count']}\t"
        f"references={summary['final_reference_count']}\t"
        f"unresolved_gaps={summary['controlled_gap_count_after_recovery']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
