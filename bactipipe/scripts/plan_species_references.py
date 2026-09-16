#!/usr/bin/env python3
"""Create an auditable canonical-reference selection proposal."""

from __future__ import annotations

import argparse
import sys

from bactipipe.reference_curation import (
    CurationPolicy,
    plan_references,
    read_ncbi_jsonl,
    write_curation_plan,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Rank NCBI Datasets assembly reports and propose one canonical "
            "isolate-derived reference per resolved bacterial species."
        )
    )
    parser.add_argument("--ncbi-reports", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--minimum-completeness", type=float, default=90.0)
    parser.add_argument("--maximum-contamination", type=float, default=5.0)
    parser.add_argument("--maximum-contigs", type=int, default=500)
    parser.add_argument("--minimum-n50", type=int, default=20_000)
    parser.add_argument("--minimum-length", type=int, default=100_000)
    args = parser.parse_args(argv)
    try:
        rows = plan_references(
            read_ncbi_jsonl(args.ncbi_reports),
            policy=CurationPolicy(
                minimum_completeness=args.minimum_completeness,
                maximum_contamination=args.maximum_contamination,
                maximum_contigs=args.maximum_contigs,
                minimum_n50=args.minimum_n50,
                minimum_length=args.minimum_length,
            ),
        )
        write_curation_plan(rows, args.output)
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe species-reference planning failed: {exc}\n")
        return 1
    selected = sum(row.decision == "selected_canonical" for row in rows)
    print(f"records={len(rows)}\tselected_species={selected}\toutput={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
