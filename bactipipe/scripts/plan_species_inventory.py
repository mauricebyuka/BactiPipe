#!/usr/bin/env python3
"""Create a bounded species-reference inventory from pinned GTDB metadata."""

from __future__ import annotations

import argparse
import sys

from bactipipe.reference_curation import (
    CurationPolicy,
    plan_gtdb_reference_inventory,
    write_inventory_plan,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Use a controlled veterinary species catalog and pinned GTDB metadata "
            "to propose a bounded set of target and same-genus near-neighbor assemblies."
        )
    )
    parser.add_argument("--gtdb-metadata", required=True)
    parser.add_argument("--controlled-catalog", required=True)
    parser.add_argument(
        "--ncbi-taxonomy-dir",
        required=True,
        help="Pinned NCBI taxdump directory containing merged.dmp",
    )
    parser.add_argument("--output", required=True, help="Auditable reference proposal TSV")
    parser.add_argument("--accessions-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--missing-controlled-output", required=True)
    parser.add_argument("--controlled-aliases-output", required=True)
    parser.add_argument("--minimum-completeness", type=float, default=90.0)
    parser.add_argument("--maximum-contamination", type=float, default=5.0)
    parser.add_argument("--maximum-contigs", type=int, default=500)
    parser.add_argument("--minimum-n50", type=int, default=20_000)
    parser.add_argument("--minimum-length", type=int, default=100_000)
    parser.add_argument("--maximum-diversity-references", type=int, default=3)
    args = parser.parse_args(argv)
    try:
        references, summary = plan_gtdb_reference_inventory(
            args.gtdb_metadata,
            args.controlled_catalog,
            policy=CurationPolicy(
                minimum_completeness=args.minimum_completeness,
                maximum_contamination=args.maximum_contamination,
                maximum_contigs=args.maximum_contigs,
                minimum_n50=args.minimum_n50,
                minimum_length=args.minimum_length,
            ),
            maximum_diversity_references=args.maximum_diversity_references,
            ncbi_taxonomy_dir=args.ncbi_taxonomy_dir,
        )
        write_inventory_plan(
            references,
            output_tsv=args.output,
            accessions_path=args.accessions_output,
            summary=summary,
            summary_path=args.summary_output,
            missing_controlled_path=args.missing_controlled_output,
            controlled_aliases_path=args.controlled_aliases_output,
        )
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe species inventory planning failed: {exc}\n")
        return 1
    counts = summary["counts"]
    print(
        f"selected_references={counts['selected_references']}\t"
        f"selected_species={counts['selected_species']}\t"
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
