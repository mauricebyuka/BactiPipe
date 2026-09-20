#!/usr/bin/env python3
"""Run the evidence-based BactiPipe species-identification stage."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bactipipe.species_identification import IdentificationPolicy, MixedSignalPolicy
from bactipipe.species_identification import (
    default_kraken_confidence,
    default_mixed_signal_policy,
)
from bactipipe.species_workflow import identification_database_kind, identify_assembly


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bactipipe identify",
        description=(
            "Identify an assembled bacterial genome using whole-genome ANI and "
            "aligned fraction, while failing closed on mixed or ambiguous evidence."
        ),
        epilog="Separate placement: bactipipe identify --place --help. Redraw saved figures: bactipipe identify --render --help.",
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--assembly", required=True, help="Draft assembly FASTA")
    parser.add_argument(
        "--database-root",
        required=True,
        help="Built BactiPipe species-identification database version",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--expected-organism",
        default="Unknown",
        help="Expected organism; Unknown is accepted",
    )
    parser.add_argument(
        "--kmerfinder-results",
        help=(
            "Optional existing legacy KmerFinder results.txt retained only for "
            "side-by-side validation; it never changes the identification"
        ),
    )
    evidence = parser.add_mutually_exclusive_group()
    evidence.add_argument(
        "--taxonomic-signals",
        help=(
            "Optional classifier-independent TSV containing species-level read count, "
            "abundance, distinct minimizers, and reference breadth"
        ),
    )
    parser.add_argument(
        "--reads",
        action="append",
        default=[],
        help="Filtered FASTQ; provide once for Nanopore or twice for paired Illumina",
    )
    parser.add_argument("--platform", choices=["Illumina", "Nanopore"], required=True)
    parser.add_argument(
        "--kraken-confidence",
        type=float,
        help="Classifier candidate threshold; defaults to 0.10 for Illumina and 0.05 for Nanopore",
    )
    parser.add_argument(
        "--gtdb-metadata",
        help=(
            "GTDB metadata TSV; when supplied, unresolved non-mixed samples receive "
            "a bac120 phylogenetic placement"
        ),
    )
    parser.add_argument(
        "--gtdbtk-conda-env",
        default="",
        help="Optional Conda environment containing GTDB-Tk",
    )
    parser.add_argument(
        "--gtdbtk-in-memory",
        action="store_true",
        help=(
            "Run pplacer without a disk-backed scratch directory. This requires "
            "substantial RAM and is intended for controlled placement runs."
        ),
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--minimum-ani", type=float, default=95.0)
    parser.add_argument("--minimum-aligned-fraction", type=float, default=50.0,
                        help="Minimum AF percent (default: 50; either query or reference)")
    parser.add_argument("--aligned-fraction-mode", choices=("either", "both"), default="either",
                        help="AF aggregation: either=max (GTDB-Tk); both=min (legacy)")
    parser.add_argument("--minimum-ani-margin", type=float, default=0.5)
    parser.add_argument("--minimum-aligned-fraction-margin", type=float, default=2.0)
    parser.add_argument("--minimum-secondary-reads", type=int)
    parser.add_argument("--minimum-secondary-distinct-minimizers", type=int)
    parser.add_argument("--minimum-secondary-read-percent", type=float)
    parser.add_argument("--minimum-secondary-reference-breadth", type=float)
    parser.add_argument("--comparable-secondary-read-ratio", type=float)
    parser.add_argument(
        "--maximum-candidate-species",
        type=int,
        default=20,
        help="Maximum SpeciesFinder taxa retained for ANI validation (default: 20)",
    )
    parser.add_argument(
        "--candidate-references-per-species",
        type=int,
        default=3,
        help="Maximum exact SpeciesFinder templates retained per taxon (default: 3)",
    )
    parser.add_argument(
        "--allow-candidate-database",
        action="store_true",
        help="Permit an inactive candidate database for controlled benchmarking",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = Path(args.outdir).resolve()
        command_prefix = (
            ("conda", "run", "-n", args.gtdbtk_conda_env, "gtdbtk")
            if args.gtdbtk_conda_env
            else ("gtdbtk",)
        )
        platform = args.platform or ""
        defaults = default_mixed_signal_policy(platform)
        kraken_confidence = args.kraken_confidence
        if kraken_confidence is None:
            kraken_confidence = default_kraken_confidence(platform)
        result = identify_assembly(
            sample_id=args.sample,
            assembly=args.assembly,
            database_root=args.database_root,
            output_directory=output,
            expected_organism=args.expected_organism,
            reads=args.reads,
            platform=platform,
            kraken_database=(
                Path(args.database_root).resolve() / "kraken2"
                if not args.taxonomic_signals
                and identification_database_kind(args.database_root)
                != "speciesfinder_kma"
                else None
            ),
            taxonomic_signals_tsv=args.taxonomic_signals,
            kmerfinder_results=args.kmerfinder_results,
            gtdb_metadata_tsv=args.gtdb_metadata,
            gtdbtk_command_prefix=command_prefix,
            gtdbtk_use_scratch=not args.gtdbtk_in_memory,
            threads=args.threads,
            kraken_confidence=kraken_confidence,
            policy=IdentificationPolicy(
                minimum_ani=args.minimum_ani,
                minimum_aligned_fraction=args.minimum_aligned_fraction,
                aligned_fraction_mode=args.aligned_fraction_mode,
                minimum_ani_margin=args.minimum_ani_margin,
                minimum_aligned_fraction_margin=(
                    args.minimum_aligned_fraction_margin
                ),
            ),
            mixed_signal_policy=MixedSignalPolicy(
                minimum_secondary_reads=(
                    args.minimum_secondary_reads
                    if args.minimum_secondary_reads is not None
                    else defaults.minimum_secondary_reads
                ),
                minimum_secondary_distinct_minimizers=(
                    args.minimum_secondary_distinct_minimizers
                    if args.minimum_secondary_distinct_minimizers is not None
                    else defaults.minimum_secondary_distinct_minimizers
                ),
                minimum_secondary_read_percent=(
                    args.minimum_secondary_read_percent
                    if args.minimum_secondary_read_percent is not None
                    else defaults.minimum_secondary_read_percent
                ),
                minimum_secondary_reference_breadth=(
                    args.minimum_secondary_reference_breadth
                    if args.minimum_secondary_reference_breadth is not None
                    else defaults.minimum_secondary_reference_breadth
                ),
                comparable_to_primary_ratio=(
                    args.comparable_secondary_read_ratio
                    if args.comparable_secondary_read_ratio is not None
                    else defaults.comparable_to_primary_ratio
                ),
            ),
            allow_candidate_database=args.allow_candidate_database,
            require_mixed_signal_assessment=True,
            maximum_candidate_species=args.maximum_candidate_species,
            candidate_references_per_species=args.candidate_references_per_species,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        sys.stderr.write(f"BactiPipe species identification failed: {exc}\n")
        return 1
    print(
        f"{result.sample_id}\t{result.status}\t{result.reported_organism}\t"
        f"{result.match_percent if result.match_percent is not None else 'N/A'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
