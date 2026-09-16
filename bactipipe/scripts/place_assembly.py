"""Explicit post-QC placement; does not regenerate or change QC reports."""
import argparse
import json
from bactipipe.iterative_phylogeny import place_iteratively


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--cache-root", required=True, help="Private persistent cache for this sample")
    parser.add_argument("--reference-root", required=True, help="Immutable GTDB-Tk data directory")
    parser.add_argument("--metadata", required=True, help="Matching GTDB species metadata TSV[.gz]")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument('--identify-first', action='store_true', help='Confirm whole-assembly ANI against GTDB before considering a tree')
    parser.add_argument('--identification-json', help='Original species_identification.json, including thresholds and competitors')
    parser.add_argument('--sample', default='query')
    parser.add_argument('--expected-organism', default='Unknown')
    args = parser.parse_args()
    options = dict(assembly=args.assembly, cache_root=args.cache_root,
        reference_root=args.reference_root, metadata=args.metadata, threads=args.threads)
    if args.identify_first:
        from pathlib import Path
        from bactipipe.species_confirmation import investigate_assembly
        evidence = json.loads(Path(args.identification_json).read_text()) if args.identification_json else {}
        result = investigate_assembly(**options, evidence=evidence, sample=args.sample, expected=args.expected_organism)
    else:
        result = place_iteratively(**options)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
